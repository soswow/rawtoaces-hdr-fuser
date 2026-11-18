#include "pipeline.h"

#include <cmath>
#include <limits>

// Motion estimation between frames (docs/ALIGNMENT.md, stage 1). Frames are
// compared as log-rate quad images, so a gain between them is an offset; a
// translation and that offset are fitted per tile by robust Gauss-Newton.

namespace rawfuse {
namespace {

constexpr int tileQuads = 64;        // Tile edge in quads (128 sites)
constexpr int coverageStride = 8;    // Sites between coverage samples
constexpr double snrCap = 100;       // Per-quad SNR limit for the matching weight
constexpr double modelError = 0.01;  // Relative error the noise model does not cover
constexpr double cauchy = 2.385;     // Cauchy constant for unit-scale residuals
constexpr double acceptStd = 0.1;    // Largest tile translation uncertainty (quads)
constexpr double acceptJump = 4;     // Largest tile departure from the global fit (quads)
constexpr double minimumSnr = 10;    // Quads darker than this carry too little to match
constexpr double levelError = 2;     // Additive error per frame (black level, dark signal), DN
constexpr double fieldTolerance = 0.1;  // Camera-motion departure from the rigid model, sites
constexpr double outlierFloor = 0.5;    // Tiles further than this from it are scene motion, sites
constexpr size_t minimumInliers = 16;   // Fewest consistent tiles for a usable estimate

// A quad is used for matching when all four sites are usable (not clipped,
// defective or outside the group's validity) and its mean is well above noise.
// Besides random noise, the variance holds an additive error shared by the
// quad: a frame's black level or dark signal can be off by a DN or two. In log
// it becomes an offset that grows toward dark quads; across a brightness
// gradient (a sky) it reads as a shift, so dark quads must weigh less.
bool matchable(const RawFrame& f, int qy, int qx, const Options& o, double& mean, double& variance) {
    double sum = 0;
    variance = 0;
    for (int k = 0; k < 4; ++k) {
        const int y = 2 * qy + k / 2, x = 2 * qx + k % 2;
        const auto s = sample(f, y, x, f.layout.index(y, x), o);
        if (s.state != SampleState::usable) return false;
        const auto& noise = f.noise[f.layout.colour(y, x)];
        sum += s.fraction;
        variance += noise.shot * std::max(0.0f, s.fraction) + noise.read + s.extra;
    }
    mean = sum / 4;
    const double level = levelError / (f.white - f.black(2 * qy, 2 * qx));
    variance = variance / 16 + level * level;
    return mean > 0 && mean * mean >= minimumSnr * minimumSnr * variance;
}

struct Level {
    int rows = 0, cols = 0;
    std::vector<float> value, weight;  // Log rate; inverse variance of it (0: undefined)
    size_t at(int r, int c) const { return size_t(r) * cols + c; }
};
using Pyramid = std::vector<Level>;

Level quadLevel(const RawFrame& f, const Options& o) {
    Level level;
    level.rows = f.layout.height / 2;
    level.cols = f.layout.width / 2;
    level.value.assign(size_t(level.rows) * level.cols, 0);
    level.weight.assign(level.value.size(), 0);
    parallelRows(level.rows, o.threads, [&](int qy) {
        for (int qx = 0; qx < level.cols; ++qx) {
            double mean, variance;
            if (!matchable(f, qy, qx, o, mean, variance)) continue;
            // E[log m] = log μ − σ²/2μ² to second order. The bias differs between
            // exposures (their SNRs differ), so it is removed before matching.
            const double noise = variance / (mean * mean);
            const double relative = noise + 1 / (snrCap * snrCap);
            level.value[level.at(qy, qx)] = float(std::log(mean / f.effective) + noise / 2);
            level.weight[level.at(qy, qx)] = float(1 / relative);
        }
    });
    return level;
}

Level reduce(const Level& in, unsigned threads) {
    Level out;
    out.rows = in.rows / 2;
    out.cols = in.cols / 2;
    out.value.assign(size_t(out.rows) * out.cols, 0);
    out.weight.assign(out.value.size(), 0);
    parallelRows(out.rows, threads, [&](int r) {
        for (int c = 0; c < out.cols; ++c) {
            double sum = 0, total = 0;
            int count = 0;
            for (int k = 0; k < 4; ++k) {
                const size_t i = in.at(2 * r + k / 2, 2 * c + k % 2);
                if (in.weight[i] <= 0) continue;
                sum += double(in.weight[i]) * in.value[i];
                total += in.weight[i];
                ++count;
            }
            if (count < 3) continue;
            out.value[out.at(r, c)] = float(sum / total);
            out.weight[out.at(r, c)] = float(total);
        }
    });
    return out;
}

Pyramid buildPyramid(const RawFrame& f, const Options& o) {
    Pyramid pyramid{quadLevel(f, o)};
    while (pyramid.size() < 7 && std::min(pyramid.back().rows, pyramid.back().cols) >= 80)
        pyramid.push_back(reduce(pyramid.back(), o.threads));
    return pyramid;
}

// Bilinear value, inverse variance and exact gradient of the interpolant.
bool interpolate(const Level& l, double y, double x, double& value, double& weight, double& gy, double& gx) {
    const int r = int(std::floor(y)), c = int(std::floor(x));
    if (r < 0 || c < 0 || r + 1 >= l.rows || c + 1 >= l.cols) return false;
    const size_t i00 = l.at(r, c), i01 = i00 + 1, i10 = i00 + l.cols, i11 = i10 + 1;
    const float w00 = l.weight[i00], w01 = l.weight[i01], w10 = l.weight[i10], w11 = l.weight[i11];
    if (w00 <= 0 || w01 <= 0 || w10 <= 0 || w11 <= 0) return false;
    const double fy = y - r, fx = x - c;
    const double v00 = l.value[i00], v01 = l.value[i01], v10 = l.value[i10], v11 = l.value[i11];
    const double top = v00 + fx * (v01 - v00), bottom = v10 + fx * (v11 - v10);
    value = top + fy * (bottom - top);
    gy = bottom - top;
    gx = (1 - fy) * (v01 - v00) + fy * (v11 - v10);
    const double variance = (1 - fy) * ((1 - fx) / w00 + fx / w01) + fy * ((1 - fx) / w10 + fx / w11);
    weight = 1 / variance;
    return true;
}

struct Params { double dy = 0, dx = 0, bias = 0; };
struct Fit {
    Params p;
    bool converged = false;
    size_t points = 0;
    double scale = 1;      // Robust spread of the normalized residuals
    double std = INFINITY; // Larger standard deviation of the translation (level units)
};

struct Normal {
    double h[3][3]{}, g[3]{}, n = 0;
    void add(const Normal& o) {
        for (int i = 0; i < 3; ++i) {
            g[i] += o.g[i];
            for (int j = 0; j < 3; ++j) h[i][j] += o.h[i][j];
        }
        n += o.n;
    }
};

bool solve3(const double h[3][3], const double g[3], double x[3]) {
    const double det = h[0][0] * (h[1][1] * h[2][2] - h[1][2] * h[2][1]) -
                       h[0][1] * (h[1][0] * h[2][2] - h[1][2] * h[2][0]) +
                       h[0][2] * (h[1][0] * h[2][1] - h[1][1] * h[2][0]);
    if (!(std::abs(det) > 1e-300)) return false;
    for (int k = 0; k < 3; ++k) {
        double m[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) m[i][j] = j == k ? g[i] : h[i][j];
        x[k] = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])) / det;
    }
    return true;
}

// Translation block of the inverse normal matrix: its larger eigenvalue.
double translationVariance(const double h[3][3]) {
    double inv[3][3];
    for (int k = 0; k < 3; ++k) {
        const double e[3] = {double(k == 0), double(k == 1), double(k == 2)};
        double x[3];
        if (!solve3(h, e, x)) return INFINITY;
        for (int i = 0; i < 3; ++i) inv[i][k] = x[i];
    }
    const double a = inv[0][0], b = inv[0][1], d = inv[1][1];
    return 0.5 * (a + d) + std::sqrt(0.25 * (a - d) * (a - d) + b * b);
}

// Catmull-Rom (Keys, a = -1/2) weights and their derivatives at phase t.
void cubicWeights(double t, double w[4], double d[4]) {
    const double t2 = t * t, t3 = t2 * t;
    w[0] = 0.5 * (-t3 + 2 * t2 - t);
    w[1] = 0.5 * (3 * t3 - 5 * t2 + 2);
    w[2] = 0.5 * (-3 * t3 + 4 * t2 + t);
    w[3] = 0.5 * (t3 - t2);
    d[0] = 0.5 * (-3 * t2 + 4 * t - 1);
    d[1] = 0.5 * (9 * t2 - 10 * t);
    d[2] = 0.5 * (-9 * t2 + 8 * t + 1);
    d[3] = 0.5 * (3 * t2 - 2 * t);
}

// Cubic value, inverse variance and gradient. Bilinear interpolation pulls a
// matched shift toward whole pixels by a few percent of its fractional part; the
// cubic kernel reduces that bias by an order of magnitude. Next to undefined
// quads, where the 4 × 4 support is incomplete, the bilinear value is used.
bool interpolateCubic(const Level& l, double y, double x, double& value, double& weight, double& gy, double& gx) {
    const int r = int(std::floor(y)), c = int(std::floor(x));
    if (r < 1 || c < 1 || r + 2 >= l.rows || c + 2 >= l.cols)
        return interpolate(l, y, x, value, weight, gy, gx);
    double wy[4], dy[4], wx[4], dx[4];
    cubicWeights(y - r, wy, dy);
    cubicWeights(x - c, wx, dx);
    double v = 0, vy = 0, vx = 0, variance = 0;
    for (int i = 0; i < 4; ++i) {
        const size_t row = l.at(r - 1 + i, c - 1);
        for (int j = 0; j < 4; ++j) {
            const float w = l.weight[row + j];
            if (w <= 0) return interpolate(l, y, x, value, weight, gy, gx);
            const double s = l.value[row + j], k = wy[i] * wx[j];
            v += k * s;
            vy += dy[i] * wx[j] * s;
            vx += wy[i] * dx[j] * s;
            variance += k * k / w;
        }
    }
    value = v;
    gy = vy;
    gx = vx;
    weight = 1 / variance;
    return true;
}

// Robust Gauss-Newton for (dy, dx, bias) on reference rows [r0, r1) and
// columns [c0, c1): mov(x + d) = ref(x) + bias.
Fit gaussNewton(const Level& ref, const Level& mov, int r0, int r1, int c0, int c1, int stride,
                Params p, int iterations, unsigned threads, bool cubic) {
    auto sampleAt = cubic ? interpolateCubic : interpolate;
    Fit fit;
    double scale = INFINITY;  // First pass without robust weights
    const int rows = (r1 - r0 + stride - 1) / stride;
    std::vector<Normal> partial(size_t(std::max(rows, 0)));
    std::vector<std::vector<float>> residuals(partial.size());
    for (int it = 0; it < iterations; ++it) {
        parallelRows(rows, threads, [&](int k) {
            const int r = r0 + k * stride;
            Normal n;
            auto& z = residuals[size_t(k)];
            z.clear();
            for (int c = c0; c < c1; c += stride) {
                const size_t i = ref.at(r, c);
                if (ref.weight[i] <= 0) continue;
                double value, weight, gy, gx;
                if (!sampleAt(mov, r + p.dy, c + p.dx, value, weight, gy, gx)) continue;
                const double res = value - ref.value[i] - p.bias;
                const double variance = 1 / double(ref.weight[i]) + 1 / weight + modelError * modelError;
                const double normalized = res / std::sqrt(variance);
                z.push_back(float(std::abs(normalized)));
                const double u = std::isfinite(scale) ? normalized / (cauchy * scale) : 0;
                const double w = 1 / ((1 + u * u) * variance);
                const double j[3] = {gy, gx, -1};
                for (int a = 0; a < 3; ++a) {
                    n.g[a] -= w * j[a] * res;
                    for (int b = 0; b < 3; ++b) n.h[a][b] += w * j[a] * j[b];
                }
                n.n += 1;
            }
            partial[size_t(k)] = n;
        });
        Normal total;
        std::vector<float> all;
        for (size_t k = 0; k < partial.size(); ++k) {
            total.add(partial[k]);
            all.insert(all.end(), residuals[k].begin(), residuals[k].end());
        }
        fit.points = size_t(total.n);
        if (all.size() < 16) return fit;
        auto middle = all.begin() + all.size() / 2;
        std::nth_element(all.begin(), middle, all.end());
        scale = std::max(1.4826 * double(*middle), 1e-3);
        double step[3];
        if (!solve3(total.h, total.g, step)) return fit;
        // Gauss-Newton steps beyond a level pixel leave the linear range.
        const double length = std::hypot(step[0], step[1]);
        const double limit = length > 1 ? 1 / length : 1;
        p.dy += limit * step[0];
        p.dx += limit * step[1];
        p.bias += limit * step[2];
        fit.p = p;
        fit.scale = scale;
        fit.std = std::sqrt(translationVariance(total.h)) * scale;
        if (it > 0 && length < 1e-3) {
            fit.converged = true;
            break;
        }
    }
    return fit;
}

// Exhaustive integer search at the coarsest level, robust to large shifts.
Params coarseSearch(const Level& ref, const Level& mov, int radius) {
    Params best;
    double bestCost = INFINITY;
    size_t mostPoints = 0;
    struct Candidate { int dy, dx; double cost, bias; size_t points; };
    std::vector<Candidate> candidates;
    for (int dy = -radius; dy <= radius; ++dy)
        for (int dx = -radius; dx <= radius; ++dx) {
            double sw = 0, swr = 0;
            std::vector<std::pair<double, double>> values;
            for (int r = 0; r < ref.rows; ++r)
                for (int c = 0; c < ref.cols; ++c) {
                    const size_t i = ref.at(r, c);
                    const int rr = r + dy, cc = c + dx;
                    if (ref.weight[i] <= 0 || rr < 0 || cc < 0 || rr >= mov.rows || cc >= mov.cols) continue;
                    const size_t j = mov.at(rr, cc);
                    if (mov.weight[j] <= 0) continue;
                    const double w = 1 / (1 / double(ref.weight[i]) + 1 / double(mov.weight[j]) +
                                          modelError * modelError);
                    const double res = double(mov.value[j]) - ref.value[i];
                    values.push_back({res, w});
                    sw += w;
                    swr += w * res;
                }
            if (values.empty()) continue;
            const double bias = swr / sw;
            double cost = 0;
            for (const auto& [res, w] : values) cost += std::min((res - bias) * (res - bias) * w, 25.0);
            candidates.push_back({dy, dx, cost / values.size(), bias, values.size()});
            mostPoints = std::max(mostPoints, values.size());
        }
    for (const auto& c : candidates) {
        if (c.points * 10 < mostPoints * 7 || c.cost >= bestCost) continue;
        bestCost = c.cost;
        best = {double(c.dy), double(c.dx), c.bias};
    }
    return best;
}

// A camera bump or tripod creep moves the image rigidly: a translation and a
// roll about the frame centre. Pan and tilt through a fisheye vary the shift
// across the frame only to second order, which per-tile values keep.
struct CameraModel {
    double ty = 0, tx = 0, roll = 0;  // Sites at the centre; roll in sites per 1000 sites
    double cy = 0, cx = 0;            // Frame centre, sites
    size_t inliers = 0;
    double rms = 0;                   // Inlier tiles' departure from the model, sites
    bool spread = false;              // Inliers cover enough of the frame to constrain it
    std::array<double, 2> at(double row, double col) const {
        return {ty + roll * (col - cx) / 1000, tx - roll * (row - cy) / 1000};
    }
};

struct PairResult {
    Params global;
    Fit globalFit;
    MotionField field;  // Displacement of the moving frame against the reference frame, sites
    CameraModel camera;
    bool confident = false;
};

// A robust rigid fit to the measured tiles separates camera motion from scene
// motion (clouds, foliage, a moving blocker, the sun): tiles that depart from
// the fit by more than their noise allows are dropped.
CameraModel rejectSceneMotion(const MotionField& field, double cy, double cx, const std::vector<double>& dy,
                              const std::vector<double>& dx, const std::vector<float>& stds,
                              std::vector<uint8_t>& known) {
    CameraModel model;
    model.cy = cy;
    model.cx = cx;
    std::vector<size_t> tiles;
    for (size_t t = 0; t < known.size(); ++t)
        if (known[t]) tiles.push_back(t);
    if (tiles.empty()) return model;
    const bool rotation = tiles.size() >= 4;
    std::vector<double> robust(tiles.size(), 1);
    auto centre = [&](size_t t, double& y, double& x) {
        y = (double(t / size_t(field.tileCols)) + 0.5) * field.tileSites;
        x = (double(t % size_t(field.tileCols)) + 0.5) * field.tileSites;
    };
    auto residual = [&](size_t k, double& ry, double& rx) {
        double y, x;
        centre(tiles[k], y, x);
        const auto d = model.at(y, x);
        ry = dy[tiles[k]] - d[0];
        rx = dx[tiles[k]] - d[1];
    };
    auto tolerance = [&](size_t k) {
        const double s = stds[tiles[k]];
        return std::sqrt(s * s + fieldTolerance * fieldTolerance);
    };
    // Start from the median translation and fit translation alone first, so the
    // majority of tiles decides which are inliers. Fitted together from the
    // start, roll can tilt the model through a cluster of moving tiles on one
    // side of the frame (beside the sun, say) and shift its centre.
    {
        std::vector<double> ys, xs;
        for (const size_t t : tiles) {
            ys.push_back(dy[t]);
            xs.push_back(dx[t]);
        }
        model.ty = median(ys);
        model.tx = median(xs);
        for (size_t k = 0; k < tiles.size(); ++k) {
            double ry, rx;
            residual(k, ry, rx);
            const double z = std::hypot(ry, rx) / tolerance(k) / cauchy;
            robust[k] = 1 / (1 + z * z);
        }
    }
    for (int it = 0; it < 20; ++it) {
        const bool roll = rotation && it >= 10;
        double h[3][3]{}, g[3]{};
        for (size_t k = 0; k < tiles.size(); ++k) {
            double y, x;
            centre(tiles[k], y, x);
            const double s = tolerance(k), w = robust[k] / (s * s);
            const double u = roll ? (x - cx) / 1000 : 0, v = roll ? -(y - cy) / 1000 : 0;
            // Rows equation: dy = ty + roll·u; columns equation: dx = tx + roll·v.
            const double jy[3] = {1, 0, u}, jx[3] = {0, 1, v};
            for (int a = 0; a < 3; ++a) {
                g[a] += w * (jy[a] * dy[tiles[k]] + jx[a] * dx[tiles[k]]);
                for (int b = 0; b < 3; ++b) h[a][b] += w * (jy[a] * jy[b] + jx[a] * jx[b]);
            }
        }
        if (!roll) h[2][2] = 1;
        double p[3];
        if (!solve3(h, g, p)) return model;
        model.ty = p[0];
        model.tx = p[1];
        model.roll = roll ? p[2] : 0;
        // The tolerance already states the spread camera motion allows, so the
        // scale is fixed: one estimated from all tiles grows with the share of
        // moving ones and lets them pull the fit.
        for (size_t k = 0; k < tiles.size(); ++k) {
            double ry, rx;
            residual(k, ry, rx);
            const double z = std::hypot(ry, rx) / tolerance(k) / cauchy;
            robust[k] = 1 / (1 + z * z);
        }
    }
    double sum = 0;
    int top = field.tileRows, bottom = -1, left = field.tileCols, right = -1;
    for (size_t k = 0; k < tiles.size(); ++k) {
        double ry, rx;
        residual(k, ry, rx);
        const double distance = std::hypot(ry, rx);
        if (distance > std::max(outlierFloor, 3 * tolerance(k))) {
            known[tiles[k]] = 0;
            continue;
        }
        ++model.inliers;
        sum += distance * distance;
        const int r = int(tiles[k] / size_t(field.tileCols)), c = int(tiles[k] % size_t(field.tileCols));
        top = std::min(top, r);
        bottom = std::max(bottom, r);
        left = std::min(left, c);
        right = std::max(right, c);
    }
    model.rms = model.inliers ? std::sqrt(sum / double(model.inliers)) : 0;
    // Inliers clustered in one place (the sky around the sun in an ND frame) are
    // usually the sun's or a cloud's own motion, not the camera's. Frames of a
    // few tiles only support a translation, from most of them.
    if (known.size() < 12) model.spread = model.inliers * 2 >= known.size();
    else model.spread = model.inliers >= std::max<size_t>(minimumInliers, known.size() / 25) &&
                        (bottom - top + 1) * 10 >= 4 * field.tileRows &&
                        (right - left + 1) * 10 >= 4 * field.tileCols;
    return model;
}

PairResult alignPair(const Pyramid& ref, const Pyramid& mov, const Options& o) {
    PairResult result;
    const size_t top = std::min(ref.size(), mov.size()) - 1;
    Params p = coarseSearch(ref[top], mov[top], 3);
    for (size_t l = top + 1; l-- > 0;) {
        const auto& a = ref[l];
        // The global fit only seeds the tiles: about 65,000 points per level do.
        const int stride = std::max(1, int(std::sqrt(double(a.rows) * a.cols / 65536)));
        const auto fit = gaussNewton(a, mov[l], 0, a.rows, 0, a.cols, stride, p, 30, o.threads, false);
        if (fit.points >= 16) {
            p = fit.p;
            if (l == 0) result.globalFit = fit;
        }
        if (l) {
            p.dy *= 2;
            p.dx *= 2;
        }
    }
    result.global = p;

    auto& field = result.field;
    const auto& base = ref[0];
    field.tileSites = 2 * tileQuads;
    field.tileRows = std::max(1, (base.rows + tileQuads - 1) / tileQuads);
    field.tileCols = std::max(1, (base.cols + tileQuads - 1) / tileQuads);
    const size_t tiles = size_t(field.tileRows) * field.tileCols;
    std::vector<double> dy(tiles, 2 * p.dy), dx(tiles, 2 * p.dx);
    std::vector<uint8_t> known(tiles, 0);
    std::vector<float> stds(tiles, INFINITY);
    parallelRows(field.tileRows, o.threads, [&](int tr) {
        for (int tc = 0; tc < field.tileCols; ++tc) {
            Params q{p.dy / 2, p.dx / 2, p.bias};
            Fit fit;
            for (int l = ref.size() > 1 ? 1 : 0; l >= 0; --l) {
                const auto& a = ref[size_t(l)];
                const int size = tileQuads >> l;
                const int r0 = tr * size, c0 = tc * size;
                const int r1 = std::min(a.rows, r0 + size), c1 = std::min(a.cols, c0 + size);
                fit = gaussNewton(a, mov[size_t(l)], r0, r1, c0, c1, 1, q, 20, 1, l == 0);
                if (fit.points < 16) break;
                q = fit.p;
                if (l) {
                    q.dy *= 2;
                    q.dx *= 2;
                }
            }
            const size_t t = size_t(tr) * field.tileCols + tc;
            const int area = std::min(tileQuads, base.rows - tr * tileQuads) *
                             std::min(tileQuads, base.cols - tc * tileQuads);
            const bool accepted = fit.converged && fit.points * 5 >= size_t(area) &&
                                  fit.std <= acceptStd &&
                                  std::hypot(fit.p.dy - p.dy, fit.p.dx - p.dx) <= acceptJump;
            if (!accepted) continue;
            dy[t] = 2 * fit.p.dy;
            dx[t] = 2 * fit.p.dx;
            stds[t] = float(2 * fit.std);
            known[t] = 1;
        }
    });
    // Tiles that depart from the camera model coherently with their neighbours
    // moved with the scene (clouds drift several sites a second at 8 mm); lone
    // ones are bad fits (a validity-mask edge, glare, the dim rim of the image
    // circle). Both are only labelled: the field holds camera motion, measured
    // where tiles agree with it and modelled elsewhere. Whether a merge should
    // follow scene motion is decided per pixel by its robustness test.
    auto inliers = known;
    result.camera = rejectSceneMotion(field, base.rows, base.cols, dy, dx, stds, inliers);
    result.confident = result.camera.spread;
    std::vector<uint8_t> scene(tiles, 0);
    for (int r = 0; r < field.tileRows; ++r)
        for (int c = 0; c < field.tileCols; ++c) {
            const size_t t = size_t(r) * field.tileCols + c;
            if (!known[t] || inliers[t]) continue;
            int agree = 0;
            for (int nr = std::max(0, r - 1); nr <= std::min(field.tileRows - 1, r + 1); ++nr)
                for (int nc = std::max(0, c - 1); nc <= std::min(field.tileCols - 1, c + 1); ++nc) {
                    const size_t n = size_t(nr) * field.tileCols + nc;
                    agree += n != t && known[n] && std::hypot(dy[n] - dy[t], dx[n] - dx[t]) <= outlierFloor;
                }
            scene[t] = agree >= 2;
        }
    for (size_t t = 0; t < tiles; ++t) {
        if (inliers[t]) continue;
        const auto d = result.camera.at((double(t / size_t(field.tileCols)) + 0.5) * field.tileSites,
                                        (double(t % size_t(field.tileCols)) + 0.5) * field.tileSites);
        dy[t] = result.confident ? d[0] : 2 * p.dy;
        dx[t] = result.confident ? d[1] : 2 * p.dx;
        stds[t] = INFINITY;
    }
    field.measured = std::move(inliers);
    field.scene = std::move(scene);
    field.dy.assign(dy.begin(), dy.end());
    field.dx.assign(dx.begin(), dx.end());
    field.std = std::move(stds);
    return result;
}

std::vector<uint8_t> coverage(const RawFrame& f, const Options& o) {
    const int rows = f.layout.height / coverageStride, cols = f.layout.width / coverageStride;
    std::vector<uint8_t> result(size_t(rows) * cols, 0);
    parallelRows(rows, o.threads, [&](int r) {
        for (int c = 0; c < cols; ++c) {
            double mean, variance;
            result[size_t(r) * cols + c] =
                matchable(f, r * coverageStride / 2, c * coverageStride / 2, o, mean, variance);
        }
    });
    return result;
}

} // namespace

std::array<double, 2> MotionField::at(double row, double col) const {
    if (dy.empty()) return {0, 0};
    const double fy = std::clamp((row + 0.5) / tileSites - 0.5, 0.0, double(tileRows - 1));
    const double fx = std::clamp((col + 0.5) / tileSites - 0.5, 0.0, double(tileCols - 1));
    const int r0 = int(fy), c0 = int(fx);
    const int r1 = std::min(r0 + 1, tileRows - 1), c1 = std::min(c0 + 1, tileCols - 1);
    const double wy = fy - r0, wx = fx - c0;
    auto lerp = [&](const std::vector<float>& v) {
        auto at = [&](int r, int c) { return double(v[size_t(r) * tileCols + c]); };
        return (1 - wy) * ((1 - wx) * at(r0, c0) + wx * at(r0, c1)) + wy * ((1 - wx) * at(r1, c0) + wx * at(r1, c1));
    };
    return {lerp(dy), lerp(dx)};
}

Alignment estimateAlignment(const FrameGroups& groups, size_t reference, const Options& o, Warnings& warnings) {
    Alignment result;
    struct Node { size_t group, frame; const RawFrame* raw; };
    std::vector<Node> nodes;
    for (size_t g = 0; g < groups.size(); ++g)
        for (size_t i = 0; i < groups[g]->size(); ++i) nodes.push_back({g, i, &(*groups[g])[i]});
    const size_t root = reference;  // Group 0 comes first
    std::vector<std::vector<uint8_t>> covered;
    for (const auto& node : nodes) covered.push_back(coverage(*node.raw, o));
    auto overlap = [&](size_t a, size_t b) {
        size_t n = 0;
        for (size_t i = 0; i < covered[a].size(); ++i) n += covered[a][i] & covered[b][i];
        return n;
    };
    const size_t minimum = std::max<size_t>(64, covered.front().size() / 100);

    result.frames.resize(nodes.size());
    for (size_t k = 0; k < nodes.size(); ++k) {
        auto& entry = result.frames[k];
        entry.group = nodes[k].group;
        entry.frame = nodes[k].frame;
        entry.coverage = std::count(covered[k].begin(), covered[k].end(), 1);
    }
    auto& rootEntry = result.frames[root];
    rootEntry.reference = rootEntry.aligned = true;
    std::vector<uint8_t> done(nodes.size(), 0);
    done[root] = 1;

    // The two most recent pyramids are kept: a chain of neighbours reuses them.
    std::vector<std::pair<size_t, std::shared_ptr<const Pyramid>>> cache;
    auto pyramid = [&](size_t k) {
        for (const auto& [id, built] : cache)
            if (id == k) return built;
        auto built = std::make_shared<const Pyramid>(buildPyramid(*nodes[k].raw, o));
        cache.emplace_back(k, built);
        if (cache.size() > 2) cache.erase(cache.begin());
        return built;
    };

    const auto& size = nodes.front().raw->layout;
    // Estimates `from` against the aligned `to` and composes it with `to`'s motion.
    auto link = [&](size_t from, size_t to, size_t shared) {
        const auto ref = pyramid(to);
        const auto mov = pyramid(from);
        auto pair = alignPair(*ref, *mov, o);
        auto& entry = result.frames[from];
        const auto& parent = result.frames[to];
        done[from] = 1;
        entry.alignedTo = to;
        entry.overlap = shared;
        entry.gainEv = parent.gainEv + pair.global.bias / std::log(2.0);
        entry.residualScale = pair.globalFit.scale;
        entry.inlierTiles = pair.camera.inliers;
        entry.modelRms = pair.camera.rms;
        const auto centre = pair.camera.at(size.height / 2.0, size.width / 2.0);
        entry.shift = {parent.shift[0] + centre[0], parent.shift[1] + centre[1]};
        entry.roll = parent.roll + pair.camera.roll;
        entry.uncertainty = std::hypot(parent.uncertainty, 2 * pair.globalFit.std);
        entry.field = std::move(pair.field);
        if (!parent.reference)
            for (size_t t = 0; t < entry.field.dy.size(); ++t) {
                entry.field.dy[t] += parent.field.dy[t];
                entry.field.dx[t] += parent.field.dx[t];
            }
        entry.aligned = pair.confident;
        entry.unconstrained = !pair.confident;
        return pair.confident;
    };
    // Grows a group's tree from its aligned frames, always pairing the frames
    // that share the most valid quads. Pairs within a bracket see the same glare
    // and scene, so they match far better than pairs across groups.
    auto grow = [&](size_t g) {
        for (;;) {
            size_t from = SIZE_MAX, to = SIZE_MAX, best = 0;
            for (size_t u = 0; u < nodes.size(); ++u) {
                if (done[u] || nodes[u].group != g) continue;
                for (size_t a = 0; a < nodes.size(); ++a) {
                    if (!result.frames[a].aligned || nodes[a].group != g) continue;
                    const size_t n = overlap(u, a);
                    if (n > best) {
                        best = n;
                        from = u;
                        to = a;
                    }
                }
            }
            if (from == SIZE_MAX || best < minimum) return;
            link(from, to, best);
        }
    };
    grow(0);
    // Each later group is tied to the earlier ones by its single best pair (up to
    // four candidates are tried in order of overlap), then grown inside itself.
    for (size_t g = 1; g < groups.size(); ++g) {
        std::vector<std::array<size_t, 3>> pairs;
        for (size_t u = 0; u < nodes.size(); ++u) {
            if (nodes[u].group != g) continue;
            for (size_t a = 0; a < nodes.size(); ++a) {
                if (!result.frames[a].aligned || nodes[a].group >= g) continue;
                const size_t n = overlap(u, a);
                if (n >= minimum) pairs.push_back({n, u, a});
            }
        }
        std::sort(pairs.rbegin(), pairs.rend());
        size_t attempts = 0;
        for (const auto& [n, u, a] : pairs) {
            if (done[u]) continue;
            if (link(u, a, n) || ++attempts == 4) break;
        }
        grow(g);
    }
    // A frame too dark or too bright to match, or whose estimate was not
    // constrained, takes the motion of the nearest exposure in its own bracket,
    // captured moments apart.
    size_t unconstrained = 0, isolated = 0;
    for (size_t k = 0; k < nodes.size(); ++k) {
        auto& entry = result.frames[k];
        if (entry.aligned) continue;
        size_t nearest = SIZE_MAX;
        auto gap = [&](size_t i) { return std::abs(int(result.frames[i].frame) - int(entry.frame)); };
        for (size_t a = 0; a < nodes.size(); ++a) {
            const auto& other = result.frames[a];
            if (!other.aligned || other.inherited || other.group != entry.group) continue;
            if (nearest == SIZE_MAX || gap(a) < gap(nearest)) nearest = a;
        }
        if (nearest == SIZE_MAX) {
            entry.unconstrained ? ++unconstrained : ++isolated;
            continue;
        }
        const auto& source = result.frames[nearest];
        entry.aligned = entry.inherited = true;
        entry.alignedTo = nearest;
        entry.shift = source.shift;
        entry.roll = source.roll;
        entry.gainEv = NAN;
        entry.uncertainty = source.uncertainty;
        entry.field = source.field;
    }
    for (size_t k = 0; k < nodes.size(); ++k) {
        auto& entry = result.frames[k];
        if (!entry.aligned || entry.reference) continue;
        std::vector<double> lengths{std::hypot(entry.shift[0], entry.shift[1])};
        for (size_t t = 0; t < entry.field.dy.size(); ++t)
            lengths.push_back(std::hypot(entry.field.dy[t], entry.field.dx[t]));
        entry.measuredTiles = size_t(std::count(entry.field.measured.begin(), entry.field.measured.end(), 1));
        entry.maxShift = *std::max_element(lengths.begin(), lengths.end());
        entry.medianShift = median(lengths);
    }
    if (unconstrained)
        warnings.push_back(std::to_string(unconstrained) + (unconstrained == 1 ? " frame has" : " frames have") +
                           " too few consistent tiles spread over the frame to tell camera motion from "
                           "scene motion (the sun, clouds); its motion is not estimated");
    if (isolated)
        warnings.push_back(std::to_string(isolated) + (isolated == 1 ? " frame shares" : " frames share") +
                           " too few valid quads with any aligned frame to estimate its motion");
    return result;
}
}
