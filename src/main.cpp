#include "pipeline.h"

#include <iostream>
#include <string_view>

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--version") {
        std::cout << "rawfuse " RAWFUSE_VERSION "\n";
        return 0;
    }
    if (argc == 2 && (std::string_view(argv[1]) == "--help" || std::string_view(argv[1]) == "-h")) {
        rawfuse::printHelp();
        return 0;
    }
    try {
        rawfuse::run(rawfuse::parseOptions(argc, argv));
        return 0;
    } catch (const std::exception& error) {
        std::string message = error.what();
        for (auto& c : message) if (c == '\n' || c == '\r') c = ' ';
        std::cerr << "rawfuse: " << message << '\n';
        return 2;
    }
}
