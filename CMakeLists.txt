cmake_minimum_required(VERSION 3.18)
project(rawfuse VERSION 0.3.0 LANGUAGES CXX)

include(CTest)
include(GNUInstallDirs)
find_package(PkgConfig REQUIRED)
pkg_check_modules(LibRaw REQUIRED IMPORTED_TARGET libraw>=0.21)
find_package(OpenImageIO CONFIG REQUIRED)
find_package(OpenEXR 3 CONFIG REQUIRED)
find_package(RAWTOACES 2.2 CONFIG REQUIRED COMPONENTS util)
if(NOT TARGET rawtoaces_util)
    message(FATAL_ERROR "RAWTOACES must export its C++ utility library")
endif()
find_package(nlohmann_json 3.9 CONFIG REQUIRED)
find_package(Threads REQUIRED)

add_library(rawfuse_pipeline STATIC
    src/options.cpp src/metadata.cpp src/raw_frame.cpp src/dark_model.cpp src/merge.cpp
    src/bracket.cpp src/nd.cpp src/groups.cpp src/finish.cpp src/colour.cpp src/writer.cpp src/pipeline.cpp)
target_compile_features(rawfuse_pipeline PUBLIC cxx_std_17)
target_compile_definitions(rawfuse_pipeline PUBLIC RAWFUSE_VERSION="${PROJECT_VERSION}")
target_link_libraries(rawfuse_pipeline PRIVATE
    PkgConfig::LibRaw OpenImageIO::OpenImageIO OpenEXR::OpenEXR
    rawtoaces_util nlohmann_json::nlohmann_json Threads::Threads)
target_include_directories(rawfuse_pipeline PUBLIC src)
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rawfuse_pipeline PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(rawfuse src/main.cpp)
target_link_libraries(rawfuse PRIVATE rawfuse_pipeline)
if(APPLE)
    set(rawfuse_relative_rpath "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
elseif(UNIX)
    set(rawfuse_relative_rpath "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
set_target_properties(rawfuse PROPERTIES
    INSTALL_RPATH "${rawfuse_relative_rpath};$<TARGET_FILE_DIR:rawtoaces_util>;$<TARGET_FILE_DIR:rawtoaces_core>"
    INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS rawfuse RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(BUILD_TESTING)
    find_package(Python3 3.8 COMPONENTS Interpreter REQUIRED)
    add_test(NAME synthetic_dng COMMAND Python3::Interpreter
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/synthetic.py $<TARGET_FILE:rawfuse>)
    add_test(NAME pipeline_regressions COMMAND Python3::Interpreter
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/review.py $<TARGET_FILE:rawfuse>)
endif()
