cmake_minimum_required(VERSION 3.20)
project(lms_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# GDAL ships a CMake package config (GDALConfig.cmake) exposing the GDAL::GDAL
# target, which carries include dirs + link libs automatically.
find_package(GDAL CONFIG REQUIRED)

# I/O smoke test: read the model tiles and print a summary.
add_executable(lms_io_test src/main.cpp)
target_include_directories(lms_io_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(lms_io_test PRIVATE GDAL::GDAL)
