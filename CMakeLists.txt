cmake_minimum_required(VERSION 3.20)
project(lms_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# GDAL ships a CMake package config (GDALConfig.cmake) exposing the GDAL::GDAL
# target, which carries include dirs + link libs automatically.
find_package(GDAL CONFIG REQUIRED)
file(GLOB_RECURSE LMS_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

add_subdirectory("/home/jack-li/Desktop/starlings" "${CMAKE_CURRENT_BINARY_DIR}/ex_opt")

add_library(lms_engine_lib ${LMS_SOURCES})
target_include_directories(lms_engine_lib PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lms_engine_lib PRIVATE GDAL::GDAL)
# I/O smoke test: read the model tiles and print a summary.
add_executable(simulate_test test/simulate.cpp)

target_include_directories(simulate_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(simulate_test PRIVATE GDAL::GDAL lms_engine_lib)
add_executable(optimize_test test/optimize.cpp)
target_include_directories(optimize_test PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(optimize_test PRIVATE GDAL::GDAL lms_engine_lib ex_opt)