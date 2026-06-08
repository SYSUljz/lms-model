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
file(GLOB_RECURSE LMS_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_library(lms_engine_lib ${LMS_SOURCES})
target_include_directories(lms_engine_lib PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lms_engine_lib PRIVATE GDAL::GDAL)
# I/O smoke test: read the model tiles and print a summary.
add_executable(lms_engine test/main.cpp)

target_include_directories(lms_engine PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lms_engine PRIVATE GDAL::GDAL lms_engine_lib)
