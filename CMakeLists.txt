cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)
get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -std=c++20 -fopenmp -O")

message(STATUS "-------------------------------")
message(STATUS "Processing Project ${PROJECT_NAME}:")

file(GLOB_RECURSE SOURCE_FILES 
     src/*.cpp 
     src/*.hpp 
     src/*.h
     libs/lodepng-master/*.cpp
     libs/lodepng-master/*.h
     libs/tinyexr/*.c
     libs/tinyexr/*.h
     )

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
add_subdirectory(libs/Walnut/)

add_executable (mrt ${SOURCE_FILES})

target_link_libraries(mrt PUBLIC walnut)

find_package (OpenImageIO REQUIRED)
target_link_libraries(mrt PRIVATE OpenImageIO::OpenImageIO)

target_include_directories (mrt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/lodepng-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/tinyexr/)
