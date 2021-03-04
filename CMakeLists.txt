cmake_minimum_required (VERSION 2.8.11)
project(RayTracer)

set(CMAKE_CXX_FLAGS "-Wall -std=c++14 -g -fopenmp -O")

set(mrt_SRCS
  image.cpp
  kdtree.cpp
  ./lodepng-master/lodepng.cpp
  main.cpp
  raytracer.cpp
  scene.cpp
  )

add_executable (mrt ${mrt_SRCS})
target_include_directories (mrt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/lodepng-master/ )


set(unit_test_SRCS
  image.cpp
  kdtree.cpp
  ./lodepng-master/lodepng.cpp
  unit-test.cpp
  raytracer.cpp
  scene.cpp
  )

add_executable (unit-test ${unit_test_SRCS})
target_include_directories (unit-test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/lodepng-master/ )


set(make_test_SRCS
  image.cpp
  kdtree.cpp
  ./lodepng-master/lodepng.cpp
  make-test.cpp
  raytracer.cpp
  scene.cpp
  )

#add_executable (make-test ${make_test_SRCS})
#target_include_directories (make-test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/glm-master/glm ${CMAKE_CURRENT_SOURCE_DIR}/lodepng-master/ )
