project(RayTracer)

set(CMAKE_CXX_FLAGS "-Wall -std=c++17 -g -fopenmp -O")

set(mrt_SRCS
  src/image.cpp
  src/kdtree.cpp
  libs/lodepng-master/lodepng.cpp
  src/main.cpp
  src/raytracer.cpp
  src/scene.cpp
  src/bsdf.cpp
  src/realtime/Camera.cpp
  src/realtime/Renderer.cpp
  src/example_scenes.cpp
  src/Light.cpp
  src/Object.cpp
  src/materials/CookTorrance.cpp
  src/materials/Blinn.cpp
  src/shapes/sphere.cpp
  src/shapes/triangle.cpp
  src/shapes/cube.cpp
  src/shapes/plane.cpp
  )

add_subdirectory(libs/Walnut/)

add_executable (mrt ${mrt_SRCS})

target_link_libraries(mrt PUBLIC walnut)

target_include_directories (mrt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/lodepng-master/ )


set(unit_test_SRCS
  src/image.cpp
  src/kdtree.cpp
  libs/lodepng-master/lodepng.cpp
  src/raytracer.cpp
  src/scene.cpp
  src/intersection.cpp
  src/bsdf.cpp
  src/realtime/Camera.cpp
  src/realtime/Renderer.cpp
  src/example_scenes.cpp
  src/Light.cpp
  src/unit-test.cpp
  src/Object.cpp
  src/shapes/sphere.cpp
  src/shapes/triangle.cpp
  src/shapes/cube.cpp
  src/shapes/plane.cpp
  )

#add_executable (unit-test ${unit_test_SRCS})
#target_link_libraries(unit-test PUBLIC walnut)
#target_include_directories (unit-test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/lodepng-master/ )


set(make_test_SRCS
  src/image.cpp
  src/kdtree.cpp
  libs/lodepng-master/lodepng.cpp
  src/make-test.cpp
  src/raytracer.cpp
  src/scene.cpp
  src/example_scenes.cpp
  )

#add_executable (make-test ${make_test_SRCS})
#target_include_directories (make-test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/glm-master/glm ${CMAKE_CURRENT_SOURCE_DIR}/lodepng-master/ )
