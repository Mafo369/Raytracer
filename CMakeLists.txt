project(RayTracer)

set(CMAKE_CXX_FLAGS "-Wall -std=c++20 -fopenmp -O2")

set(mrt_SRCS
  src/image.cpp
  src/kdtree.cpp
  libs/lodepng-master/lodepng.cpp
  src/main.cpp
  src/raytracer.cpp
  src/scene.cpp
  src/bsdf.cpp
  src/example_scenes.cpp
  src/Light.cpp
  src/Object.cpp
  src/materials/CookTorrance.cpp
  src/materials/Blinn.cpp
  src/shapes/sphere.cpp
  src/shapes/triangle.cpp
  src/shapes/cube.cpp
  src/shapes/plane.cpp
  src/sampling/sampling.cpp
  src/sampling/sampler.cpp
  src/sampling/stratified.cpp
  src/Material.cpp
  src/mediums/Fog.cpp
  libs/tinyexr/miniz.c
  src/integrator.cpp
  src/LightDistribution.cpp
  src/application.cpp
  )

#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
#add_subdirectory(libs/Walnut/)
find_package(Vulkan REQUIRED FATAL_ERROR) # error
find_package(glfw3 REQUIRED FATAL_ERROR)

add_executable (mrt ${mrt_SRCS})
target_link_libraries(mrt PUBLIC vulkan glfw)

#find_package (OpenImageIO REQUIRED)
#target_link_libraries(mrt PRIVATE OpenImageIO::OpenImageIO)



target_include_directories (mrt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/lodepng-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/tinyexr/ ${VULKAN_INCLUDE_DIRS})


set(unit_test_SRCS
  src/image.cpp
  src/kdtree.cpp
  libs/lodepng-master/lodepng.cpp
  src/raytracer.cpp
  src/scene.cpp
  src/intersection.cpp
  src/bsdf.cpp
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
