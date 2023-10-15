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

add_subdirectory(libs/imgui)

add_executable (mrt ${mrt_SRCS})
target_link_libraries(mrt PUBLIC vulkan glfw imgui)

#find_package (OpenImageIO REQUIRED)
#target_link_libraries(mrt PRIVATE OpenImageIO::OpenImageIO)

target_include_directories (mrt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/glm-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/lodepng-master/ ${CMAKE_CURRENT_SOURCE_DIR}/libs/tinyexr/ ${VULKAN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/libs/imgui/)


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

#--------------------------------------------------------------------------------------------------
# Shaders in project
#
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/utilities.cmake)

set(SHD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/shaders)

# Local shaders
file(GLOB SHD_SRC # Shader sources
    ${SHD_DIR}/*.vert 
    ${SHD_DIR}/*.frag 
    ${SHD_DIR}/*.rgen 
    ${SHD_DIR}/*.rchit 
    ${SHD_DIR}/*.rahit 
    ${SHD_DIR}/*.rmiss 
    ${SHD_DIR}/*.comp
    )
    
file(GLOB SHD_HDR  # Shader headers
    ${SHD_DIR}/*.glsl 
    ${SHD_DIR}/*.h
    )



# Compiling shaders to Spir-V header
compile_glsl(
    SOURCE_FILES ${SHD_SRC}
    HEADER_FILES ${SHD_HDR}
    DST "${CMAKE_CURRENT_SOURCE_DIR}/_autogen" 
    VULKAN_TARGET "vulkan1.3"
    HEADER ON
    DEPENDENCY ${VULKAN_BUILD_DEPENDENCIES}
    FLAGS -I${SHD_DIR} -g
    )

target_sources(mrt PRIVATE ${GLSL_SOURCES} ${GLSL_HEADERS})
source_group("Shaders/src"   FILES ${GLSL_SOURCES} )
source_group("Shaders/hdr"   FILES ${GLSL_HEADERS} )
