#define TINYOBJLOADER_IMPLEMENTATION

#include <iostream>
#include <string>

#include "scene.h"
#include "scene_types.h"
#include <string.h>
#include <algorithm>

//#include <glm/gtx/norm.hpp>

#include "Light.h"

Object *initSphere(point3 center, float radius, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = SPHERE;
    ret->geom.sphere.center = point3(0,0,0);
    ret->geom.sphere.radius = 0.25;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>
glm::vec3 extractScale(const glm::mat4 &m)
{
    // length2 returns length squared i.e. v·v
    // no square root involved
    return glm::vec3(glm::length2( glm::vec3(m[0]) ),
                     glm::length2( glm::vec3(m[1]) ),
                     glm::length2( glm::vec3(m[2]) ));
}

Object *initSphere(glm::mat4 transform, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = SPHERE;
    ret->geom.sphere.center = transform * vec4(0,0,0,1);
    glm::vec3 scalesSq = extractScale(transform);
    float const maxScaleSq = *std::max_element(&scalesSq[0], &scalesSq[0] + scalesSq.length());  // length gives the dimension here i.e. 3
    // one sqrt when you know the largest of the three
    float const largestScale = std::sqrt(maxScaleSq);
    ret->geom.sphere.radius = largestScale;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    ret->transform = transform;
    ret->invTransform = glm::inverse(ret->transform);
    return ret;
}

Object *initPlane(vec3 normal, float d, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = PLANE;
    ret->geom.plane.normal = normalize(normal);
    ret->geom.plane.dist = d;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object *initCube(vec3 min, vec3 max, Material mat) {
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = CUBE;
    ret->geom.cube.min = min;
    ret->geom.cube.max = max;
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object *initTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], Material mat){
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = TRIANGLE;
    ret->geom.triangle.p1 = p1;
    ret->geom.triangle.p2 = p2;
    ret->geom.triangle.p3 = p3;
    ret->geom.triangle.normal = n;
    ret->geom.triangle.tex[0] = t[0];
    ret->geom.triangle.tex[1] = t[1];
    ret->geom.triangle.tex[2] = t[2];
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

Object *initSmoothTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], vec3 n1, vec3 n2, vec3 n3, Material mat){
    Object *ret;
    ret = (Object *)malloc(sizeof(Object));
    ret->geom.type = TRIANGLE;
    ret->geom.triangle.p1 = p1;
    ret->geom.triangle.p2 = p2;
    ret->geom.triangle.p3 = p3;
    ret->geom.triangle.n1 = n1;
    ret->geom.triangle.n2 = n2;
    ret->geom.triangle.n3 = n3;
    ret->geom.triangle.normal = n;
    ret->geom.triangle.tex[0] = t[0];
    ret->geom.triangle.tex[1] = t[1];
    ret->geom.triangle.tex[2] = t[2];
    memcpy(&(ret->mat), &mat, sizeof(Material));
    return ret;
}

void freeObject(Object *obj) {
    free(obj);
}

Light *initPointLight(point3 position, color3 color) {
    Light *light = new PointLight(position, color);
    return light;
}

void freeLight(Light *light) {
    delete light;
}

Scene * initScene() {
    return new Scene;
}

void freeScene(Scene *scene) {
    std::for_each(scene->objects.begin(), scene->objects.end(), freeObject);
    std::for_each(scene->lights.begin(), scene->lights.end(), freeLight);
    delete scene;
}

void setCamera(Scene *scene, point3 position, point3 at, vec3 up, float fov, float aspect) {
    scene->cam.fov = fov;
    scene->cam.aspect = aspect;
    scene->cam.position = position;
    scene->cam.zdir = normalize(at-position);
    scene->cam.xdir = normalize(cross(up, scene->cam.zdir));
    scene->cam.ydir = normalize(cross(scene->cam.zdir, scene->cam.xdir));
    scene->cam.center = 1.f / tanf ((scene->cam.fov * glm::pi<float>() / 180.f) * 0.5f) * scene->cam.zdir;
}

void addObject(Scene *scene, Object *obj) {
    scene->objects.push_back(obj);
}

void addLight(Scene *scene, Light *light) {
    scene->lights.push_back(light);
}

void setSkyColor(Scene *scene, color3 c) {
    scene->skyColor = c;
}

//https://github.com/tinyobjloader/tinyobjloader
void readObjToTriangleMesh(const char *file, tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes, std::vector<tinyobj::material_t> &materials){
  std::string inputfile = file;
  std::string err;

  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, inputfile.c_str(), "../assets/");

  if(!err.empty()){
      std::cout << err << std::endl;
  }

  if(!ret)
    exit(1);

  std::cout << "Model vertices: " << attrib.vertices.size() << std::endl;
}
