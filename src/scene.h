#pragma once

#include "defines.h"
#include "tiny_obj_loader.h"

//#include "Object.h"
#include "textures.hpp"

class Light;

// SCENE
//typedef struct scene_s Scene;
////typedef struct object_s Object;
////typedef struct light_s Light;
//typedef struct camera_s Camera;

typedef struct camera_s {
  point3 position; //! eye position
    vec3 zdir; //! view direction
    vec3 xdir; //! right direction
    vec3 ydir; //! up direction
    point3 center; //! center of the image plane
  float fov;  //! field of view
  float aspect; //! aspect ratio (typically use WIDTH/HEIGHT of the computed image
} Camera;


//typedef struct object_s {
//    glm::mat4 transform;
//    glm::mat4 invTransform;
//    Geometry geom;
//    Material mat;
//} Object;

class Object;

typedef std::vector<Object*> Objects;
typedef std::vector<Light*> Lights;

typedef struct scene_s {
  Lights lights; //! the scene have several lights
  Objects objects; //! the scene have several objects
  Camera cam; //! the scene have one camera
  color3 skyColor; //! the sky color, could be extended to a sky function ;)
} Scene;

typedef struct material_s Material;


//! create a new sphere structure
Object *initSphere(point3 center, float radius, Material mat);
Object* initSphere(Material mat, glm::mat4 transform = glm::mat4(1.f));
Object *initCube(Material mat, glm::mat4 transform = glm::mat4(1.f));
Object* initPlane(vec3 normal, float d, Material mat);
Object* initTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], Material mat);
Object *initSmoothTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], vec3 n1, vec3 n2, vec3 n3, Material mat);

//! release memory for the object obj
void freeObject(Object *obj);

//! init a new light at position with a give color (no special unit here for the moment)
Light* initPointLight(point3 position, color3 color);

//! release memory for the light
void freeLight(Light *);

// allocate the momery for the scene
Scene *initScene();
void freeScene(Scene *scene);

void setCamera(Scene *scene, point3 position, vec3 at, vec3 up, float fov, float aspect);

//! take ownership of obj freeScene will free obj) ... typically use addObject(scene, initPlane()
void addObject(Scene *scene, Object *obj);

//! take ownership of light : freeScene will free light) ... typically use addObject(scene, initLight()
void addLight(Scene *scene, Light *light);

void setSkyColor(Scene *scene, color3 c);

void readObjToTriangleMesh(const char *file, tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes, std::vector<tinyobj::material_t> &materials);
