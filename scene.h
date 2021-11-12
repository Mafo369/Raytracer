#ifndef __SCENE_H__
#define __SCENE_H__

#include "defines.h"
#include "tiny_obj_loader.h"

// SCENE
typedef struct scene_s Scene;
typedef struct object_s Object;
typedef struct light_s Light;
typedef struct camera_s Camera;

enum matType {LAMBERTIAN=1, METAL};

typedef struct material_s {
  float IOR;	//! Index of refraction (for dielectric)
  float roughness; //! 0.001 - 0.01 : very smooth finish with slight imperfections. 0.1 : relatively rough. 0.3-0.7 extremely rough 
  color3 specularColor;	//! Specular "albedo"
  color3 diffuseColor;	//! Base color
  matType type;
  float fuzz;
} Material;

enum Etype {SPHERE=1, PLANE, TRIANGLE};

//! create a new sphere structure
Object* initSphere(point3 center, float radius, Material mat);
Object* initPlane(vec3 normal, float d, Material mat);
Object* initTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], Material mat);

//! release memory for the object obj
void freeObject(Object *obj);

//! init a new light at position with a give color (no special unit here for the moment)
Light* initLight(point3 position, color3 color);

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


#endif
