#ifndef __SCENE_H__
#define __SCENE_H__

#include "defines.h"

// SCENE
typedef struct scene_s Scene;
typedef struct object_s Object;
typedef struct light_s Light;
typedef struct camera_s Camera;

typedef struct material_s {
  float IOR;	//! Index of refraction (for dielectric)
  float roughness; //! 0.001 - 0.01 : very smooth finish with slight imperfections. 0.1 : relatively rough. 0.3-0.7 extremely rough 
  color3 specularColor;	//! Specular "albedo"
  color3 diffuseColor;	//! Base color
} Material;

enum Etype {SPHERE=1, PLANE};


//! create a new sphere structure
Object* initSphere(point3 center, float radius, Material mat);
Object* initPlane(vec3 normal, float d, Material mat);

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


#endif
