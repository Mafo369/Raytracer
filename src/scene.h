#pragma once

#include "defines.h"
#include "tiny_obj_loader.h"

//#include "Object.h"
#include "textures.hpp"

class Camera;
class Light;
class Object;
class Material;

typedef std::vector<Object*> Objects;
typedef std::vector<Light*> Lights;

typedef struct scene_s {
  Lights lights; //! the scene have several lights
  Objects objects; //! the scene have several objects
  Camera* cam; //! the scene have one camera
  color3 skyColor; //! the sky color, could be extended to a sky function ;)
  texture* m_skyTexture;
} Scene;

//! create a new sphere structure
Object* initSphere(std::shared_ptr<Material> mat, Transform transform);
Object *initCube(std::shared_ptr<Material> mat, Transform transform);
Object* initPlane(vec3 normal, float d, std::shared_ptr<Material> mat);
Object* initTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], std::shared_ptr<Material> mat);
Object *initSmoothTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], vec3 n1, vec3 n2, vec3 n3, Transform transform, std::shared_ptr<Material> mat);

//! release memory for the object obj
void freeObject(Object *obj);

//! init a new light at position with a give color (no special unit here for the moment)
Light* initPointLight(point3 position, color3 color);
Light *initDirectLight(vec3 direction, color3 color);

Light* initAmbientLight(color3 color);

//! release memory for the light
void freeLight(Light *);

// allocate the momery for the scene
Scene *initScene();
void freeScene(Scene *scene);

void setCameraFOV(Scene *scene, point3 position, vec3 at, vec3 up, float fov, float width, float height, float aperture = 0.01, float dist_to_focus = 1);
void setSimpleCamera(Scene *scene, point3 position, vec3 at, vec3 up, float fov, float w, float h, float radius, float distFocus);

//! take ownership of obj freeScene will free obj) ... typically use addObject(scene, initPlane()
void addObject(Scene *scene, Object *obj);

//! take ownership of light : freeScene will free light) ... typically use addObject(scene, initLight()
void addLight(Scene *scene, Light *light);

void setSkyColor(Scene *scene, color3 c);

void readObjToTriangleMesh(const char *file, tinyobj::attrib_t &attrib, std::vector<tinyobj::shape_t> &shapes, std::vector<tinyobj::material_t> &materials);
