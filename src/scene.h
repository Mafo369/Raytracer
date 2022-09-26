#pragma once

#include "defines.h"
#include "tiny_obj_loader.h"

#include "textures.hpp"

class Light {
  public:
    virtual ~Light() = default;
    Light() = default;
    vec3 getColor() { return m_color; } 
    vec3 getPosition() { return m_position; }
  protected:
    color3 m_color;
    vec3 m_position;
};

class PointLight : public Light {
  public:
    PointLight(vec3 position , color3 color);
    ~PointLight();
  private:
};

// SCENE
typedef struct scene_s Scene;
typedef struct object_s Object;
//typedef struct light_s Light;
typedef struct camera_s Camera;

enum Mtype {DIFFUSE=1, DIELECTRIC};

typedef struct material_s {
  float IOR;	//! Index of refraction (for dielectric)
  float roughness; //! 0.001 - 0.01 : very smooth finish with slight imperfections. 0.1 : relatively rough. 0.3-0.7 extremely rough 
  color3 specularColor;	//! Specular "albedo"
  color3 diffuseColor;	//! Base color
  Mtype mtype = DIFFUSE;
  
  texture* m_texture = nullptr;
  
} Material;

enum Etype {SPHERE=1, PLANE, TRIANGLE};


//! create a new sphere structure
Object* initSphere(point3 center, float radius, Material mat);
Object* initPlane(vec3 normal, float d, Material mat);
Object* initTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], Material mat);
Object *initSmoothTriangle(vec3 p1, vec3 p2, vec3 p3, vec3 n, vec2 t[3], vec3 n1, vec3 n2, vec3 n3, Material mat);

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

//class AreaLight : public Light {
//  public:
//    // AreaLight Interface
//    AreaLight(int nbSamples);
//    virtual color3 L(const Intersection &intr, const vec3 &w) const = 0;
//};
//
//class DiffuseAreaLight : public AreaLight {
//  public:
//    // DiffuseAreaLight Public Methods
//    DiffuseAreaLight(const color3 &Lemit, int nSamples,
//                     const std::shared_ptr<Object> &shape,
//                     bool twoSided);
//    color3 L(const Intersection &intr, const vec3 &w) const {
//        return (twoSided || dot(intr.normal, w) > 0) ? Lemit : color3(0.f);
//    }
//    color3 Power() const;
//    color3 Sample_Li(const Intersection &ref, const point2 &u, vec3 *wo,
//                       float *pdf) const;
//    float Pdf_Li(const Intersection &, const vec3 &, Scene*, KdTree*) const;
//    color3 Sample_Le(const point2 &u1, const point2 &u2, float time,
//                       Ray *ray, vec3 *nLight, float *pdfPos,
//                       float *pdfDir) const;
//    void Pdf_Le(const Ray &, const vec3 &, float *pdfPos,
//                float *pdfDir) const;
//
//  protected:
//    // DiffuseAreaLight Protected Data
//    const color3 Lemit;
//    std::shared_ptr<Object> shape;
//    // Added after book publication: by default, DiffuseAreaLights still
//    // only emit in the hemimsphere around the surface normal.  However,
//    // this behavior can now be overridden to give emission on both sides.
//    const bool twoSided;
//    const float area;
//};
