#ifndef __RAYTRACER_H__
#define __RAYTRACER_H__

#include "defines.h"
#include "image.h"
#include "scene.h"
#include "ray.h"



//! An intersection contains all the information to shade an intersection point.
typedef struct intersection_s { 
  vec3 normal; //! the normal of the intersection point
  point3 position; //! the intersection point
  Material *mat; //! the material of th intersected object

  bool isTex;
  vec2 t;

  bool front_face;
} Intersection;

float length(vec3 v);
float random_float();
float random_float(float min, float max);
float clamp(float x, float min, float max);
vec3 random_in_unit_disk();
float degrees_to_radians(float degrees);
vec3 unit_vector(vec3 v);
bool hit(Ray *r, float t_min, float t_max,Intersection *rec, Object *obj);
/// test the ray intersection against each object of the scene, the nearest intersection
// is stored in the parameter intersection
// Possible intersection are considered only between ray->tmin and ray->tmax
// ray->tmax is updated during this process
bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection );
bool intersectCylinder (Ray *ray, Intersection *intersection, Object *cylinder);
bool intersectPlane(Ray *ray, Intersection *intersection, Object *plane);
bool intersectSphere(Ray *ray, Intersection *intersection, Object *sphere);
bool intersectTriangle(Ray *ray, Intersection *intersection, Object *sphere);

void renderImage(Image *img, Scene *scene);

float RDM_Beckmann(float NdotH, float alpha);
float RDM_Fresnel(float LdotH, float extIOR, float intIOR);
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m);
color3 RDM_bsdf_d(Material *m);
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m);

#endif
