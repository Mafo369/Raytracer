#pragma once

#include "defines.h"
#include "image.h"
#include "scene.h"
#include "ray.h"


//! An intersection contains all the information to shade an intersection point.
typedef struct intersection_s { 
  vec3 normal; //! the normal of the intersection point
  point3 position; //! the intersection point
  Material *mat; //! the material of th intersected object

  vec2 t;
  float u;
  float v;

  bool isOutside;
} Intersection;

#include "kdtree.h"

/// test the ray intersection against each object of the scene, the nearest intersection
// is stored in the parameter intersection
// Possible intersection are considered only between ray->tmin and ray->tmax
// ray->tmax is updated during this process
bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection );
bool intersectCylinder (Ray *ray, Intersection *intersection, Object *cylinder);
bool intersectPlane(Ray *ray, Intersection *intersection, Object *plane);
bool intersectSphere(Ray *ray, Intersection *intersection, Object *sphere);
bool intersectTriangle(Ray *ray, Intersection *intersection, Object *sphere);

void renderImage(RenderImage *img, Scene *scene);
color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree);
color3 trace_ray_4multisampling(Scene *scene, KdTree *tree, int indexI, int indexJ, vec3 dx,
                               vec3 dy, vec3 ray_delta_x, vec3 ray_delta_y);
color3 trace_ray_multisampling(Scene *scene, KdTree *tree, int indexI, int indexJ, vec3 dx,
                               vec3 dy, vec3 ray_delta_x, vec3 ray_delta_y);

float RDM_Beckmann(float NdotH, float alpha);
float RDM_Fresnel(float LdotH, float extIOR, float intIOR);
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m);
color3 RDM_bsdf_d(Material *m);
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m);

