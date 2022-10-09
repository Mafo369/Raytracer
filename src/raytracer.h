#pragma once

#include "defines.h"
#include "image.h"
#include "scene.h"
#include "ray.h"


//! An intersection contains all the information to shade an intersection point.
typedef struct intersection_s { 
  vec3 normal; //! the normal of the intersection point
  point3 position; //! the intersection point
  std::shared_ptr<Material> mat; //! the material of th intersected object

  vec2 t;
  float u;
  float v;

  glm::mat4 transform;

  bool isOutside;
  int face = -1;
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
color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection);
