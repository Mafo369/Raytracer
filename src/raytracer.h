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
  Transform transform;

  bool isOutside;
  int face = -1;

  float u;
  float v;

  vec3 dpdv;
  vec3 dn[2];
  vec3 dpdu;
  float dudx = 0;
  float dvdx = 0;
  float dudy = 0;
  float dvdy = 0;
  vec3 dpdx;
  vec3 dpdy;

  bool parametric = true;
  bool hit = false;

  // From PBR..

  bool SolveLinearSystem2x2(const float A[2][2],
          const float B[2], float *x0, float *x1) {
      float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
      if (std::abs(det) < 1e-10f)
          return false;
      *x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
      *x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
      if (std::isnan(*x0) || std::isnan(*x1))
          return false;
      return true;
  }


  void computeDifferentials(Ray* ray){
    if(parametric) {
      float d = -dot(normal, position);
      float tx = (-dot(normal, vec3(ray->dox)) - d) / dot(normal, ray->ddx);
      float ty = (-dot(normal, vec3(ray->doy)) - d) / dot(normal, ray->ddy);
      point3 px = ray->dox + tx * ray->ddx;
      point3 py = ray->doy + ty * ray->ddy;

      dpdx = px - position;
      dpdy = py - position;

      int dim[2];
      if (std::abs(normal.x) > std::abs(normal.y) && std::abs(normal.x) > std::abs(normal.z)) {
           dim[0] = 1; dim[1] = 2;    
       } else if (std::abs(normal.y) > std::abs(normal.z)) {
           dim[0] = 0; dim[1] = 2;    
       } else {
           dim[0] = 0; dim[1] = 1;
       }

      float A[2][2] = { { dpdu[dim[0]], dpdv[dim[0]] },
                        { dpdu[dim[1]], dpdv[dim[1]] } };
      float Bx[2] = { px[dim[0]] - position[dim[0]], px[dim[1]] - position[dim[1]] };
      float By[2] = { py[dim[0]] - position[dim[0]], py[dim[1]] - position[dim[1]] };

      if(!SolveLinearSystem2x2(A, Bx, &dudx, &dvdx))
        dudx = dvdx = 0;
      if(!SolveLinearSystem2x2(A, By, &dudy, &dvdy))
        dudy = dvdy = 0;
    }

    dudx *= 2.f;
    dvdx *= 2.f;
    dudy *= 2.f;
    dvdy *= 2.f;
  }

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

class Sampler;

void renderImage(RenderImage *img, Scene *scene);
color3 trace_ray(Scene *scene, Ray *ray, KdTree *tree, Intersection* intersection, bool show_lights=true, Sampler* sampler=nullptr);
