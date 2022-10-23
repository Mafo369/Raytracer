#include "sphere.h"

bool Sphere::intersect(Ray *ray, Intersection *intersection) const {
  Ray transformedRay = transformRay(ray);

  vec3 oc = transformedRay.orig - point3(0.f, 0.f, 0.f);
  float a = dot(transformedRay.dir, transformedRay.dir);
  float b = 2.f * dot(transformedRay.dir, oc);
  float c = dot(oc, oc) - 1.f;

  float delta = b * b - 4.f * a * c;

  if (delta == 0) {
    // Une solution
    float t = -b / (2.f * a);
    if (t >= ray->tmin && t <= ray->tmax) {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = mat;
      vec3 objectPoint = transform.transformTo(intersection->position);
      vec3 objectNormal = objectPoint - vec3(0.f, 0.f, 0.f);
      vec3 normal = transform.vectorTransformFrom(objectNormal);
      intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
      intersection->normal = normalize(normal);

      float pi = M_PI;
      auto theta = glm::atan(objectPoint.x, objectPoint.z);
      auto vec = glm::vec3(objectPoint.x, objectPoint.y, objectPoint.z);
      auto radius = glm::length(vec);

      auto phi = acos(objectPoint.y / radius);
      auto raw_u = theta / (2 * pi);
      vec3 uv = vec3(1 - (raw_u + 0.5), 1 - phi / pi, 0);

      intersection->u = uv.x;
      intersection->v = uv.y;

      vec3 d = normalize(transformedRay.dir);
      float _t = length(t * transformedRay.dir);
      vec3 dDx = ray->ddx;
      vec3 dDy = ray->ddy;

      vec3 dtx =
          -(ray->dox + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
      vec3 dty =
          -(ray->doy + _t * dot(dDy, objectNormal) / dot(d, objectNormal));

      // delta hit point on plane
      vec3 dXx = ray->dox + _t * dDx + dtx * d;
      vec3 dXy = ray->doy + _t * dDy + dty * d;

      intersection->dn[0] = dXx / radius;
      intersection->dn[1] = dXy / radius;
      intersection->duv[0] = dXx;
      intersection->duv[1] = dXy;

      ray->tmax = t;
      return true;
    }
  } else if (delta > 0) {

    // Deux solutions
    float t1 = (-b + sqrtf(delta)) / (2 * a);
    float t2 = (-b - sqrtf(delta)) / (2 * a);
    float t;
    if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax &&
        t2 >= transformedRay.tmin && t2 <= transformedRay.tmax) {
      t = std::min(t1, t2);
    } else if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax) {
      t = t1;
    } else if (t2 >= transformedRay.tmin && t2 <= transformedRay.tmax) {
      t = t2;
    } else {
      return false;
    }
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = mat;
    vec3 objectPoint = transformedRay.orig + (t * transformedRay.dir);
    vec3 objectNormal = normalize(objectPoint - vec3(0, 0, 0));
    vec3 normal = transform.vectorTransformFrom(objectNormal);
    intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
    intersection->normal = normalize(normal);

    float pi = M_PI;
    auto theta = glm::atan(objectPoint.x, objectPoint.z);
    auto vec = glm::vec3(objectPoint.x, objectPoint.y, objectPoint.z);
    auto radius = glm::length(vec);

    auto phi = acos(objectPoint.y / radius);
    auto raw_u = theta / (2 * pi);
    vec3 uv = vec3(1 - (raw_u + 0.5), 1 - phi / pi, 0);

    intersection->u = uv.x;
    intersection->v = uv.y;

    if(!intersection->isOutside){
      objectNormal = -objectNormal;
    }

    vec3 d = normalize(transformedRay.dir);
    float _t = length(t * transformedRay.dir);
    vec3 dDx = ray->ddx;
    vec3 dDy = ray->ddy;

    vec3 dtx =
        -(ray->dox + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
    vec3 dty =
        -(ray->doy + _t * dot(dDy, objectNormal) / dot(d, objectNormal));

    // delta hit point on plane
    vec3 dXx = ray->dox + _t * dDx + dtx * d;
    vec3 dXy = ray->doy + _t * dDy + dty * d;

    intersection->dn[0] = dXx / radius;
    intersection->dn[1] = dXy / radius;
    intersection->duv[0] = dXx;
    intersection->duv[1] = dXy;

    ray->tmax = t;
    return true;
  }
  return false;
}
