#include "sphere.h"

bool Sphere::intersect(Ray *ray, Intersection *intersection) const {
  Ray transformedRay = transformRay(ray);

  vec3 oc = transformedRay.orig - point3(0.f,0.f,0.f);
  float a = dot(transformedRay.dir, transformedRay.dir);
  float b = 2.f * dot(transformedRay.dir, oc);
  float c = dot(oc, oc) - 1.f;

  float delta = b * b - 4.f * a * c;

  if (delta == 0)
  {
    //Une solution
    float t = -b / (2.f*a);
    if (t >= ray->tmin && t <= ray->tmax)
    {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = mat;
      vec3 objectPoint = invTransform * vec4(intersection->position, 1);
      vec3 objectNormal = objectPoint - vec3(0.f,0.f,0.f);
      glm::mat4 normalMatrix = glm::transpose(invTransform);
      vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
      vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
      intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
      intersection->transform = transform;
      intersection->normal = normalize(normal);

      float pi = M_PI;
      auto theta = acos(-normal.y);
      auto phi = atan2(-normal.z, normal.x) + pi;

      intersection->u = phi / (2*pi);
      intersection->v = theta / pi;

      ray->tmax = t;
      return true;
    }
  }
  else if (delta > 0)
  {

    //Deux solutions
    float t1 = (-b + sqrtf(delta)) / (2*a);
    float t2 = (-b - sqrtf(delta)) / (2*a);
    float t;
    if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax && t2 >= transformedRay.tmin && t2 <= transformedRay.tmax)
    {
      t = std::min(t1, t2);
    }
    else if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax)
    {
      t = t1;
    }
    else if (t2 >= transformedRay.tmin && t2 <= transformedRay.tmax)
    {
      t = t2;
    }
    else
    {
      return false;
    }
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = mat;
    vec3 objectPoint = invTransform * vec4(intersection->position, 1);
    vec3 objectNormal = objectPoint - vec3(0,0,0);
    glm::mat4 normalMatrix = glm::transpose(invTransform);
    vec3 normal = normalMatrix * vec4(objectNormal, 0);
    intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
    intersection->transform = transform;
    intersection->normal = normalize(normal);

    float pi = M_PI;
    auto theta = acos(-normal.y);
    auto phi = atan2(-normal.z, normal.x) + pi;
    intersection->u = phi / (2*pi);
    intersection->v = theta / pi;


    ray->tmax = t;
    return true;
  }
  return false;
}
