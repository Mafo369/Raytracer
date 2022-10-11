#include "triangle.h"

bool Triangle::intersect(Ray *ray, Intersection *intersection) const {
  Ray transformedRay = transformRay(ray);
  vec3 v1v2 = geom.triangle.p2 - geom.triangle.p1;
  vec3 v1v3 = geom.triangle.p3 - geom.triangle.p1;

  vec3 cross_rayDir_v1v3 = cross(transformedRay.dir, v1v3);

  float det = dot(v1v2, cross_rayDir_v1v3);

  if (det > -acne_eps && det < acne_eps)
    return false;

  float inv_det = 1.f / det;

  vec3 o_minus_p1 = transformedRay.orig - geom.triangle.p1;

  float u = dot(o_minus_p1, cross_rayDir_v1v3) * inv_det;

  if (u < 0.f || u > 1.f)
    return false;

  vec3 cross_oMinusp1_v1v2 = cross(o_minus_p1, v1v2);

  float v = dot(transformedRay.dir, cross_oMinusp1_v1v2) * inv_det;

  if (v < 0.f || u + v > 1.f)
    return false;

  float t = dot(v1v3, cross_oMinusp1_v1v2) * inv_det;

  if (t >= transformedRay.tmin && t <= transformedRay.tmax)
  {
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = mat;
    intersection->isOutside = !(det < 0.f);

    vec3 objectNormal = geom.triangle.n2 * u + geom.triangle.n3 * v + geom.triangle.n1 * (1.f - u - v);
    glm::mat4 normalMatrix = glm::transpose(transform.getInvTransform());
    vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
    vec3 normal = vec3(normal4.x, normal4.y, normal4.z);

    auto uv = geom.triangle.tex[1] * u + geom.triangle.tex[2] * v + geom.triangle.tex[0] * (1.f - u - v);
    intersection->u = abs(fmod(uv.x, 1.0));
    intersection->v = abs(fmod(uv.y, 1.0));

    intersection->normal = normalize(normal);
    ray->tmax = t;
    return true;
  }
  return false;
}
