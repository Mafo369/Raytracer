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
    //glm::mat4 normalMatrix = glm::transpose(transform.getInvTransform());
    //vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
    //vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
    vec3 normal = transform.vectorTransformFrom(objectNormal);
    intersection->normal = normalize(normal);

    auto uv = geom.triangle.tex[1] * u + geom.triangle.tex[2] * v + geom.triangle.tex[0] * (1.f - u - v);
    intersection->u = abs(fmod(uv.x, 1.0));
    intersection->v = abs(fmod(uv.y, 1.0));
    //std::cout << "uv: " << intersection->u << " " << intersection->v << std::endl;
    if(mat->m_texture != nullptr){
      intersection->dn[0] = vec3(0);
      intersection->dn[1] = vec3(0);

      vec3 d = normalize(transformedRay.dir);
      float _t = length(t * transformedRay.dir);
      vec3 dDx = ray->ddx;
      vec3 dDy = ray->ddy;
      if(dDx.x > 1.f || dDx.y > 1.f || dDx.z > 1.f){
        //std::cout << "ddx: " << glm::to_string(ray->ddx) <<  std::endl;
      }
      if(dDy.x > 1.f || dDy.y > 1.f || dDy.z > 1.f){
        //std::cout << "ddy: " << glm::to_string(ray->ddy) <<  std::endl;
      }

      float dtx = -(0 + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
      float dty = -(0 + _t * dot(dDy, objectNormal) / dot(d, objectNormal));

      //std::cout << "dt" << dtx << " " << dty << std::endl;
                                              
      // delta hit point on plane
      vec3 dXx = 0.f +_t* dDx + dtx * d;
      vec3 dXy = 0.f +_t* dDy + dty * d;

      //ray->dox = x + dXx;
      //ray->doy = x + dXy;
      
      intersection->duv[0] = dXx;
      intersection->duv[1] = dXy;
      //std::cout << "dxx: " << glm::to_string(dXx) <<  std::endl;
      //std::cout << "dxy: " << glm::to_string(dXy) <<  std::endl;
      //std::cout << "TRIANGLE" << std::endl;
      //std::cout << glm::to_string(dXx) << std::endl;
      //std::cout << glm::to_string(dXy) << std::endl;
    }

    ray->tmax = t;
    return true;
  }
  return false;
}
