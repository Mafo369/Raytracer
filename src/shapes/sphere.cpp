#include "sphere.h"

#define PHIMAX (2.f * M_PI)
#define THETAMAX M_PI
#define THETAMIN 0.f

//bool Sphere::intersect(Ray *ray, Intersection *intersection) const {
//  Ray transformedRay = transformRay(ray);
//
//  vec3 oc = transformedRay.orig - point3(0.f, 0.f, 0.f);
//  float a = dot(transformedRay.dir, transformedRay.dir);
//  float b = 2.f * dot(transformedRay.dir, oc);
//  float c = dot(oc, oc) - 1.f;
//
//  float delta = b * b - 4.f * a * c;
//
//  if (delta == 0) {
//    // Une solution
//    float t = -b / (2.f * a);
//    if (t >= ray->tmin && t <= ray->tmax) {
//      intersection->position = ray->orig + (t * ray->dir);
//      intersection->mat = mat;
//      vec3 objectPoint = transform.transformTo(intersection->position);
//      vec3 objectNormal = objectPoint - vec3(0.f, 0.f, 0.f);
//      vec3 normal = transform.vectorTransformFrom(objectNormal);
//      intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
//      intersection->normal = normalize(normal);
//
//      float pi = M_PI;
//      auto theta = glm::atan(objectPoint.x, objectPoint.z);
//      auto vec = glm::vec3(objectPoint.x, objectPoint.y, objectPoint.z);
//      auto radius = glm::length(vec);
//
//      auto phi = acos(objectPoint.y / radius);
//      auto raw_u = theta / (2 * pi);
//      vec3 uv = vec3(1 - (raw_u + 0.5), 1 - phi / pi, 0);
//
//      intersection->u = uv.x;
//      intersection->v = uv.y;
//
//      vec3 d = normalize(transformedRay.dir);
//      float _t = length(t * transformedRay.dir);
//      vec3 dDx = ray->ddx;
//      vec3 dDy = ray->ddy;
//
//      vec3 dtx =
//          -(ray->dox + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
//      vec3 dty =
//          -(ray->doy + _t * dot(dDy, objectNormal) / dot(d, objectNormal));
//
//      // delta hit point on plane
//      vec3 dXx = ray->dox + _t * dDx + dtx * d;
//      vec3 dXy = ray->doy + _t * dDy + dty * d;
//
//      intersection->dn[0] = dXx / radius;
//      intersection->dn[1] = dXy / radius;
//      intersection->duv[0] = dXx;
//      intersection->duv[1] = dXy;
//
//      ray->tmax = t;
//      return true;
//    }
//  } else if (delta > 0) {
//
//    // Deux solutions
//    float t1 = (-b + sqrtf(delta)) / (2 * a);
//    float t2 = (-b - sqrtf(delta)) / (2 * a);
//    float t;
//    if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax &&
//        t2 >= transformedRay.tmin && t2 <= transformedRay.tmax) {
//      t = std::min(t1, t2);
//    } else if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax) {
//      t = t1;
//    } else if (t2 >= transformedRay.tmin && t2 <= transformedRay.tmax) {
//      t = t2;
//    } else {
//      return false;
//    }
//    intersection->position = ray->orig + (t * ray->dir);
//    intersection->mat = mat;
//    vec3 objectPoint = transformedRay.orig + (t * transformedRay.dir);
//    vec3 objectNormal = normalize(objectPoint - vec3(0, 0, 0));
//    vec3 normal = transform.vectorTransformFrom(objectNormal);
//    intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
//    intersection->normal = normalize(normal);
//
//    float pi = M_PI;
//    auto theta = glm::atan(objectPoint.x, objectPoint.z);
//    auto vec = glm::vec3(objectPoint.x, objectPoint.y, objectPoint.z);
//    auto radius = glm::length(vec);
//
//    auto phi = acos(objectPoint.y / radius);
//    auto raw_u = theta / (2 * pi);
//    vec3 uv = vec3(1 - (raw_u + 0.5), 1 - phi / pi, 0);
//
//    intersection->u = uv.x;
//    intersection->v = uv.y;
//
//    if(!intersection->isOutside){
//      objectNormal = -objectNormal;
//    }
//
//    vec3 d = normalize(transformedRay.dir);
//    float _t = length(t * transformedRay.dir);
//    vec3 dDx = ray->ddx;
//    vec3 dDy = ray->ddy;
//
//    vec3 dtx =
//        -(ray->dox + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
//    vec3 dty =
//        -(ray->doy + _t * dot(dDy, objectNormal) / dot(d, objectNormal));
//
//    // delta hit point on plane
//    vec3 dXx = ray->dox + _t * dDx + dtx * d;
//    vec3 dXy = ray->doy + _t * dDy + dty * d;
//
//    intersection->dn[0] = dXx / radius;
//    intersection->dn[1] = dXy / radius;
//    intersection->duv[0] = dXx;
//    intersection->duv[1] = dXy;
//
//    ray->tmax = t;
//    return true;
//  }
//  return false;
//}

//bool SolveLinearSystem2x2(const float A[2][2],
//        const float B[2], float *x0, float *x1) {
//    float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
//    if (std::abs(det) < 1e-10f)
//        return false;
//    *x0 = (A[1][1] * B[0] - A[0][1] * B[1]) / det;
//    *x1 = (A[0][0] * B[1] - A[1][0] * B[0]) / det;
//    if (std::isnan(*x0) || std::isnan(*x1))
//        return false;
//    return true;
//}

void computeDifferentials(vec3 objectPoint, vec3 objectNormal, Ray* ray, float radius, float theta,
                          float* dudx, float* dudy, float* dvdx, float* dvdy, vec3* dndu, vec3* dndv, vec3* dpdu, vec3* dpdv){
    float phiMax = PHIMAX;
    float thetaMax = THETAMAX;
    float thetaMin = THETAMIN;

    float zRadius = std::sqrt(objectPoint.x * objectPoint.x + objectPoint.y * objectPoint.y);
    float invZRadius = 1 / zRadius;
    float cosPhi = objectPoint.x * invZRadius;
    float sinPhi = objectPoint.y * invZRadius;
    *dpdu = vec3(-phiMax * objectPoint.y, phiMax * objectPoint.x, 0);
    *dpdv = (thetaMax - thetaMin) *
             vec3(objectPoint.z * cosPhi, objectPoint.z * sinPhi,
             -radius * std::sin(theta));

    vec3 d2Pduu = -phiMax * phiMax * vec3(objectPoint.x, objectPoint.y, 0);
    vec3 d2Pduv = (thetaMax - thetaMin) * objectPoint.z * phiMax *
                      vec3(-sinPhi, cosPhi, 0.);
    vec3 d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                      vec3(objectPoint.x, objectPoint.y, objectPoint.z);
    float E = dot(*dpdu, *dpdu);
    float F = dot(*dpdu, *dpdv);
    float G = dot(*dpdv, *dpdv);
    vec3 N = normalize(cross(*dpdu, *dpdv));
    float e = dot(N, d2Pduu);
    float f = dot(N, d2Pduv);
    float g = dot(N, d2Pdvv);

    float invEGF2 = 1 / (E * G - F * F);
    *dndu = vec3((f * F - e * G) * invEGF2 * *dpdu + 
                 (e * F - f * E) * invEGF2 * *dpdv);
    *dndv = vec3((g * F - f * G) * invEGF2 * *dpdu + 
                 (f * F - g * E) * invEGF2 * *dpdv);
}

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
      auto vec = glm::vec3(objectPoint.x, objectPoint.y, objectPoint.z);
      auto radius = glm::length(vec);
      
      auto phi = std::atan2(objectPoint.y, objectPoint.x);
      if(phi < 0) phi += 2.f * pi;

      float u = 1.f-phi / PHIMAX;
      float theta = std::acos(clamp(objectPoint.z / radius, -1.f, 1.f ));
      float v = (theta - THETAMIN) / (THETAMAX - THETAMIN);

      intersection->u = u;
      intersection->v = v;

      float dudx, dvdx, dudy, dvdy;
      vec3 dndu, dndv;
      vec3 dpdu, dpdv;
      computeDifferentials(objectPoint, objectNormal, ray, radius, theta, 
                            &dudx, &dudy, &dvdx, &dvdy, &dndu, &dndv, &dpdu, &dpdv);
      intersection->dn[0] = transform.getTransform()* dndu;
      intersection->dn[1] = transform.getTransform()* dndv;
      intersection->dpdu = transform.getTransform() * dpdu;
      intersection->dpdv = transform.getTransform() * dpdv;

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
    auto vec = glm::vec3(objectPoint.x, objectPoint.y, objectPoint.z);
    auto radius = glm::length(vec);
    
    auto phi = std::atan2(objectPoint.y, objectPoint.x);
    if(phi < 0) phi += 2.f * pi;

    float u = 1.f-phi / PHIMAX;
    float theta = std::acos(clamp(objectPoint.z / radius, -1.f, 1.f ));
    float v = (theta - THETAMIN) / (THETAMAX - THETAMIN);

    intersection->u = u;
    intersection->v = v;

    float dudx, dvdx, dudy, dvdy;
    vec3 dndu, dndv;
    vec3 dpdu, dpdv;
    computeDifferentials(objectPoint, objectNormal, ray, radius, theta, 
                          &dudx, &dudy, &dvdx, &dvdy, &dndu, &dndv, &dpdu, &dpdv);
    intersection->dn[0] = transform.vectorTransformFrom(dndu);
    intersection->dn[1] = transform.vectorTransformFrom(dndv);
    intersection->dpdu = transform.getTransform() * dpdu;
    intersection->dpdv = transform.getTransform() * dpdv;

    ray->tmax = t;
    return true;
  }
  return false;
}
