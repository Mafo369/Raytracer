#include "Object.h"

Object::Object( std::shared_ptr<Material> material, Transform transform ) :
    transform( transform ), mat( material ) {}

Object::~Object() {}

Ray Object::transformRay( Ray* ray ) const {
    vec3 origin = transform.transformTo( ray->orig );
    // vec3 direction = transform.vectorTransformTo(ray->dir);
    vec3 direction = transform.getInvTransform() * vec3( ray->dir );
    return Ray(origin, direction, ray->pixel, ray->tmin, ray->tmax, ray->depth );
}

float Object::pdf(const Intersection& inter, const vec3& wi) const {
  auto ray = Ray(inter.position + wi * acne_eps, wi, vec2(0), 0, 10000, 0);
  Intersection isectLight;
  if(!intersect(&ray, &isectLight)) return 0;

  float pdf = glm::length_sq(isectLight.position - inter.position) / (abs(dot(isectLight.normal, -wi) * Area()));
  if(std::isinf(pdf)) pdf = 0.f;
  return pdf;
}
