#include "Object.h"

Object::Object(std::shared_ptr<Material> material, Transform transform)
  : transform(transform), mat(material) {}

Object::~Object() {}

Ray Object::transformRay(Ray* ray) const {
  Ray transformedRay;
  vec3 origin = transform.transformTo(ray->orig);
  //vec3 direction = transform.vectorTransformTo(ray->dir);
  vec3 direction = transform.getInvTransform() * vec3(ray->dir);
  rayInit(&transformedRay, origin, direction, ray->pixel,ray->tmin, ray->tmax, ray->depth); 
  return transformedRay;
}
