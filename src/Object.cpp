#include "Object.h"

Object::Object(Material mat, glm::mat4 transform)
  : transform(transform), invTransform(glm::inverse(transform)), mat(mat) {}

Object::~Object() {}

Ray Object::transformRay(Ray* ray) const {
  Ray transformedRay;
  vec4 origin = invTransform * vec4(ray->orig, 1);
  vec4 direction = invTransform * vec4(ray->dir, 0);
  rayInit(&transformedRay, vec3(origin.x, origin.y, origin.z), vec3(direction.x, direction.y, direction.z), ray->tmin, ray->tmax, ray->depth); 
  return transformedRay;
}
