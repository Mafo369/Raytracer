#include "plane.h"

bool Plane::intersect(Ray *ray, Intersection *intersection) const {
  //! \todo : compute intersection of the ray and the plane object

  vec3 n = geom.plane.normal;
  float dn = dot(ray->dir, n);

  if (dn == 0.0)
  {
    //Pas des solutions
    return false;
  }
  else
  {
    float t = -(dot(ray->orig, n) + geom.plane.dist) / dn;
    if (t >= ray->tmin && t <= ray->tmax)
    {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(mat);
      intersection->normal = n;
      intersection->u = abs(fmod(intersection->position.x, 1.0f));
      intersection->v = abs(fmod(intersection->position.z, 1.0f));
      //if(intersection->position.x < 0 && intersection->position.z > 0){
      //  intersection->u = 1+fmod(intersection->position.x / 6, 1.0f);
      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
      //}
      //else if(intersection->position.x < 0 && intersection->position.z < 0){
      //  intersection->u = abs(fmod(intersection->position.x / 6, 1.0f));
      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
      //}
      //else if(intersection->position.x > 0 && intersection->position.z < 0){
      //  intersection->u = abs(fmod(intersection->position.x / 6, 1.0f));
      //  intersection->v = 1+fmod(intersection->position.z / 6, 1.0f);
      //}
      //else
      //{
      //  intersection->u = 1+fmod(intersection->position.x / 6, 1.0f);
      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
      //}
      ray->tmax = t;
      return true;
    }
  }

  return false;
}
