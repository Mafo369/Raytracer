#include "plane.h"
#include <glm/gtc/matrix_transform.hpp>

//bool Plane::intersect(Ray *ray, Intersection *intersection) const {
//  //! \todo : compute intersection of the ray and the plane object
//
//  vec3 n = geom.plane.normal;
//  float dn = dot(ray->dir, n);
//
//  if (dn == 0.0)
//  {
//    //Pas des solutions
//    return false;
//  }
//  else
//  {
//    float t = -(dot(ray->orig, n) + geom.plane.dist) / dn;
//    if (t >= ray->tmin && t <= ray->tmax)
//    {
//      intersection->position = ray->orig + (t * ray->dir);
//      intersection->mat = mat;
//      intersection->normal = n;
//      intersection->u = abs(fmod(intersection->position.x, 1.0f));
//      intersection->v = abs(fmod(intersection->position.z, 1.0f));
//      //if(intersection->position.x < 0 && intersection->position.z > 0){
//      //  intersection->u = 1+fmod(intersection->position.x / 6, 1.0f);
//      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
//      //}
//      //else if(intersection->position.x < 0 && intersection->position.z < 0){
//      //  intersection->u = abs(fmod(intersection->position.x / 6, 1.0f));
//      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
//      //}
//      //else if(intersection->position.x > 0 && intersection->position.z < 0){
//      //  intersection->u = abs(fmod(intersection->position.x / 6, 1.0f));
//      //  intersection->v = 1+fmod(intersection->position.z / 6, 1.0f);
//      //}
//      //else
//      //{
//      //  intersection->u = 1+fmod(intersection->position.x / 6, 1.0f);
//      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
//      //}
//      ray->tmax = t;
//      return true;
//    }
//  }
//
//  return false;
//}

bool Plane::intersect(Ray *ray, Intersection *intersection) const {
  Ray transformedRay = transformRay(ray);

  // in obj space
	float rayPz = transformedRay.orig.z;
  float rayDz = transformedRay.dir.z;

	float t = -rayPz / rayDz;
	// hit the opposite face, or not the closest one
	if (t < 0 || t > ray->tmax) return false;
	// x is the hit point in the unit plane's plane
	vec3 x = transformedRay.orig + t * transformedRay.dir;
	if (x.x < -1 || x.x > 1 || x.y < -1 || x.y > 1) {
		return false;
	}

	// Set hit info
  intersection->position = ray->orig + (t * ray->dir);
  intersection->mat = mat;
  vec3 objectNormal = vec3(0.f,0.f,1.f);
  glm::mat4 normalMatrix = glm::transpose(transform.getInvTransform());
  vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
  vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
  intersection->isOutside = true;
  intersection->normal = normalize(normal);
  ray->tmax = t;

  vec3 uvw = vec3((1.f + x.x) * .5f, (1.f+x.y) * .5f, 0);
  if(mat->m_texture != nullptr){
    vec3 p = uvw - vec3(0,0,0);
    auto uvt = mat->m_texture->m_transform.transformTo(p);
    intersection->u = uvt.x;
    intersection->v = uvt.y;
  
    Ray tdifx = transformRay(ray->difx);
    Ray tdify = transformRay(ray->dify);
    float tdx = -tdifx.orig.z / tdifx.dir.z;
    float tdy = -tdify.orig.z / tdify.dir.z;
    point3 dxi = (tdifx.orig + (tdx * tdifx.dir));
    point3 dyi = (tdify.orig + (tdy * tdify.dir));
    intersection->dxi = dxi;
    intersection->dyi = dyi;

	  //// Ray Differential
    vec3 duvw[2];
    duvw[0] = vec3(0, 0, 0);
    duvw[1] = vec3(0, 0, 0);

    duvw[0].x = 2.0 * (dxi.x - x.x);
    duvw[0].y = 2.0 * (dxi.y - x.y);
    duvw[0].z = 0;

    duvw[1].x = 2.0 * (dyi.x - x.x);
    duvw[1].y = 2.0 * (dyi.y - x.y);
    duvw[1].z = 0;

    intersection->duv[0] = mat->m_texture->m_transform.transformTo(duvw[0] + uvw) - uvt;
    intersection->duv[1] = mat->m_texture->m_transform.transformTo(duvw[1] + uvw) - uvt;
  }else
  {
    intersection->u = uvw.x;
    intersection->v = uvw.y;
  }


	return true;
}
