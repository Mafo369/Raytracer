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
  //std::cout << t <<"   " << - glm::dot(transformedRay.orig, glm::vec3(0,0,1.f)) / (glm::dot(transformedRay.dir, glm::vec3(0,0,1.f)))<< std::endl;
	// hit the opposite face, or not the closest one
	if (t < 0 || t > ray->tmax) return false;
	// x is the hit point in the unit plane's plane
	vec3 x = transformedRay.orig + t * transformedRay.dir;
	if (x.x < -1 || x.x > 1 || x.y < -1 || x.y > 1) {
		return false;
	}

	// Set hit info
  ray->tmax = t;
  intersection->position = ray->orig + t * ray->dir;
  intersection->mat = mat;

  vec3 objectNormal = vec3(0.f,0.f,1.f);
  glm::mat4 normalMatrix = glm::transpose(transform.getInvTransform());
  vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
  vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
  intersection->isOutside = true;
  intersection->normal = normalize(normal);

  vec3 uvw = vec3((1.f + x.x) * .5f, (1.f+x.y) * .5f, 0);
  if(mat->m_texture != nullptr){
    vec3 p = uvw - vec3(0,0,0);
    auto uvt = mat->m_texture->m_transform.transformTo(p);
    intersection->u = uvt.x;
    intersection->v = uvt.y;

    vec3 ndir = normalize(transformedRay.dir);
    float st = length(t * transformedRay.dir);
    vec3 ddx = (dot(ndir, ndir) * ray->ddx - dot(ndir, ray->ddx) * ndir) / pow(dot(ndir, ndir), 1.5f);
    vec3 ddy = (dot(ndir, ndir) * ray->ddy - dot(ndir, ray->ddy) * ndir) / pow(dot(ndir, ndir), 1.5f);

    float dtx = -(0 + st * dot(ddx, objectNormal) / dot(ndir, objectNormal));
    float dty = -(0 + st * dot(ddy, objectNormal) / dot(ndir, objectNormal));

    vec3 dpx = 0.f + st * ddx + dtx * ndir;
    vec3 dpy = 0.f + st * ddy + dty * ndir;

    intersection->duv[0] = mat->m_texture->m_transform.transformTo((dpx / 2.f) + uvw) - uvt;
    intersection->duv[1] = mat->m_texture->m_transform.transformTo((dpy / 2.f) + uvw) - uvt;

    //vec3 dox = transform.transformTo(normalize(ray->dox + t*ray->ddx));
    //vec3 doy = transform.transformTo(normalize(ray->doy + t*ray->ddy));
    //float dtx = -dox.z / transformedRay.dir.z;
    //float dty = -doy.z / transformedRay.dir.z;

    //vec3 new_dox = (dox + dtx * transformedRay.dir);
    //vec3 new_doy = (doy + dty * transformedRay.dir);

    //ray->dox = new_dox;
    //ray->doy = new_doy;

    //intersection->duv[0] = mat->m_texture->m_transform.transformTo(new_dox);
    //intersection->duv[1] = mat->m_texture->m_transform.transformTo(new_doy);
  }else
  {
    intersection->u = uvw.x;
    intersection->v = uvw.y;
  }


	return true;
}
