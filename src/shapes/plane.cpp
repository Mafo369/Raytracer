#include "plane.h"
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

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
  vec3 normal = transform.vectorTransformFrom(objectNormal);
  //glm::mat4 normalMatrix = glm::transpose(transform.getInvTransform());
  //vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
  //vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
  intersection->isOutside = true;
  intersection->normal = normalize(normal);

  vec3 uvw = vec3((1.f + x.x) * .5f, (1.f+x.y) * .5f, 0);
  if(mat->m_texture != nullptr){
    vec3 p = uvw - vec3(0,0,0);
    //auto uvt = mat->m_texture->m_transform.transformTo(p);
    intersection->u = uvw.x;
    intersection->v = uvw.y;
    //intersection->duv[0] = vec3(0,0,0);
    //intersection->duv[1] = vec3(0,0,0);

    //float DdotN = dot(ray->dir, intersection->normal);

    //float dtdx = -dot(ray->dox + t * ray->ddx, intersection->normal) / DdotN;
    //float dtdy = -dot(ray->doy + t * ray->ddy, intersection->normal) / DdotN;

    //ray->dox += t * ray->ddx + dtdx * ray->dir; 
    //ray->doy += t * ray->ddy + dtdy * ray->dir; 

    //vec3 doxT = transform.transformTo(intersection->position + (ray->dox * ray->dXPixel));
    //vec3 doyT = transform.transformTo(intersection->position + (ray->doy * ray->dYPixel));

    //vec3 uvwx = vec3((1.f + doxT.x) * .5f, (1.f+doxT.y) * .5f, 0);
    //vec3 uvwy = vec3((1.f + doyT.x) * .5f, (1.f+doyT.y) * .5f, 0);

    //intersection->duv[0] = mat->m_texture->m_transform.transformTo(uvwx+uvw) - uvt;
    //intersection->duv[1] = mat->m_texture->m_transform.transformTo(uvwy+uvw) - uvt;

    //intersection->u = uvwy.x;
    //intersection->v = uvwy.y;

    //vec3 dox = ray->dox + t*ray->ddx;
    //vec3 doy = ray->doy + t*ray->ddy;
    //float dtx = -(t * dot(dox, intersection->normal)) / dot(ray->dir, intersection->normal);
    //float dty = -(t * dot(doy, intersection->normal)) / dot(ray->dir, intersection->normal);
    //assert(dtx == dtx);
    //assert(dty == dty);

    //vec3 new_dox = (dox + dtx * ray->dir);
    //vec3 new_doy = (doy + dty * ray->dir);

    //ray->dox = new_dox;
    //ray->doy = new_doy;

    //vec3 ddx = new_dox * ray->dXPixel;
    //vec3 ddy = new_doy * ray->dYPixel;

    //vec3 px = intersection->position + ddx;
    //vec3 py = intersection->position + ddy;
    ////std::cout <<"inter: " << glm::to_string(intersection->position) << std::endl;
    ////std::cout << glm::to_string(px) << std::endl;

    //vec3 uvwx = vec3(glm::mod(px.x, 1.f), glm::mod(px.y, 1.f), 0.f) - uvt;
    //vec3 uvwy = vec3(glm::mod(py.x, 1.f), glm::mod(py.y, 1.f), 0.f) - uvt;
    ////std::cout << glm::to_string(uvwx) << std::endl;

    //intersection->duv[0] = uvwx;
    //intersection->duv[1] = uvwy;
    //std::cout << glm::to_string(intersection->duv[0]) << std::endl;
    //
    vec3 d = normalize(transformedRay.dir);
		float _t = length(t * transformedRay.dir);
    float det = pow(dot(d,d), 1.5f);
		vec3 dDx = (dot(d,d) * ray->dXPixel - dot(d, ray->dXPixel) *	d) / det;
		vec3 dDy = (dot(d,d) * ray->dYPixel - dot(d, ray->dYPixel) *	d) / det;

		float dtx = -(0 + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
		float dty = -(0 + _t * dot(dDy, objectNormal) / dot(d, objectNormal));
																					  
		// delta hit point on plane
		vec3 dXx = 0.f +_t* dDx + dtx * d;
		vec3 dXy = 0.f +_t* dDy + dty * d;

		intersection->duv[0] = dXx / 2.f;
		intersection->duv[1] = dXy / 2.f;

  }else
  {
    intersection->u = uvw.x;
    intersection->v = uvw.y;
  }


	return true;
}
