#include "plane.h"
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

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
    intersection->u = uvw.x;
    intersection->v = uvw.y;

    auto doxT = transform.transformTo(ray->dox);
    auto doyT = transform.transformTo(ray->doy);
    auto ddxT = transform.getInvTransform() * ray->ddx;
    auto ddyT = transform.getInvTransform() * ray->ddy;

    float tx = -doxT.z / ddxT.z; 
    float ty = -doyT.z / ddyT.z; 
    vec3 xx = doxT + tx * ddxT;
    vec3 xy = doyT + ty * ddyT;

    vec3 uvwx = vec3((1.f + xx.x) * .5f, (1.f+xx.y) * .5f, 0);
    vec3 uvwy = vec3((1.f + xy.x) * .5f, (1.f+xy.y) * .5f, 0);

    intersection->dn[0] = vec3(0);
    intersection->dn[1] = vec3(0);
    intersection->dpdx = transform.getTransform() * (xx - x);
    intersection->dpdy = transform.getTransform() * (xy - x);
    intersection->dpdu = vec3(0);
    intersection->dpdv = vec3(0);
    intersection->dudx = (uvwx.x - uvw.x);
    intersection->dvdx = (uvwx.y - uvw.y);
    intersection->dudy = (uvwy.x - uvw.x);
    intersection->dvdy = (uvwy.y - uvw.y);
    intersection->parametric = false;

    //vec3 dpdu = (uvwx - uvw);
    //vec3 dpdv = (uvwy - uvw);
    //intersection->dpdu = transform.getTransform() * dpdu;
    //intersection->dpdv = transform.getTransform() * dpdv;

    //vec3 d = normalize(transformedRay.dir);
		//float _t = length(t * transformedRay.dir);

    //vec3 dDx = ray->ddx;
    //vec3 dDy = ray->ddy;

		//vec3 dtx = -(ray->dox + _t * dot(dDx, objectNormal) / dot(d, objectNormal));
		//vec3 dty = -(ray->doy + _t * dot(dDy, objectNormal) / dot(d, objectNormal));
		//																			  
		//// delta hit point on plane
		//vec3 dXx = ray->dox +_t* dDx + dtx * d;
		//vec3 dXy = ray->doy +_t* dDy + dty * d;

    ////ray->dox = x + dXx;
    ////ray->doy = x + dXy;

		////intersection->duv[0] = dXx;
		////intersection->duv[1] = dXy;

  }else
  {
    intersection->u = uvw.x;
    intersection->v = uvw.y;
  }


	return true;
}
