#include "plane.h"

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
  glm::mat4 normalMatrix = glm::transpose(invTransform);
  vec4 normal4 = normalMatrix * vec4(objectNormal, 0.f);
  vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
  intersection->isOutside = true;
  intersection->transform = transform;
  intersection->normal = normalize(normal);
  ray->tmax = t;

  intersection->u = abs(fmod(intersection->position.x, 1.0f));
  intersection->v = abs(fmod(intersection->position.z, 1.0f));

	// Set uv Info
	//vec3 uvw;
	//uvw.x = (1 + intersection.p.x) / 2.f;
	//uvw.y = (1 + intersection.p.y) / 2.f;
	//intersection.uvw = uvw;
	// Ray Differential
	//vec3 duvw[2];
	//duvw[0] = Vec3f(0, 0, 0);
	//duvw[1] = Vec3f(0, 0, 0);

	// dx = dd
	//{
	//	Vec3f d = ray.dir.GetNormalized();
	//	float _t = (t * ray.dir).Length();
	//	Vec3f dDx = (d.Dot(d) * dd_x - d.Dot(dd_x) *	d) / pow(d.Dot(d), 1.5f);
	//	Vec3f dDy = (d.Dot(d) * dd_y - d.Dot(dd_y) *	d) / pow(d.Dot(d), 1.5f);

	//	float dtx = -(0 + _t * dDx.Dot(intersection.N) / d.Dot(intersection.N));
	//	float dty = -(0 + _t * dDy.Dot(intersection.N) / d.Dot(intersection.N));
	//																				  
	//	// delta hit point on plane
	//	Vec3f dXx = 0 +_t* dDx + dtx * d;
	//	Vec3f dXy = 0 +_t* dDy + dty * d;

	//	duvw[0] = dXx / 2.f;
	//	duvw[1] = dXy / 2.f;
	//}

	//intersection.duvw[0] = duvw[0];
	//intersection.duvw[1] = duvw[1];

	return true;
}
