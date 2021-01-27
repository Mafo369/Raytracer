

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "defines.h"
#include "ray.h"
#include "scene.h"
#include "raytracer.h"
#include "image.h"

#include "expected.h"

void validTest(const char *desc, bool value, bool expected){
  printf("%s \t: [%s]\n", desc, value == expected ? "OK":"fail");
}






// from http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-box-intersection/
bool intersectAabb(Ray *theRay,  vec3 min, vec3 max) {
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    vec3 bounds[2] = {min, max};
    tmin = (bounds[theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
    tmax = (bounds[1-theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
    tymin = (bounds[theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
    tymax = (bounds[1-theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    tzmin = (bounds[theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
    tzmax = (bounds[1-theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
    if ((tmin > tzmax) || (tzmin > tmax)) return false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;
    if (tmin > theRay->tmin) theRay->tmin = tmin;
    if (tmax < theRay->tmax) theRay->tmax = tmax;
    return tmin>0 || tmax>0;
}



int main(void){

  Material dummy;
  Object *plane1 = initPlane(vec3(0,0,1), 0, dummy);
  Object *plane2 = initPlane(vec3(1,1,1), 2, dummy);
  Object *sphere1 = initSphere(vec3(0,0,0), 1, dummy);
  Object *sphere2 = initSphere(vec3(1, 1, 1), 0.5, dummy);

  Ray r;
  Intersection dummyInter;

  rayInit(&r, point3(0,0,0), vec3(1,0,0)); validTest("r0 to sphere1", intersectSphere(&r, &dummyInter, sphere1), true);
  rayInit(&r, point3(0,0,0), vec3(1,0,0)); validTest("r0 to sphere2", intersectSphere(&r, &dummyInter, sphere2), false);
  rayInit(&r, point3(0,0,0), vec3(1,0,0)); validTest("r0 to plane1", intersectPlane(&r, &dummyInter, plane1), false);
  rayInit(&r, point3(0,0,0), vec3(1,0,0)); validTest("r0 to plane2", intersectPlane(&r, &dummyInter, plane2), false);

  rayInit(&r, point3(0,0,0), vec3(1,0,0), 0.01, 0.02); validTest("r1 to sphere1", intersectSphere(&r, &dummyInter, sphere1), false);
  rayInit(&r, point3(0,0,0), vec3(1,0,0), 0.01, 0.02); validTest("r1 to sphere2", intersectSphere(&r, &dummyInter, sphere2), false);
  rayInit(&r, point3(0,0,0), vec3(1,0,0), 0.01, 0.02); validTest("r1 to plane1", intersectPlane(&r, &dummyInter, plane1), false);
  rayInit(&r, point3(0,0,0), vec3(1,0,0), 0.01, 0.02); validTest("r1 to plane2", intersectPlane(&r, &dummyInter, plane2), false);

  rayInit(&r, point3(2,2,2), vec3(-1,-1,0-1)); validTest("r2 to sphere1", intersectSphere(&r, &dummyInter, sphere1), true);
  rayInit(&r, point3(2,2,2), vec3(-1,-1,0-1)); validTest("r2 to sphere2", intersectSphere(&r, &dummyInter, sphere2), true);
  rayInit(&r, point3(2,2,2), vec3(-1,-1,0-1)); validTest("r2 to plane1", intersectPlane(&r, &dummyInter, plane1), true);
  rayInit(&r, point3(2,2,2), vec3(-1,-1,0-1)); validTest("r2 to plane2", intersectPlane(&r, &dummyInter, plane2), true);

  freeObject(plane1);
  freeObject(plane2);
  freeObject(sphere1);
  freeObject(sphere2);

  bool beckmann=true;
  for(int i=0; i<beckmannExpectedCount; i++){
    beckmann &= abs(beckmannExpected[i].res - RDM_Beckmann(beckmannExpected[i].NdotH, beckmannExpected[i].alpha))<0.0001f;
  }
  printf("RDM_Beckmann \t: [%s]\n",  beckmann ? "OK":"fail");

  bool fresnel=true;
  for(int i=0; i<fresnelExpectedCount; i++){
    fresnel &= abs(fresnelExpected[i].res - RDM_Fresnel(fresnelExpected[i].LdotH, fresnelExpected[i].extIOR,fresnelExpected[i].intIOR))<0.0001f;
  }
  printf("RDM_Fresnel \t: [%s]\n",  fresnel ? "OK":"fail");


  rayInit(&r, point3(0,-10,0), vec3(0.00,1,0.000));
  printf("intersect aabb : %d\n", intersectAabb(&r, point3(-1, -1, -1), point3(1,1,1)));

  return 0;
}
