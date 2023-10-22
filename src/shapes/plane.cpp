#include "plane.h"
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

bool Plane::intersect( Ray* ray, Intersection* intersection ) const {
    vec3 origin = ray->orig;
    vec3 dir = ray->dir;
    transformRay( origin, dir );

    // in obj space
    float rayPz = origin.z;
    float rayDz = dir.z;

    float t = -rayPz / rayDz;
    // std::cout << t <<"   " << - glm::dot(origin, glm::vec3(0,0,1.f)) /
    // (glm::dot(dir, glm::vec3(0,0,1.f)))<< std::endl;
    // hit the opposite face, or not the closest one
    if ( t < 0 || t > ray->tmax ) return false;
    // x is the hit point in the unit plane's plane
    vec3 x = origin + t * dir;
    if ( x.x < -1 || x.x > 1 || x.y < -1 || x.y > 1 ) { return false; }

    // Set hit info
    ray->tmax              = t;
    intersection->position = ray->orig + t * ray->dir;
    intersection->mat      = mat;

    if ( ray->shadow ) return true;

    vec3 objectNormal       = vec3( 0.f, 0.f, 1.f );
    vec3 normal             = transform.vectorTransformFrom( objectNormal );
    intersection->isOutside = true;
    intersection->normal    = normalize( normal );

    vec3 uvw        = vec3( ( 1.f + x.x ) * .5f, ( 1.f + x.y ) * .5f, 0 );
    intersection->u = uvw.x;
    intersection->v = uvw.y;

    if ( ray->hasDifferentials ) {
        auto doxT = transform.transformTo( ray->dox );
        auto doyT = transform.transformTo( ray->doy );
        auto ddxT = transform.getInvTransform() * ray->ddx;
        auto ddyT = transform.getInvTransform() * ray->ddy;

        float tx = -doxT.z / ddxT.z;
        float ty = -doyT.z / ddyT.z;
        vec3 xx  = doxT + tx * ddxT;
        vec3 xy  = doyT + ty * ddyT;

        vec3 uvwx = vec3( ( 1.f + xx.x ) * .5f, ( 1.f + xx.y ) * .5f, 0 );
        vec3 uvwy = vec3( ( 1.f + xy.x ) * .5f, ( 1.f + xy.y ) * .5f, 0 );

        intersection->dn[0]      = vec3( 0 );
        intersection->dn[1]      = vec3( 0 );
        intersection->dpdx       = transform.getTransform() * ( xx - x );
        intersection->dpdy       = transform.getTransform() * ( xy - x );
        intersection->dpdu       = vec3( 0 );
        intersection->dpdv       = vec3( 0 );
        intersection->dudx       = ( uvwx.x - uvw.x );
        intersection->dvdx       = ( uvwx.y - uvw.y );
        intersection->dudy       = ( uvwy.x - uvw.x );
        intersection->dvdy       = ( uvwy.y - uvw.y );
        intersection->parametric = false;
    }

    return true;
}
