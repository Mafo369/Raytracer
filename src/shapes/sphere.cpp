#include "sphere.h"
#include "../sampling/sampling.h"

#define PHIMAX ( 2.f * M_PI )
#define THETAMAX M_PI
#define THETAMIN 0.f

void computeDifferentials( vec3 objectPoint,
                           vec3 objectNormal,
                           Ray* ray,
                           float radius,
                           float theta,
                           float* dudx,
                           float* dudy,
                           float* dvdx,
                           float* dvdy,
                           vec3* dndu,
                           vec3* dndv,
                           vec3* dpdu,
                           vec3* dpdv ) {
    float phiMax   = PHIMAX;
    float thetaMax = THETAMAX;
    float thetaMin = THETAMIN;

    // From PBR...

    float zRadius    = std::sqrt( objectPoint.x * objectPoint.x + objectPoint.y * objectPoint.y );
    float invZRadius = 1 / zRadius;
    float cosPhi     = objectPoint.x * invZRadius;
    float sinPhi     = objectPoint.y * invZRadius;
    *dpdu            = vec3( -phiMax * objectPoint.y, phiMax * objectPoint.x, 0 );
    *dpdv            = ( thetaMax - thetaMin ) *
            vec3( objectPoint.z * cosPhi, objectPoint.z * sinPhi, -radius * std::sin( theta ) );

    vec3 d2Pduu = -phiMax * phiMax * vec3( objectPoint.x, objectPoint.y, 0 );
    vec3 d2Pduv = ( thetaMax - thetaMin ) * objectPoint.z * phiMax * vec3( -sinPhi, cosPhi, 0. );
    vec3 d2Pdvv = -( thetaMax - thetaMin ) * ( thetaMax - thetaMin ) *
                  vec3( objectPoint.x, objectPoint.y, objectPoint.z );
    float E = dot( *dpdu, *dpdu );
    float F = dot( *dpdu, *dpdv );
    float G = dot( *dpdv, *dpdv );
    vec3 N  = normalize( cross( *dpdu, *dpdv ) );
    float e = dot( N, d2Pduu );
    float f = dot( N, d2Pduv );
    float g = dot( N, d2Pdvv );

    float invEGF2 = 1 / ( E * G - F * F );
    *dndu = vec3( ( f * F - e * G ) * invEGF2 * *dpdu + ( e * F - f * E ) * invEGF2 * *dpdv );
    *dndv = vec3( ( g * F - f * G ) * invEGF2 * *dpdu + ( f * F - g * E ) * invEGF2 * *dpdv );
}

bool Sphere::intersect( Ray* ray, Intersection* intersection ) const {
    vec3 origin = ray->orig;
    vec3 dir = ray->dir;
    transformRay( origin, dir );

    vec3 oc = origin - point3( 0.f, 0.f, 0.f );
    float a = dot( dir, dir );
    float b = 2.f * dot( dir, oc );
    float c = dot( oc, oc ) - 1.f;

    float delta = b * b - 4.f * a * c;

    if ( delta == 0 ) {
        // Une solution
        float t = -b / ( 2.f * a );
        if ( t >= ray->tmin && t <= ray->tmax ) {
            ray->tmax              = t;
            intersection->position = ray->orig + ( t * ray->dir );
            intersection->mat      = mat;
            if ( ray->shadow ) return true;
            vec3 objectPoint        = transform.transformTo( intersection->position );
            vec3 objectNormal       = objectPoint - vec3( 0.f, 0.f, 0.f );
            vec3 normal             = transform.vectorTransformFrom( objectNormal );
            intersection->isOutside = dot( dir, objectNormal ) < 0;
            intersection->normal    = normalize( normal );

            float pi    = M_PI;
            auto vec    = glm::vec3( objectPoint.x, objectPoint.y, objectPoint.z );
            auto radius = glm::length( vec );

            auto phi = std::atan2( objectPoint.y, objectPoint.x );
            if ( phi < 0 ) phi += 2.f * pi;

            float u     = 1.f - phi / PHIMAX;
            float theta = std::acos( clamp( objectPoint.z / radius, -1.f, 1.f ) );
            float v     = ( theta - THETAMIN ) / ( THETAMAX - THETAMIN );

            intersection->u = u;
            intersection->v = v;

            if ( ray->hasDifferentials ) {
                float dudx, dvdx, dudy, dvdy;
                vec3 dndu, dndv;
                vec3 dpdu, dpdv;
                computeDifferentials( objectPoint,
                                      objectNormal,
                                      ray,
                                      radius,
                                      theta,
                                      &dudx,
                                      &dudy,
                                      &dvdx,
                                      &dvdy,
                                      &dndu,
                                      &dndv,
                                      &dpdu,
                                      &dpdv );
                intersection->dn[0] = transform.vectorTransformFrom( dndu );
                intersection->dn[1] = transform.vectorTransformFrom( dndv );
                intersection->dpdu  = transform.getTransform() * dpdu;
                intersection->dpdv  = transform.getTransform() * dpdv;
            }

            return true;
        }
    }
    else if ( delta > 0 ) {

        // Deux solutions
        float t1 = ( -b + sqrtf( delta ) ) / ( 2.f * a );
        float t2 = ( -b - sqrtf( delta ) ) / ( 2.f * a );
        float t;
        if ( t1 >= ray->tmin && t1 < ray->tmax && t2 >= ray->tmin &&
             t2 < ray->tmax ) {
            t = std::min( t1, t2 );
        }
        else if ( t1 >= ray->tmin && t1 < ray->tmax ) {
            t = t1;
        }
        else if ( t2 >= ray->tmin && t2 < ray->tmax ) {
            t = t2;
        }
        else {
            return false;
        }
        ray->tmax              = t;
        intersection->position = ray->orig + ( t * ray->dir );
        intersection->mat      = mat;
        if ( ray->shadow ) return true;

        vec3 objectPoint        = origin + ( t * dir );
        vec3 objectNormal       = normalize( objectPoint - vec3( 0, 0, 0 ) );
        vec3 normal             = transform.vectorTransformFrom( objectNormal );
        intersection->isOutside = dot( dir, objectNormal ) < 0;
        intersection->normal    = normalize( normal );

        float pi    = M_PI;
        auto vec    = glm::vec3( objectPoint.x, objectPoint.y, objectPoint.z );
        auto radius = glm::length( vec );

        auto phi = std::atan2( objectPoint.y, objectPoint.x );
        if ( phi < 0 ) phi += 2.f * pi;

        float u     = 1.f - phi / PHIMAX;
        float theta = std::acos( clamp( objectPoint.z / radius, -1.f, 1.f ) );
        float v     = ( theta - THETAMIN ) / ( THETAMAX - THETAMIN );

        intersection->u = u;
        intersection->v = v;

        if ( ray->hasDifferentials ) {
            float dudx, dvdx, dudy, dvdy;
            vec3 dndu, dndv;
            vec3 dpdu, dpdv;
            computeDifferentials( objectPoint,
                                  objectNormal,
                                  ray,
                                  radius,
                                  theta,
                                  &dudx,
                                  &dudy,
                                  &dvdx,
                                  &dvdy,
                                  &dndu,
                                  &dndv,
                                  &dpdu,
                                  &dpdv );
            intersection->dn[0] = transform.vectorTransformFrom( dndu );
            intersection->dn[1] = transform.vectorTransformFrom( dndv );
            intersection->dpdu  = transform.getTransform() * dpdu;
            intersection->dpdv  = transform.getTransform() * dpdv;
        }
        return true;
    }
    return false;
}

Intersection Sphere::sample( const point2& u ) const {
    point3 pObj = point3( 0 ) + geom.sphere.radius * UniformSampleSphere( u );
    Intersection it;
    it.normal = normalize( transform.getTransform() * vec3( pObj.x, pObj.y, pObj.z ) );
    pObj *= geom.sphere.radius / length( pObj - point3( 0 ) );
    it.position = transform.transformFrom( pObj );
    return it;
}

Intersection Sphere::sample( const Intersection& inter, const point2& u, float* pdf ) const {
    point3 pCenter = geom.sphere.center;
    vec3 wc        = normalize( pCenter - inter.position );
    vec3 wcX, wcY;
    CoordinateSystem( wc, &wcX, &wcY );

    vec3 normal    = inter.isOutside ? inter.normal : -inter.normal;
    point3 pOrigin = inter.position + ( acne_eps * normal );
    auto dist      = length( pOrigin - pCenter );
    auto radius2   = geom.sphere.radius * geom.sphere.radius;
    if ( ( dist * dist ) <= radius2 ) return sample( u );

    dist               = length( inter.position - pCenter );
    float sinThetaMax2 = radius2 / ( dist * dist );
    float cosThetaMax  = sqrt( max( 0.f, 1.f - sinThetaMax2 ) );
    float cosTheta     = ( 1.f - u[0] ) + u[0] * cosThetaMax;
    float sinTheta     = sqrt( max( 0.f, 1.f - cosTheta * cosTheta ) );
    float phi          = u[1] * 2.f * Pi;

    float ds = dist * cosTheta - sqrt( max( 0.f, radius2 - dist * dist * sinTheta * sinTheta ) );
    float cosAlpha = ( dist * dist + radius2 - ds * ds ) / ( 2.f * dist * geom.sphere.radius );
    float sinAlpha = sqrt( max( 0.f, 1.f - cosAlpha * cosAlpha ) );

    vec3 nObj   = SphericalDirection( sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc );
    point3 pObj = geom.sphere.radius * point3( nObj.x, nObj.y, nObj.z );

    Intersection it;
    pObj *= geom.sphere.radius / length( pObj - point3( 0 ) );

    it.position = transform.transformFrom( pObj );
    it.normal   = transform.getTransform() * nObj;
    *pdf        = 1 / ( 2 * Pi * ( 1 - cosThetaMax ) );
    return it;
}

float Sphere::pdf( const Intersection& inter, const vec3& wi ) const {
    point3 pCenter = transform.transformFrom( point3( 0 ) );
    point3 pOrigin = inter.position + ( acne_eps * inter.normal );
    if ( length_sq( pOrigin - pCenter ) <= geom.sphere.radius * geom.sphere.radius )
        return 1.f; // TODO: return default object pdf

    float sinThetaMax2 =
        geom.sphere.radius * geom.sphere.radius / length_sq( inter.position - pCenter );
    float cosThetaMax = sqrt( max( 0.f, 1.f - sinThetaMax2 ) );
    return UniformConePdf( cosThetaMax );
}
