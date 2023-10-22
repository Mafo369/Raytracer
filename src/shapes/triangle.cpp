#include "triangle.h"
#include "../sampling/sampling.h"

bool Triangle::intersect( Ray* ray, Intersection* intersection ) const {
    vec3 origin = ray->orig;
    vec3 dir = ray->dir;
    transformRay( origin, dir );

    vec3 v1v2          = geom.triangle.p2 - geom.triangle.p1;
    vec3 v1v3          = geom.triangle.p3 - geom.triangle.p1;

    vec3 cross_rayDir_v1v3 = cross( dir, v1v3 );

    float det = dot( v1v2, cross_rayDir_v1v3 );

    if ( det > -acne_eps && det < acne_eps ) return false;

    float inv_det = 1.f / det;

    vec3 o_minus_p1 = origin - geom.triangle.p1;

    float u = dot( o_minus_p1, cross_rayDir_v1v3 ) * inv_det;

    if ( u < 0.f || u > 1.f ) return false;

    vec3 cross_oMinusp1_v1v2 = cross( o_minus_p1, v1v2 );

    float v = dot( dir, cross_oMinusp1_v1v2 ) * inv_det;

    if ( v < 0.f || u + v > 1.f ) return false;

    float t = dot( v1v3, cross_oMinusp1_v1v2 ) * inv_det;

    if ( t >= ray->tmin && t <= ray->tmax ) {
        ray->tmax               = t;
        intersection->position  = ray->orig + ( t * ray->dir );
        intersection->mat       = mat;
        intersection->isOutside = !( det < 0.f );
        if ( ray->shadow ) return true;

        vec3 objectNormal    = normalize( geom.triangle.n2 * u + geom.triangle.n3 * v +
                                       geom.triangle.n1 * ( 1.f - u - v ) );
        vec3 normal          = transform.vectorTransformFrom( objectNormal );
        intersection->normal = normalize( normal );

        auto uv = geom.triangle.tex[1] * u + geom.triangle.tex[2] * v +
                  geom.triangle.tex[0] * ( 1.f - u - v );
        intersection->u = uv.x;
        intersection->v = uv.y;

        // From PBR...
        if ( ray->hasDifferentials ) {
            vec3 dpdu, dpdv;
            vec3 dndu, dndv;
            vec2 duv02 = geom.triangle.tex[0] - geom.triangle.tex[2],
                 duv12 = geom.triangle.tex[1] - geom.triangle.tex[2];
            vec3 dp02  = geom.triangle.p1 - geom.triangle.p3,
                 dp12  = geom.triangle.p2 - geom.triangle.p3;

            float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
            if ( determinant == 0 ) {
                CoordinateSystem( normalize( cross( v1v3, v1v2 ) ), &dpdu, &dpdv );
                dndu = dndv = vec3( 0, 0, 0 );
            }
            else {
                float invdet = 1 / determinant;
                dpdu         = ( duv12[1] * dp02 - duv02[1] * dp12 ) * invdet;
                dpdv         = ( -duv12[0] * dp02 + duv02[0] * dp12 ) * invdet;
                vec3 dn1     = geom.triangle.n1 - geom.triangle.n3;
                vec3 dn2     = geom.triangle.n2 - geom.triangle.n3;
                dndu         = ( duv12[1] * dn1 - duv02[1] * dn2 ) * invdet;
                dndv         = ( -duv12[0] * dn1 + duv02[0] * dn2 ) * invdet;
            }

            intersection->dn[0] = transform.vectorTransformFrom( dndu );
            intersection->dn[1] = transform.vectorTransformFrom( dndv );
            intersection->dpdu  = transform.getTransform() * dpdu;
            intersection->dpdv  = transform.getTransform() * dpdv;
        }

        return true;
    }
    return false;
}

Intersection Triangle::sample( const Intersection& inter, const point2& u, float* pdf ) const {
    vec2 b           = UniformSampleTriangle( u );
    const point3& p0 = geom.triangle.p1;
    const point3& p1 = geom.triangle.p2;
    const point3& p2 = geom.triangle.p3;
    Intersection it;
    it.position = b[0] * p0 + b[1] * p1 + ( 1.f - b[0] - b[1] ) * p2;
    it.normal   = normalize( geom.triangle.n2 * b[0] + geom.triangle.n3 * b[1] +
                           geom.triangle.n1 * ( 1.f - b[0] - b[1] ) );
    *pdf        = 1 / Area();
    return it;
}

float Triangle::Area() const {
    const point3& p0 = geom.triangle.p1;
    const point3& p1 = geom.triangle.p2;
    const point3& p2 = geom.triangle.p3;
    return 0.5 * length( cross( p1 - p0, p2 - p0 ) );
}
