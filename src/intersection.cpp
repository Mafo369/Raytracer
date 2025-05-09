#include "intersection.h"
#include "ray.h"

Intersection::Intersection() {}

Intersection::~Intersection() {}

bool Intersection::SolveLinearSystem2x2( const float A[2][2],
                                         const float B[2],
                                         float* x0,
                                         float* x1 ) {
    float det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if ( std::abs( det ) < 1e-10f ) return false;
    *x0 = ( A[1][1] * B[0] - A[0][1] * B[1] ) / det;
    *x1 = ( A[0][0] * B[1] - A[1][0] * B[0] ) / det;
    if ( std::isnan( *x0 ) || std::isnan( *x1 ) ) return false;
    return true;
}

void Intersection::computeDifferentials( Ray* ray ) {
    if ( parametric ) {
        float d   = -dot( normal, position );
        float tx  = ( -dot( normal, vec3( ray->dox ) ) - d ) / dot( normal, ray->ddx );
        float ty  = ( -dot( normal, vec3( ray->doy ) ) - d ) / dot( normal, ray->ddy );
        point3 px = ray->dox + tx * ray->ddx;
        point3 py = ray->doy + ty * ray->ddy;

        dpdx = px - position;
        dpdy = py - position;

        int dim[2];
        if ( std::abs( normal.x ) > std::abs( normal.y ) &&
             std::abs( normal.x ) > std::abs( normal.z ) ) {
            dim[0] = 1;
            dim[1] = 2;
        }
        else if ( std::abs( normal.y ) > std::abs( normal.z ) ) {
            dim[0] = 0;
            dim[1] = 2;
        }
        else {
            dim[0] = 0;
            dim[1] = 1;
        }

        float A[2][2] = { { dpdu[dim[0]], dpdv[dim[0]] }, { dpdu[dim[1]], dpdv[dim[1]] } };
        float Bx[2]   = { px[dim[0]] - position[dim[0]], px[dim[1]] - position[dim[1]] };
        float By[2]   = { py[dim[0]] - position[dim[0]], py[dim[1]] - position[dim[1]] };

        if ( !SolveLinearSystem2x2( A, Bx, &dudx, &dvdx ) ) dudx = dvdx = 0;
        if ( !SolveLinearSystem2x2( A, By, &dudy, &dvdy ) ) dudy = dvdy = 0;
    }

    dudx *= 2.f;
    dvdx *= 2.f;
    dudy *= 2.f;
    dvdy *= 2.f;
}
