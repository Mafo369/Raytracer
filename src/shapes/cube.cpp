#include "cube.h"

void Cube::check_axis( float origin, float direction, float& tmin, float& tmax ) const {
    float tmin_numerator = ( -1 - origin );
    float tmax_numerator = ( 1 - origin );

    if ( abs( direction ) >= acne_eps ) {
        tmin = tmin_numerator / direction;
        tmax = tmax_numerator / direction;
    }
    else {
        tmin = tmin_numerator * INFINITY;
        tmax = tmax_numerator * INFINITY;
    }

    if ( tmin > tmax ) std::swap( tmin, tmax );
}

vec3 Cube::computeCubeNormal( vec3 position ) const {
    auto absPoint = glm::abs( position );
    auto maxc     = std::max( std::max( absPoint.x, absPoint.y ), absPoint.z );
    if ( maxc == absPoint.x ) { return normalize( vec3( position.x, 0, 0 ) ); }
    else if ( maxc == absPoint.y ) {
        return normalize( vec3( 0, position.y, 0 ) );
    }
    return normalize( vec3( 0, 0, position.z ) );
}

vec2 Cube::cube_uv_front( point3 point ) const {
    float u = abs( fmod( point.x + 1, 2.0f ) / 2.0f );
    float v = abs( fmod( point.y + 1, 2.0f ) / 2.0f );
    return vec2( u, v );
}

vec2 Cube::cube_uv_back( point3 point ) const {
    float u = abs( fmod( 1 - point.x, 2.0f ) / 2.0f );
    float v = abs( fmod( point.y + 1, 2.0f ) / 2.0f );
    return vec2( u, v );
}

vec2 Cube::cube_uv_left( point3 point ) const {
    float u = abs( fmod( point.z + 1, 2.0f ) / 2.0f );
    float v = abs( fmod( point.y + 1, 2.0f ) / 2.0f );
    return vec2( u, v );
}

vec2 Cube::cube_uv_right( point3 point ) const {
    float u = abs( fmod( 1 - point.z, 2.0f ) / 2.0f );
    float v = abs( fmod( point.y + 1, 2.0f ) / 2.0f );
    return vec2( u, v );
}

vec2 Cube::cube_uv_up( point3 point ) const {
    float u = abs( fmod( point.x + 1, 2.0f ) / 2.0f );
    float v = abs( fmod( 1 - point.z, 2.0f ) / 2.0f );
    return vec2( u, v );
}

vec2 Cube::cube_uv_down( point3 point ) const {
    float u = abs( fmod( point.x + 1, 2.0f ) / 2.0f );
    float v = abs( fmod( point.z + 1, 2.0f ) / 2.0f );
    return vec2( u, v );
}

bool Cube::intersect( Ray* ray, Intersection* intersection ) const {
    vec3 origin = ray->orig;
    vec3 dir = ray->dir;
    transformRay( origin, dir );

    float xtmin, xtmax;
    check_axis( origin.x, dir.x, xtmin, xtmax );
    float ytmin, ytmax;
    check_axis( origin.y, dir.y, ytmin, ytmax );
    float ztmin, ztmax;
    check_axis( origin.z, dir.z, ztmin, ztmax );

    float tmin = std::max( std::max( xtmin, ytmin ), ztmin );
    float tmax = std::min( std::min( xtmax, ytmax ), ztmax );

    if ( tmin < 0 ) tmin = tmax;

    if ( tmin > tmax || tmin < ray->tmin || tmin > ray->tmax ) return false;

    intersection->position = ray->orig + ( tmin * ray->dir );
    intersection->mat      = mat;

    vec3 objectPoint        = transform.transformTo( intersection->position );
    vec3 objectNormal       = computeCubeNormal( objectPoint );
    glm::mat4 normalMatrix  = glm::transpose( transform.getInvTransform() );
    vec3 normal             = normalMatrix * vec4( objectNormal, 1 );
    intersection->isOutside = dot( dir, objectNormal ) < 0;
    intersection->normal    = normalize( normal );

    vec2 uv;
    int face;
    if ( objectNormal.x < 0 ) {
        uv   = cube_uv_left( objectPoint );
        face = 0;
    }
    else if ( objectNormal.x > 0 ) {
        uv   = cube_uv_right( objectPoint );
        face = 1;
    }
    else if ( objectNormal.z > 0 ) {
        uv   = cube_uv_front( objectPoint );
        face = 2;
    }
    else if ( objectNormal.z < 0 ) {
        uv   = cube_uv_back( objectPoint );
        face = 3;
    }
    else if ( objectNormal.y > 0 ) {
        uv   = cube_uv_up( objectPoint );
        face = 4;
    }
    else {
        uv   = cube_uv_down( objectPoint );
        face = 5;
    }
    intersection->u    = uv.x;
    intersection->v    = uv.y;
    intersection->face = face;

    ray->tmax = tmin;
    return true;
}
