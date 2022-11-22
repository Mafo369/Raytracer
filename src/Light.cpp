#include "Light.h"
#include "raytracer.h"
#include "scene.h"

bool Light::is_shadowed( vec3 lightPosition, vec3 normal, vec3 point, Scene* scene, KdTree* tree ) {
    Intersection temp_inter;
    Ray ray;
    vec3 dir = normalize( lightPosition - point );
    // rayInit(&ray, point + (acne_eps * normal), dir, vec2(0,0),0.f, 100);
    rayInit( &ray,
             point + ( acne_eps * normal ),
             dir,
             vec2( 0, 0 ),
             0.f,
             distance( point + ( acne_eps * normal ), lightPosition ) );
    ray.shadow = true;
    ray.dox    = vec3( 0.f );
    ray.doy    = vec3( 0.f );
    ray.ddx    = vec3( 0.f );
    ray.ddy    = vec3( 0.f );
    if ( intersectKdTree( scene, tree, &ray, &temp_inter ) ) { return true; }
    return false;
}

PointLight::PointLight( vec3 position, color3 color, float size ) {
    m_position  = position;
    m_color     = color;
    m_size      = size;
    m_shadowMin = 5;
    m_shadowMax = 10;
    m_samples.push_back( position );
}

PointLight::~PointLight() {}

float PointLight::intensityAt( vec3 point,
                               Scene* scene,
                               KdTree* tree,
                               vec3 view,
                               Intersection* intersection ) {
    // distance square
    // float distanceSquare = glm::length_sq(m_position - point);
    // float distance = sqrt(distanceSquare);
    // float fallFactor = 1.0f / (1.0f + 0.09f * distance + 0.032f * distanceSquare);
    if ( m_size == 0.0 ) {
        return is_shadowed( m_position, intersection->normal, intersection->position, scene, tree )
                   ? 0.f
                   : 1.f /** fallFactor*/;
    }
    else {
        // to detect if we are in the penumbra
        bool penumbra = false;

        // keep track of running shadow variables
        int count;
        float mean = 0.0;

        // calculate random rotation for Halton sequence on our sphere of confusion
        float rotate = uniform01( engine ) * 2.0 * M_PI;

        // cast our minmum number of shadow rays
        for ( count = 0; count < m_shadowMin; count++ ) {

            // calculate (partially) randomized shadow ray
            vec3 p = getLightPoint( point, count, rotate );

            // cast shadow ray
            int val = is_shadowed( p, intersection->normal, intersection->position, scene, tree )
                          ? 0.f
                          : 1.f;

            // update our mean shadow value
            mean = ( (float)( mean * count + val ) ) / ( (float)( count + 1 ) );

            // check if we are in penumbra
            if ( mean != 0.0 && mean != 1.0 ) penumbra = true;
        }

        // continue casting more shadow rays, if in penumbra
        if ( penumbra ) {

            // continue casting shadow rays
            for ( count = m_shadowMin; count < m_shadowMax; count++ ) {

                // calculate (partially) randomized shadow ray
                vec3 p = getLightPoint( point, count, rotate );

                // cast shadow ray
                int val =
                    is_shadowed( p, intersection->normal, intersection->position, scene, tree )
                        ? 0.f
                        : 1.f;

                // update our mean shadow value
                mean = ( (float)( mean * count + val ) ) / ( (float)( count + 1 ) );
            }
        }
        // mean = mean * fallFactor;

        // return our final shaded intensity
        return mean;
    }
}

// Halton sequence generator, with an index (how deep) & a base
float Halton( int index, int base ) {

    // initial value
    float r = 0.0;

    // iterate through (using the base number) to find the value in the sequence
    float f = 1.0 / (float)base;
    for ( int i = index; i > 0; i /= base ) {
        r += f * ( i % base );
        f /= (float)base;
    }

    // return the Halton sequence value
    return r;
}

vec3 PointLight::getLightPoint( point3 p, int c, float r ) {
    vec3 dir     = normalize( m_position - p );
    vec3 v0      = vec3( 0, 1, 0 );
    float dotV0D = dot( v0, dir );
    if ( dotV0D < 0.5 && dotV0D > -0.5 ) v0 = vec3( 0, 0, 1 );
    vec3 v1 = normalize( cross( v0, dir ) );

    float diskRad;
    if ( c < 4 )
        diskRad = 1.0 * m_size;
    else
        diskRad = sqrt( Halton( c - 4, 2 ) ) * m_size;

    float diskRot;
    if ( c == 0 )
        diskRot = 0.0;
    else if ( c == 1 )
        diskRot = 1.0 * M_PI;
    else if ( c == 2 )
        diskRot = 0.5 * M_PI;
    else if ( c == 3 )
        diskRot = 1.5 * M_PI;
    else
        diskRot = Halton( c - 4, 3 ) * 2.0 * M_PI;

    // compute our semi-random position inside the disk
    vec3 pos =
        m_position + ( v0 * diskRad * cos( diskRot + r ) ) + ( v1 * diskRad * sin( diskRot + r ) );
    return pos;
}

vec3 PointLight::getDirection( point3 p ) {
    return p - m_position;
}

AmbientLight::AmbientLight( vec3 position, color3 color ) {
    m_position = position;
    m_color    = color;
    m_samples.push_back( position );
    m_ambient = true;
}

float AmbientLight::intensityAt( vec3 point,
                                 Scene* scene,
                                 KdTree* tree,
                                 vec3 view,
                                 Intersection* intersection ) {
    return 0.0;
}

vec3 AmbientLight::getDirection( point3 p ) {
    return vec3( 0, 0, 0 );
}

AmbientLight::~AmbientLight() {}

AreaLight::AreaLight( vec3 corner,
                      vec3 full_uvec,
                      int usteps,
                      vec3 full_vvec,
                      int vsteps,
                      vec3 color ) {
    m_corner   = corner;
    m_usteps   = usteps;
    m_vsteps   = vsteps;
    uvec       = full_uvec / float( usteps );
    vvec       = full_vvec / float( vsteps );
    m_color    = color;
    nbSamples  = usteps * vsteps;
    m_position = corner + ( full_uvec / 2.f ) + ( full_vvec / 2.f );

    for ( int v = 0; v < m_vsteps; v++ ) {
        for ( int u = 0; u < m_usteps; u++ ) {
            m_samples.push_back( pointOnLight( u, v ) );
        }
    }
}

vec3 AreaLight::getDirection( point3 p ) {
    return vec3( 0, 0, 0 );
}

float AreaLight::intensityAt( vec3 point,
                              Scene* scene,
                              KdTree* tree,
                              vec3 view,
                              Intersection* intersection ) {
    float intensity = 0.0f;
    for ( auto& sample : m_samples ) {
        if ( !is_shadowed( sample, intersection->normal, intersection->position, scene, tree ) ) {
            intensity += 1.0;
        }
    }
    return intensity / float( nbSamples );
}
point3 AreaLight::pointOnLight( float u, float v ) {
    return m_corner + uvec * ( u + uniform01( engine ) ) + vvec * ( v + uniform01( engine ) );
}

AreaLight::~AreaLight() {}

void AreaLight::setup( Scene* scene ) {
    addObject( scene, m_t1 );
    addObject( scene, m_t2 );
}

DirectLight::DirectLight( vec3 direction, color3 color ) : Light() {
    m_position = normalize( direction );
    m_color    = color;
    m_samples.push_back( m_position );
}

float DirectLight::intensityAt( vec3 point,
                                Scene* scene,
                                KdTree* tree,
                                vec3 view,
                                Intersection* intersection ) {
    return is_shadowed( m_position, intersection->normal, intersection->position, scene, tree )
               ? 0.f
               : 1.f;
}

DirectLight::~DirectLight() {}

vec3 DirectLight::getDirection( point3 p ) {
    return m_position;
}

color3
ShapeLight::sample_Li( const Intersection& inter, const point2& u, vec3* wi, float* pdf ) const {
    Intersection pShape = m_shape->sample( inter, u );
    *wi                 = normalize( pShape.position - inter.position );
    *pdf                = m_shape->pdf( inter, *wi );
    return L( pShape, -*wi );
}

float ShapeLight::pdf_Li( const Intersection& it, const vec3& wi ) const {
    return m_shape->pdf( it, wi );
}
