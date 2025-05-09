#include "Object.h"

Object::Object( size_t materialIndex, Transform transform ) :
    transform( transform ), m_MaterialIndex( materialIndex ) {}

Object::~Object() {}

void Object::transformRay( vec3& origin, vec3& direction ) const {
    origin    = transform.transformTo( origin );
    direction = transform.getInvTransform() * vec3( direction );
}

float Object::pdf( const Intersection& inter, const vec3& wi ) const {
    auto ray = Ray( inter.position + wi * acne_eps, wi, 0, 10000, 0 );
    Intersection isectLight;
    if ( !intersect( &ray, &isectLight ) ) return 0;

    float pdf = glm::length_sq( isectLight.position - inter.position ) /
                ( abs( dot( isectLight.normal, -wi ) * Area() ) );
    if ( std::isinf( pdf ) ) pdf = 0.f;
    return pdf;
}
