#pragma once
#include "../Object.h"

class Plane : public Object
{
  public:
    Plane( size_t materialIndex, Transform transform ) : Object( materialIndex, transform ) {}
    bool intersect( Ray* ray, Intersection* intersection ) const override;
    Intersection sample( const Intersection& inter, const point2& u, float* pdf ) const override {
        return Intersection();
    }
    Intersection sample( const point2& u ) const override { return Intersection(); }
    float pdf( const Intersection& inter, const vec3& wi ) const override { return 0.f; }
};
