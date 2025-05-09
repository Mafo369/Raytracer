#pragma once
#include "../Object.h"

class Triangle : public Object
{
  public:
    Triangle( size_t materialIndex, Transform transform ) : Object( materialIndex, transform ) {}
    bool intersect( Ray* ray, Intersection* intersection ) const override;
    Intersection sample( const Intersection& inter, const point2& u, float* pdf ) const override;
    Intersection sample( const point2& u ) const override { return Intersection(); }
    float Area() const override;
};
