#pragma once

#include "Material.h"
#include <memory>

class Intersection
{
  public:
    Intersection();
    ~Intersection();

    vec3 normal;          //! the normal of the intersection point
    point3 position;      //! the intersection point
    size_t materialIndex; //! the material of th intersected object

    vec2 t;
    Transform transform;

    bool isOutside;
    int face = -1;

    float u;
    float v;

    vec3 dpdv;
    vec3 dn[2];
    vec3 dpdu;
    float dudx = 0;
    float dvdx = 0;
    float dudy = 0;
    float dvdy = 0;
    vec3 dpdx;
    vec3 dpdy;

    bool parametric = true;
    bool hit        = false;

    bool SolveLinearSystem2x2( const float A[2][2], const float B[2], float* x0, float* x1 );

    void computeDifferentials( Ray* ray );

  private:
};
