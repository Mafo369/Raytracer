#pragma once

#include "defines.h"

enum Etype {SPHERE=1, PLANE, TRIANGLE, CUBE};

typedef struct geometry_s {
  Etype type; //! what kind of geometry we have, this value allows to determine which part of the union is valid;
    //anonymous union of structures to stores object data
    union {
        struct {
            // sphere
            vec3 center;
            float radius;
        } sphere;
        struct {
            // plan
            vec3 normal;
            float dist;
        } plane;
        struct {
          //triangle
          vec3 p1;
          vec3 p2;
          vec3 p3;
          vec3 n1;
          vec3 n2;
          vec3 n3;
          vec3 normal;
          vec2 tex[3];
        }triangle;
        struct
        {
          vec3 min;
          vec3 max;
        }cube;
        
    };
} Geometry;

