#pragma once

#include "defines.h"
#include "scene.h"
#include <vector>

typedef struct camera_s {
  point3 position; //! eye position
    vec3 zdir; //! view direction
    vec3 xdir; //! right direction
    vec3 ydir; //! up direction
    point3 center; //! center of the image plane
  float fov;  //! field of view
  float aspect; //! aspect ratio (typically use WIDTH/HEIGHT of the computed image
} Camera;


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

typedef struct object_s {
    glm::mat4 transform;
    glm::mat4 invTransform;
    Geometry geom;
    Material mat;
} Object;

typedef std::vector<Object*> Objects;
typedef std::vector<Light*> Lights;

typedef struct scene_s {
  Lights lights; //! the scene have several lights
  Objects objects; //! the scene have several objects
  Camera cam; //! the scene have one camera
  color3 skyColor; //! the sky color, could be extended to a sky function ;)
} Scene;

