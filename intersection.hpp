#pragma once

#include "scene.h"
#include "scene_types.h"
#include "raytracer.h"

const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj);

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj);

//Moller-Trumbore Ray Triangle Intersection
bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj);

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection);
