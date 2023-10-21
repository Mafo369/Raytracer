#pragma once

typedef struct s_kdtree KdTree;

#include "defines.h"
#include "ray.h"
#include "scene.h"

class Intersection;

bool intersectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection);
KdTree* initKdTree(Scene *scene);
