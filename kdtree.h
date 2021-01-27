#ifndef __KDTREE_H__
#define __KDTREE_H__
#include "defines.h"
#include "ray.h"
#include "raytracer.h"

typedef struct s_kdtree KdTree;

bool intersectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection);
KdTree*  initKdTree(Scene *scene);
#endif
