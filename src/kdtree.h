#pragma once

#include "defines.h"
#include "ray.h"
#include "scene.h"

class Intersection;

typedef struct s_kdtreeNode KdTreeNode;

struct s_kdtreeNode {
    bool leaf;                //! is this node a leaf ?
    int axis;                 //! axis index of the split, if not leaf
    float split;              //! position of the split
    int depth;                //! depth in the tree
    std::vector<int> objects; //! index of objects, if leaf
    KdTreeNode* left;         //! ptr to left child
    KdTreeNode* right;        //! ptr to right child
    vec3 min;                 //! min pos of node bounding box
    vec3 max;                 //! max pos of node bounding box
};

class KdTree
{
  public:
    int depthLimit;
    size_t objLimit;
    KdTreeNode* root;

    std::vector<int> outOfTree;
    std::vector<int> inTree;
};

bool intersectKdTree( Scene* scene, KdTree* tree, Ray* ray, Intersection* intersection );
KdTree* initKdTree( Scene* scene );
