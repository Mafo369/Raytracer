#include "kdtree.h"
#include "defines.h"
#include "scene.h"
#include "scene_types.h"
#include <stdio.h>

#include <vector>
#include <stack>

#include <limits>

typedef struct s_kdtreeNode KdTreeNode;

struct s_kdtreeNode {
  bool leaf; //! is this node a leaf ?
  int axis;//! axis index of the split, if not leaf
  float split;//!position of the split
  int depth; //!depth in the tree
  std::vector<int> objects;//! index of objects, if leaf
  KdTreeNode* left;//!ptr to left child
  KdTreeNode* right;//! ptr to right child
  vec3 min;//! min pos of node bounding box
  vec3 max;//! max pos of node bounding box
};

KdTreeNode * initNode(bool l, int a, int d) {
    KdTreeNode *ret = new KdTreeNode();
    ret->leaf = l;
    ret->axis = a;
    ret->depth = d;
    ret->left = NULL;
    ret->right = NULL;
    return ret;
}

typedef struct s_stackNode {
    float tmin;
    float tmax;
    KdTreeNode *node;
} StackNode;

struct s_kdtree {
    int depthLimit;
    size_t objLimit;
    KdTreeNode *root;

    std::vector<int> outOfTree;
    std::vector<int> inTree;
};

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node);

KdTree*  initKdTree(Scene *scene) {

  //!\todo compute scene bbox, store object in outOfTree or inTree depending on type
  KdTree *tree = new KdTree();

  for(int i = 0; i < scene->objects.size(); i++){
    if(scene->objects[i]->geom.type == PLANE){
      tree->outOfTree.push_back(i)
    }else if(scene->objects[i]->geom.type == SPHERE){
      tree->inTree.push_back(i);
    }
  }

  tree->depthLimit = 20;
  tree->objLimit = scene->objects.size() * sizeof(Object);

  KdTreeNode *root = initNode(false, 0, 0);
    
  std::vector<float> x_vector;
  std::vector<float> y_vector;
  std::vector<float> z_vector;

  for(int i = 0; i < node->objects.size(); i++){
    float x = scene->objects[i]->geom.position.x;
    float y = scene->objects[i]->geom.position.y;
    float z = scene->objects[i]->geom.position.z;
  
    x_vector.push_back(x);
    y_vector.push_back(y);
    z_vector.push_back(z);   
  }
  float xmin = std::min_element(x_vector.begin(),x_vector.end());
  float ymin = std::min_element(y_vector.begin(),y_vector.end());
  float zmin = std::min_element(z_vector.begin(),z_vector.end());
  float xmax = std::max_element(x_vector.begin(),x_vector.end());
  float ymax = std::max_element(y_vector.begin(),y_vector.end());
  float zmax = std::max_element(z_vector.begin(),z_vector.end());

  root->min = vec3(xmin, ymin, zmin);
  root->max = vec3(xmax, ymax, zmax); 
  
  tree->root = root;
  return tree;

}


//from http://blog.nuclex-games.com/tutorials/collision-detection/static-sphere-vs-aabb/
bool intersectSphereAabb(vec3 sphereCenter, float sphereRadius, vec3 aabbMin, vec3 aabbMax) {
    vec3 closestPointInAabb = min(max(sphereCenter, aabbMin), aabbMax);
    vec3 seg = closestPointInAabb -  sphereCenter;
    float distanceSquared = dot(seg, seg);
    // The AABB and the sphere overlap if the closest point within the rectangle is
    // within the sphere's radius
    return distanceSquared < (sphereRadius * sphereRadius);
}


void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node) {

  //!\todo generate children, compute split position, move objets to children and subdivide if needed.
  KdTreeNode *new_node = new KdTreeNode();
  int d = (node->depth) % 3;
  
  float xmin, ymin, zmin;
  float xmax, ymax, zmax;
  
  vec3 min;
  vec3 max;
  for(int i = 0; i < node->objects.size(); i++){
    float x = node->objects[i]->geom.position.x;
    float y = node->objects[i]->geom.position.y;
    float z = node->objects[i]->geom.position.z;

    if(xmin > x){
      xmin = x;
    }
    if(ymin > y){
      ymin = y;
    }
    if(zmin > z){
      zmin = z;
    }
    if(x > xmax){
      xmax = x;
    }
    if(y > ymax){
      ymax = y;
    }
    if(z > zmax){
      zmax = z;
    }
  }


}

bool traverse(Scene * scene, KdTree * tree, std::stack<StackNode> *stack, StackNode currentNode, Ray * ray, Intersection *intersection) {

    //! \todo traverse kdtree to find intersection

    return false;
}



// from http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-box-intersection/
static bool intersectAabb(Ray *theRay,  vec3 min, vec3 max) {
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    vec3 bounds[2] = {min, max};
    tmin = (bounds[theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
    tmax = (bounds[1-theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
    tymin = (bounds[theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
    tymax = (bounds[1-theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
    if ((tmin > tymax) || (tymin > tmax)) return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;
    tzmin = (bounds[theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
    tzmax = (bounds[1-theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
    if ((tmin > tzmax) || (tzmin > tmax)) return false;
    if (tzmin > tmin) tmin = tzmin;
    if (tzmax < tmax) tmax = tzmax;
    if (tmin > theRay->tmin) theRay->tmin = tmin;
    if (tmax < theRay->tmax) theRay->tmax = tmax;
    return true;
}


bool intersectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection) {
    bool hasIntersection = false;

    //!\todo call vanilla intersection on non kdtree object, then traverse the tree to compute other intersections

    return hasIntersection;
}
