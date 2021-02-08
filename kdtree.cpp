#include "kdtree.h"
#include "defines.h"
#include "scene.h"
#include "scene_types.h"
#include <stdio.h>

#include <vector>
#include <stack>

#include <algorithm>
#include <limits>
#include <string.h>

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
static bool intersectAabb(Ray *theRay,  vec3 min, vec3 max);

KdTree*  initKdTree(Scene *scene) {

  //!\todo compute scene bbox, store object in outOfTree or inTree depending on type
  KdTree *tree = new KdTree();

  for(size_t i = 0; i < scene->objects.size(); i++){
    if(scene->objects[i]->geom.type == PLANE){
      tree->outOfTree.push_back(i);
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

  float rad = scene->objects[0]->geom.sphere.radius;

  for(size_t i = 0; i < scene->objects.size(); i++){
    float temp_rad = scene->objects[i]->geom.sphere.radius;
    if(temp_rad > rad){
      rad = temp_rad;
    }
    float x = scene->objects[i]->geom.sphere.center.x;
    float y = scene->objects[i]->geom.sphere.center.y;
    float z = scene->objects[i]->geom.sphere.center.z;

  
    x_vector.push_back(x);
    y_vector.push_back(y);
    z_vector.push_back(z);   
  }
  float xmin = *std::min_element(x_vector.begin(),x_vector.end());
  float ymin = *std::min_element(y_vector.begin(),y_vector.end());
  float zmin = *std::min_element(z_vector.begin(),z_vector.end());
  float xmax = *std::max_element(x_vector.begin(),x_vector.end());
  float ymax = *std::max_element(y_vector.begin(),y_vector.end());
  float zmax = *std::max_element(z_vector.begin(),z_vector.end());

  root->min = vec3(xmin, ymin, zmin)-rad;
  root->max = vec3(xmax, ymax, zmax)+rad; 

  for(size_t i = 0; i < tree->inTree.size();i++){
    root->objects.push_back(tree->inTree[i]);
  }
  
  tree->root = root;
  subdivide(scene, tree, tree->root);
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


float surfaceArea(vec3 min, vec3 max){
  float dx = max.x - min.x;
  float dy = max.y - min.y;
  float dz = max.z - min.z;

  return 2*dx*dy + 2*dx*dz + 2*dy*dz;
}

float prob_hit(vec3 min1, vec3 max1, vec3 min2, vec3 max2){
  return surfaceArea(min2, max2)/surfaceArea(min1, max1); 
}

float lambda(int nl, int nr, float pl, float pr){
  if((nl == 0 || nr == 0) && !(pl == 1 || pr == 1))
    return 0.8f;
  return 1.f;
}

#define COST_TRAVERSE 1.0
#define COST_INTERSECT 1.5

float cost(int nl, int nr, float pl, float pr){
  return (lambda(nl, nr, pl, pr) * (COST_TRAVERSE + COST_INTERSECT * (pl * nl + pr * nr)));
}

void splitBox(int d, vec3 min, vec3 max, float split, vec3 &min_vl, vec3 &max_vl, vec3 &min_vr, vec3 &max_vr){
  min_vl = min;
  max_vl = max;
  min_vr = min;
  max_vr = max;
  if(d == 0){
    max_vl.x = split;
    min_vr.x = split;
  }else if(d == 1){
    max_vl.y = split;
    min_vr.y = split;
  }else{
    max_vl.z = split;
    min_vr.z = split;
  }
}

void sah(int d, float split, vec3 min, vec3 max, int nl, int nr, int np, float &cp){
  cp = INFINITY;

  vec3 min_vl, max_vl;
  vec3 min_vr, max_vr;
  splitBox(d, min, max, split, min_vl, max_vl, min_vr, max_vr);
  float pl, pr;
  pl = prob_hit(min_vl, max_vl, min, max);
  pr = prob_hit(min_vr, max_vr, min, max);
  if(pl == 0 || pr == 0)
    return ;
  
  float cpl, cpr;
  cpl = cost(nl+np, nr, pl, pr);
  cpr = cost(nl, nr+np, pl, pr);

}

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node) {

  //!\todo generate children, compute split position, move objets to children and subdivide if needed.
  
  if(node->depth >= tree->depthLimit || node->objects.size() <= 1){
    node->leaf = true;
    return ;
  }
  
  int d = (node->depth) % 3; // Dimension to split
  KdTreeNode *node_left = initNode(false, d, node->depth + 1);
  KdTreeNode *node_right = initNode(false, d, node->depth + 1);
  float split;

  node_left->min = node->min;
  node_left->max = node->max;
  node_right->min = node->min;
  node_right->max = node->max;
  
  if(d == 0){ // x split
    split = (node->max.x - node->min.x)/2 + node->min.x;
    node_left->max.x = split;
    node_right->min.x = split;
  }else if(d == 1){ // y split
    split = (node->max.y - node->min.y)/2 + node->min.y;
    node_left->max.y = split;
    node_right->min.y = split;
  }else{ // z split
    split = (node->max.z - node->min.z)/2 + node->min.z;
    node_left->max.z = split;
    node_right->min.z = split;
  }

  node->split = split;
  node->left = node_left;
  node->right = node_right;
  
  size_t check = 0;
  for(size_t i = 0; i < node->objects.size(); i++){
    if(intersectSphereAabb(scene->objects[node->objects[i]]->geom.sphere.center, scene->objects[node->objects[i]]->geom.sphere.radius,node_left->min, node_left->max)){
      node_left->objects.push_back(node->objects[i]);
      check++;
    }
    if(intersectSphereAabb(scene->objects[node->objects[i]]->geom.sphere.center, scene->objects[node->objects[i]]->geom.sphere.radius,node_right->min, node_right->max)){
      node_right->objects.push_back(node->objects[i]);
      check++;
    }
    if(check < i){
      printf("wrong check: node.x=%f node.y=%f node.z=%f\n", scene->objects[node->objects[i]]->geom.sphere.center.x, scene->objects[node->objects[i]]->geom.sphere.center.y,
          scene->objects[node->objects[i]]->geom.sphere.center.z);
    }
  }
  if(check < node->objects.size()-1){
    printf("CHECK WRONG: check=%lu size=%lu\n", check, node->objects.size()-1);
    return ;
  }
  node->objects.clear();
   
  subdivide(scene, tree, node_left);
  subdivide(scene, tree, node_right); 
}

bool traverse(Scene * scene, KdTree * tree, std::stack<StackNode> *stack, StackNode currentNode, Ray * ray, Intersection *intersection) {

    //! \todo traverse kdtree to find intersection

    bool hasIntersection = false;
    float dist;

    if(!stack->empty()){
      if(currentNode.node->leaf){ // Leaf -> Iterate through objects and find Intersection (just like in vanilla)
        for(size_t i = 0; i < currentNode.node->objects.size(); i++){
          Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
          if(intersectSphere(ray, temp, scene->objects[currentNode.node->objects[i]])){
            float temp_dist = ray->tmax;
            if(hasIntersection){
              if(temp_dist < dist){
                dist = temp_dist;
                *intersection = *temp;
              }
            }
            else{
              hasIntersection = true;
              *intersection = *temp;
              dist = temp_dist;
            }
          }
          free(temp);
        }
        if(hasIntersection){ // If we find intersection we return true
          return true;
        }
        else{ //Else if no intersection in leaf and stack is not empty, we use stack element as currentNode
          if(!stack->empty()){
            currentNode = stack->top();
            stack->pop();
            return traverse(scene, tree, stack, currentNode, ray, intersection);
          }
        }
      }
      else{ // No leaf -> we traverse both child nodes
        Ray *ray_left = new Ray();
        Ray *ray_right = new Ray();
        rayInit(ray_left, ray->orig, ray->dir);
        rayInit(ray_right, ray->orig, ray->dir);
        bool intersect_left = intersectAabb(ray_left, currentNode.node->left->min, currentNode.node->left->max);
        bool intersect_right = intersectAabb(ray_right, currentNode.node->right->min, currentNode.node->right->max);

        StackNode c_node; //New current node
        if(intersect_left && intersect_right){ // Both rays found intersection
          StackNode s_node; // New stack node
          if(ray_right->tmin < ray_left->tmin){ //Right closer -> right new currentNode and left new stackNode
            s_node.node = currentNode.node->left;
            s_node.tmin = ray_left->tmin;
            s_node.tmax = ray_left->tmax;
            stack->push(s_node);
            
            c_node.node = currentNode.node->right;
            c_node.tmin = ray_right->tmin;
            c_node.tmax = ray_right->tmax;
            currentNode = c_node;
          }
          else{ //Left close -> left new currentNode and right new stackNode
            s_node.node = currentNode.node->right;
            s_node.tmin = ray_right->tmin;
            s_node.tmax = ray_right->tmax;
            stack->push(s_node);
            
            c_node.node = currentNode.node->left;
            c_node.tmin = ray_left->tmin;
            c_node.tmax = ray_left->tmax;
            currentNode = c_node;
          }
        }
        else{ //Only one ray found intersection or no intersection
          if(intersect_left){ // Left new currentNode
            c_node.node = currentNode.node->left;
            c_node.tmin = ray_left->tmin;
            c_node.tmax = ray_left->tmax;
            currentNode = c_node;
          }
          else if(intersect_right){ // Right new currentNode
            c_node.node = currentNode.node->right;
            c_node.tmin = ray_right->tmin;
            c_node.tmax = ray_right->tmax;
            currentNode = c_node;
          }
          else{ // No intersection -> pop stack to use as currentNode
            if(!stack->empty()){
              currentNode = stack->top();
              stack->pop();
            }
          }
        }
        delete ray_left;
        delete ray_right;
      }
      // We traverse the tree with new entries
      return traverse(scene, tree, stack, currentNode, ray, intersection);
    }

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
    float dist;

    Ray *ray_backup = new Ray(); //Ray backup -> we'll use it to find plane intersections
    rayInit(ray_backup, ray->orig, ray->dir);

    if(intersectAabb(ray, tree->root->min, tree->root->max)){ // If ray hits biggest bbox we traverse tree to find sphere intersections
      std::stack<StackNode> stack;                            
      StackNode startNode;
      startNode.node = tree->root;
      startNode.tmin = ray->tmin;
      startNode.tmax = ray->tmax;
      stack.push(startNode);
      hasIntersection = traverse(scene, tree, &stack, startNode, ray, intersection);
    }
    if(hasIntersection){ // If sphere intersection, use tmax as distance reference
      dist = ray->tmax;
    }
    for(size_t i = 0; i < tree->outOfTree.size(); i++){ // Iterate through plane objects to find intersection
      Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
      if(intersectPlane(ray_backup, temp, scene->objects[tree->outOfTree[i]])){
        float temp_dist = ray_backup->tmax;
        if(hasIntersection){
          if(temp_dist < dist){
            ray->tmax = temp_dist;
            ray->tmin = 0;
            dist = temp_dist;
            *intersection = *temp;
          }
        }
        else{
          ray->tmax = temp_dist;
          ray->tmin = 0;
          hasIntersection = true;
          *intersection = *temp;
          dist = temp_dist;
        }
    }
    free(temp);
  }
  delete ray_backup;
  return hasIntersection;
}
