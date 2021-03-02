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

#include <iostream>

#define COST_TRAVERSE 1.0
#define COST_INTERSECT 1.5

#define STARTING 0
#define LYING 1
#define ENDING 2

#define BOTH 3
#define LEFT 4
#define RIGHT 5

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

typedef struct event_t{
  int s;
  float b;
  int k;
  int type;
}Event;

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node);

KdTree*  initKdTree(Scene *scene) {

  //!\todo compute scene bbox, store object in outOfTree or inTree depending on type
  KdTree *tree = new KdTree();

  for(size_t i = 0; i < scene->objects.size(); i++){
    if(scene->objects[i]->geom.type == PLANE){
      tree->outOfTree.push_back(i);
    }else if(scene->objects[i]->geom.type == SPHERE){
      tree->inTree.push_back(i);
    }else if(scene->objects[i]->geom.type == TRIANGLE){
      tree->inTree.push_back(i);
    }
  }

  tree->depthLimit = 30;
  tree->objLimit = scene->objects.size() * sizeof(Object);

  KdTreeNode *root = initNode(false, 0, 0);
    
  std::vector<float> x_vector;
  std::vector<float> y_vector;
  std::vector<float> z_vector;

  for(size_t i = 0; i < tree->inTree.size(); i++){
    if(scene->objects[tree->inTree[i]]->geom.type == SPHERE){
      float rad = scene->objects[tree->inTree[i]]->geom.sphere.radius;

      float x = scene->objects[tree->inTree[i]]->geom.sphere.center.x;
      float y = scene->objects[tree->inTree[i]]->geom.sphere.center.y;
      float z = scene->objects[tree->inTree[i]]->geom.sphere.center.z;

      x_vector.push_back(x+rad);
      x_vector.push_back(x-rad);
      y_vector.push_back(y+rad);
      y_vector.push_back(y-rad);
      z_vector.push_back(z+rad);
      z_vector.push_back(z-rad);
    }else if(scene->objects[tree->inTree[i]]->geom.type == TRIANGLE){
      std::vector<float> px;
      std::vector<float> py;
      std::vector<float> pz;

      px.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p1.x);
      px.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p2.x);
      px.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p3.x);

      py.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p1.y);
      py.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p2.y);
      py.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p3.y);

      pz.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p1.z);
      pz.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p2.z);
      pz.push_back(scene->objects[tree->inTree[i]]->geom.triangle.p3.z);
  
      for(int j = 0; j < 3; j++){
        x_vector.push_back(px[j]);
        y_vector.push_back(py[j]);
        z_vector.push_back(pz[j]); 
      }
    }

  }
  float xmin = *std::min_element(x_vector.begin(),x_vector.end());
  float ymin = *std::min_element(y_vector.begin(),y_vector.end());
  float zmin = *std::min_element(z_vector.begin(),z_vector.end());
  float xmax = *std::max_element(x_vector.begin(),x_vector.end());
  float ymax = *std::max_element(y_vector.begin(),y_vector.end());
  float zmax = *std::max_element(z_vector.begin(),z_vector.end());

  root->min = vec3(xmin, ymin, zmin);
  root->max = vec3(xmax, ymax, zmax); 

  for(size_t i = 0; i < tree->inTree.size();i++){
    root->objects.push_back(tree->inTree[i]);
  }
  
  tree->root = root;
  subdivide(scene, tree, tree->root);
  return tree;
}

#define axistest_x01(a, b, fa, fb)			   \
	p0_ = a*v0.y - b*v0.z;			       	   \
	p2_ = a*v2.y - b*v2.z;			       	   \
        if(p0_<p2_) {min=p0_; max=p2_;} else {min=p2_; max=p0_;} \
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return 0;

#define axistest_x2(a, b, fa, fb)			   \
	p0_ = a*v0.y - b*v0.z;			           \
	p1_ = a*v1.y - b*v1.z;			       	   \
        if(p0_<p1_) {min=p0_; max=p1_;} else {min=p1_; max=p0_;} \
	rad = fa * boxhalfsize.y + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return 0;

#define axistest_y02(a, b, fa, fb)			   \
	p0_ = -a*v0.x + b*v0.z;		      	   \
	p2_ = -a*v2.x + b*v2.z;	       	       	   \
        if(p0_<p2_) {min=p0_; max=p2_;} else {min=p2_; max=p0_;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return 0;

#define axistest_y1(a, b, fa, fb)			   \
	p0_ = -a*v0.x + b*v0.z;		      	   \
	p1_ = -a*v1.x + b*v1.z;	     	       	   \
        if(p0_<p1_) {min=p0_; max=p1_;} else {min=p1_; max=p0_;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.z;   \
	if(min>rad || max<-rad) return 0;

#define axistest_z12(a, b, fa, fb)			   \
	p1_ = a*v1.x - b*v1.y;			           \
	p2_ = a*v2.x - b*v2.y;			       	   \
        if(p2_<p1_) {min=p2_; max=p1_;} else {min=p1_; max=p2_;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   \
	if(min>rad || max<-rad) return 0;

#define axistest_z0(a, b, fa, fb)			   \
	p0_ = a*v0.x - b*v0.y;				   \
	p1_ = a*v1.x - b*v1.y;			           \
        if(p0_<p1_) {min=p0_; max=p1_;} else {min=p1_; max=p0_;} \
	rad = fa * boxhalfsize.x + fb * boxhalfsize.y;   \
	if(min>rad || max<-rad) return 0;

#define findminmax(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

int planeBoxOverlap(vec3 normal, vec3 vert, vec3 maxbox){
  int q;
  vec3 vmin,vmax;
  float v;
  for(q=0;q<=2;q++){
    v=vert[q];					
    if(normal[q]>0.0f){
      vmin[q]=-maxbox[q] - v;	
      vmax[q]= maxbox[q] - v;	
    }else{
      vmin[q]= maxbox[q] - v;	
      vmax[q]=-maxbox[q] - v;	
    }
  }
  if(dot(normal,vmin)>0.0f) return 0;
  if(dot(normal,vmax)>=0.0f) return 1;
  return 0;
}

bool intersectTriangleAabb2(vec3 p1, vec3 p2, vec3 p3, vec3 normal, vec3 aabbMin, vec3 aabbMax){
  float min,max,p0_,p1_,p2_,rad,fex,fey,fez;

  vec3 boxcenter = aabbMin + (aabbMax * 0.5f);
  vec3 boxhalfsize = (aabbMax - aabbMin) * 0.5f;
  vec3 v0 = p1 - boxcenter ;
  vec3 v1 = p2 - boxcenter;
  vec3 v2 = p3 - boxcenter;

  vec3 e0 = v1-v0;
  vec3 e1 = v2-v1;
  vec3 e2 = v0-v2;

  fex = fabsf(e0.x);
  fey = fabsf(e0.y);
  fez = fabsf(e0.z);

  axistest_x01(e0.z, e0.y, fez, fey);
  axistest_y02(e0.z, e0.x, fez, fex);
  axistest_z12(e0.y, e0.x, fey, fex);

  fex = fabsf(e1.x);
  fey = fabsf(e1.y);
  fez = fabsf(e1.z);

  axistest_x01(e1.z, e1.y, fez, fey);
  axistest_y02(e1.z, e1.x, fez, fex);
  axistest_z0(e1.y, e1.x, fey, fex);

  fex = fabsf(e2.x);
  fey = fabsf(e2.y);
  fez = fabsf(e2.z);

  axistest_x2(e2.z, e2.y, fez, fey);

  axistest_y1(e2.z, e2.x, fez, fex);

  axistest_z12(e2.y, e2.x, fey, fex);

  findminmax(v0.x,v1.x,v2.x,min,max);

  if(min>boxhalfsize.x || max<-boxhalfsize.x) return 0;

  findminmax(v0.y,v1.y,v2.y,min,max);

  if(min>boxhalfsize.y || max<-boxhalfsize.y) return 0;

  findminmax(v0.z,v1.z,v2.z,min,max);

  if(min>boxhalfsize.z || max<-boxhalfsize.z) return 0;

  /* Bullet 2: */

   /*  test if the box intersects the plane of the triangle */

   /*  compute plane equation of triangle: normal*x+d=0 */

   vec3 n = cross(e0,e1);

   // -NJMP- (line removed here)

   if(!planeBoxOverlap(n,v0,boxhalfsize)) return 0;	// -NJMP-


   return 1;   /* box and triangle overlaps */

}

//SAT (Seperating Axis Therom)
bool intersectTriangleAabb1(vec3 p1, vec3 p2, vec3 p3, vec3 normal, vec3 aabbMin, vec3 aabbMax){
  vec3 v0 = p1;
  vec3 v1 = p2;
  vec3 v2 = p3;

  if(aabbMin.x <= v0.x && aabbMin.y <= v0.y && aabbMin.z <= v0.z && aabbMax.x >= v0.x && aabbMax.y >= v0.y && aabbMax.z >= v0.z){
    if(aabbMin.x <= v1.x && aabbMin.y <= v1.y && aabbMin.z <= v1.z && aabbMax.x >= v1.x && aabbMax.y >= v1.y && aabbMax.z >= v1.z){
      if(aabbMin.x <= v2.x && aabbMin.y <= v2.y && aabbMin.z <= v2.z && aabbMax.x >= v2.x && aabbMax.y >= v2.y && aabbMax.z >= v2.z){
        return true;
      }
    }
  }

  vec3 c = aabbMin + (aabbMax * 0.5f);
  vec3 e = aabbMax - c;

  v0 -= c;
  v1 -= c;
  v2 -= c;

  vec3 f0 = normalize(v1 -v0);
  vec3 f1 = normalize(v2 - v1);
  vec3 f2 = normalize(v0 - v2);

  vec3 u0 = vec3(1.f, 0.f, 0.f);
  vec3 u1 = vec3(0.f, 1.f ,0.f);
  vec3 u2 = vec3(0.f, 0.f, 1.f);

  vec3 axis_u0_f0 = normalize(cross(u0, f0));
  vec3 axis_u0_f1 = normalize(cross(u0, f1));
  vec3 axis_u0_f2 = normalize(cross(u0, f2));

  vec3 axis_u1_f0 = normalize(cross(u1, f0));
  vec3 axis_u1_f1 = normalize(cross(u1, f1));
  vec3 axis_u1_f2 = normalize(cross(u1, f2));

  vec3 axis_u2_f0 = normalize(cross(u2, f0));
  vec3 axis_u2_f1 = normalize(cross(u2, f1));
  vec3 axis_u2_f2 = normalize(cross(u2, f2));

  // Testing axis: axis_u0_f0
  // Project all 3 vertices of the triangle onto the Seperating axis
  float p0_ = dot(v0, axis_u0_f0);
  float p1_ = dot(v1, axis_u0_f0);
  float p2_ = dot(v2, axis_u0_f0);

  float r = e.x * abs(dot(u0, axis_u0_f0)) +
            e.y * abs(dot(u1, axis_u0_f0)) + 
            e.z * abs(dot(u2, axis_u0_f0));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r){
    printf("u0f0\n");
    return false;
  }

  // Testing axis: axis_u0_f1
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u0_f1);
  p1_ = dot(v1, axis_u0_f1);
  p2_ = dot(v2, axis_u0_f1);

  r = e.x * abs(dot(u0, axis_u0_f1)) +
      e.y * abs(dot(u1, axis_u0_f1)) + 
      e.z * abs(dot(u2, axis_u0_f1));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r){
    printf("u0f1\n");
    return false;
  }

  // Testing axis: axis_u0_f2
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u0_f2);
  p1_ = dot(v1, axis_u0_f2);
  p2_ = dot(v2, axis_u0_f2);

  r = e.x * abs(dot(u0, axis_u0_f2)) +
      e.y * abs(dot(u1, axis_u0_f2)) + 
      e.z * abs(dot(u2, axis_u0_f2));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r){
    printf("u0f2\n");
    return false;
  }

  // Testing axis: axis_u1_f0
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u1_f0);
  p1_ = dot(v1, axis_u1_f0);
  p2_ = dot(v2, axis_u1_f0);

  r = e.x * abs(dot(u0, axis_u1_f0)) +
      e.y * abs(dot(u1, axis_u1_f0)) + 
      e.z * abs(dot(u2, axis_u1_f0));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  // Testing axis: axis_u1_f1
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u1_f1);
  p1_ = dot(v1, axis_u1_f1);
  p2_ = dot(v2, axis_u1_f1);

  r = e.x * abs(dot(u0, axis_u1_f1)) +
      e.y * abs(dot(u1, axis_u1_f1)) + 
      e.z * abs(dot(u2, axis_u1_f1));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  // Testing axis: axis_u1_f2
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u1_f2);
  p1_ = dot(v1, axis_u1_f2);
  p2_ = dot(v2, axis_u1_f2);

  r = e.x * abs(dot(u0, axis_u1_f2)) +
      e.y * abs(dot(u1, axis_u1_f2)) + 
      e.z * abs(dot(u2, axis_u1_f2));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;


  // Testing axis: axis_u2_f0
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u2_f0);
  p1_ = dot(v1, axis_u2_f0);
  p2_ = dot(v2, axis_u2_f0);

  r = e.x * abs(dot(u0, axis_u2_f0)) +
      e.y * abs(dot(u1, axis_u2_f0)) + 
      e.z * abs(dot(u2, axis_u2_f0));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;


  // Testing axis: axis_u2_f1
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u2_f1);
  p1_ = dot(v1, axis_u2_f1);
  p2_ = dot(v2, axis_u2_f1);

  r = e.x * abs(dot(u0, axis_u2_f1)) +
      e.y * abs(dot(u1, axis_u2_f1)) + 
      e.z * abs(dot(u2, axis_u2_f1));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  // Testing axis: axis_u2_f2
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, axis_u2_f2);
  p1_ = dot(v1, axis_u2_f2);
  p2_ = dot(v2, axis_u2_f2);

  r = e.x * abs(dot(u0, axis_u2_f2)) +
      e.y * abs(dot(u1, axis_u2_f2)) + 
      e.z * abs(dot(u2, axis_u2_f2));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  // Testing axis: u0
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, u0);
  p1_ = dot(v1, u0);
  p2_ = dot(v2, u0);

  r = e.x * abs(dot(u0, u0)) +
      e.y * abs(dot(u1, u0)) + 
      e.z * abs(dot(u2, u0));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  // Testing axis: u1
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, u1);
  p1_ = dot(v1, u1);
  p2_ = dot(v2, u1);

  r = e.x * abs(dot(u0, u1)) +
      e.y * abs(dot(u1, u1)) + 
      e.z * abs(dot(u2, u1));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  // Testing axis: u2
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, u2);
  p1_ = dot(v1, u2);
  p2_ = dot(v2, u2);

  r = e.x * abs(dot(u0, u2)) +
      e.y * abs(dot(u1, u2)) + 
      e.z * abs(dot(u2, u2));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  vec3 triangleNormal = cross(f0, f1);
  // Testing axis: triangleNormal
  // Project all 3 vertices of the triangle onto the Seperating axis
  p0_ = dot(v0, triangleNormal);
  p1_ = dot(v1, triangleNormal);
  p2_ = dot(v2, triangleNormal);

  r = e.x * abs(dot(u0, triangleNormal)) +
      e.y * abs(dot(u1, triangleNormal)) + 
      e.z * abs(dot(u2, triangleNormal));

  if(std::max(-std::max({p0_, p1_, p2_}), std::min({p0_, p1_, p2_})) > r)
    return false;

  return true;
}

bool intersectTriangleAabb(vec3 p1, vec3 p2, vec3 p3, vec3 normal, vec3 aabbMin, vec3 aabbMax){
  if(aabbMin.x <= p1.x && aabbMin.y <= p1.y && aabbMin.z <= p1.z && aabbMax.x >= p1.x && aabbMax.y >= p1.y && aabbMax.z >= p1.z)
    return true;
  if(aabbMin.x <= p2.x && aabbMin.y <= p2.y && aabbMin.z <= p2.z && aabbMax.x >= p2.x && aabbMax.y >= p2.y && aabbMax.z >= p2.z)
    return true;
  if(aabbMin.x <= p3.x && aabbMin.y <= p3.y && aabbMin.z <= p3.z && aabbMax.x >= p3.x && aabbMax.y >= p3.y && aabbMax.z >= p3.z)
    return true;
  return false;
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

float p_VSub_V(vec3 min1, vec3 max1, vec3 min2, vec3 max2){
  return surfaceArea(min1, max1)/surfaceArea(min2, max2); 
}

float lambda(int nl, int nr, float pl, float pr){
  if((nl == 0 || nr == 0))
    return 0.8f;
  return 1.f;
}

float cost(int nl, int nr, float pl, float pr){
  return (lambda(nl, nr, pl, pr) * (COST_TRAVERSE + COST_INTERSECT * (pl * nl + pr * nr)));
}

void splitBox(int d, vec3 min, vec3 max, float split, vec3 &min_vl, vec3 &max_vl, vec3 &min_vr, vec3 &max_vr){
  min_vl = min;
  max_vl = max;
  min_vr = min;
  max_vr = max;
  max_vl[d] = split;
  min_vr[d] = split;
}

void sah(vec3 min, vec3 max, float p, int nl, int nr, int np, int k,float &c){
  c = INFINITY;
  vec3 min_vl, max_vl;
  vec3 min_vr, max_vr;

  splitBox(k, min, max, p, min_vl, max_vl, min_vr, max_vr);
  float pl, pr;
  pl = p_VSub_V(min_vl, max_vl, min, max);
  pr = p_VSub_V(min_vr, max_vr, min, max);
  if(pl == 0 || pr == 0)
    return ;

  float cpl, cpr;
  cpl = cost(nl+np, nr, pl, pr);
  cpr = cost(nl, nr+np, pl, pr);

  if (cpl < cpr){
    c = cpl;
  }else{
    c = cpr;
  }
}

void clipSphereToBox(Scene *scene, int sphere, vec3 min, vec3 max, vec3 &minb, vec3 &maxb){
  float radius = scene->objects[sphere]->geom.sphere.radius;
  minb = vec3(scene->objects[sphere]->geom.sphere.center.x - radius, scene->objects[sphere]->geom.sphere.center.y - radius, scene->objects[sphere]->geom.sphere.center.z - radius);
  maxb = vec3(scene->objects[sphere]->geom.sphere.center.x + radius, scene->objects[sphere]->geom.sphere.center.y + radius, scene->objects[sphere]->geom.sphere.center.z + radius);

  for(int k = 0; k < 3; k++){
    if(min[k] > minb[k]){
      minb[k] = min[k];
    }
    if(max[k] < maxb[k]){
      maxb[k] = max[k];
    }
  }
}

void clipTriangleToBox(Scene *scene, int triangle, vec3 min, vec3 max, vec3 &minb, vec3 &maxb){
  vec3 p1 = scene->objects[triangle]->geom.triangle.p1;
  vec3 p2 = scene->objects[triangle]->geom.triangle.p2;
  vec3 p3 = scene->objects[triangle]->geom.triangle.p3;

  float xmin = fmin(p1.x, fmin(p2.x, p3.x));
  float ymin = fmin(p1.y, fmin(p2.y, p3.y));
  float zmin = fmin(p1.z, fmin(p2.z, p3.z));

  float xmax = fmax(p1.x, fmax(p2.x, p3.x));
  float ymax = fmax(p1.y, fmax(p2.y, p3.y));
  float zmax = fmax(p1.z, fmax(p2.z, p3.z));
  
  minb = vec3(xmin, ymin, zmin);
  maxb = vec3(xmax, ymax, zmax);

  for(int k = 0; k < 3; k++){
    if(min[k] > minb[k]){
      minb[k] = min[k];
    }
    if(max[k] < maxb[k]){
      maxb[k] = max[k];
    }
  }
}

bool isPlanar(vec3 min, vec3 max){
  float dx = max.x - min.x;
  float dy = max.y - min.y;
  float dz = max.z - min.z;
  return dx <= 0.01 || dy <= 0.01 || dz <= 0.01;
  //return dx <= 0.f || dy <= 0.f || dz <= 0.f;
}

bool comp_events(Event i, Event j){
  return i.b < j.b;
}

//Incremental sweep to find p
void findPlane(Scene *scene, KdTreeNode *node, float &p_, float &k_, float &c_){
  c_ = INFINITY;
  p_ = 0;

  for(int k = 0; k < 3; k++){
    std::vector<Event> events;
    for(size_t i = 0; i < node->objects.size(); i++){
      vec3 minb, maxb;
      if(scene->objects[node->objects[i]]->geom.type == SPHERE){
        clipSphereToBox(scene, node->objects[i], node->min, node->max, minb, maxb);
      }else{
        clipTriangleToBox(scene, node->objects[i], node->min, node->max, minb, maxb);
      }
      if(isPlanar(minb, maxb)){
        Event e;
        e.s = node->objects[i];
        e.b = minb[k];
        e.k = k;
        e.type = LYING;
        events.push_back(e);
      }else{
        Event e1;
        e1.s = node->objects[i];
        e1.b = minb[k];
        e1.k = k;
        e1.type = STARTING;
        Event e2;
        e2.s = node->objects[i];
        e2.b = maxb[k];
        e2.k = k;
        e2.type = ENDING;
        events.push_back(e1);
        events.push_back(e2);
      }
    }
    sort(events.begin(), events.end(), comp_events);

    int nl = 0, np = 0, nr = node->objects.size();
    for(size_t i = 0; i < events.size(); i++){
      float p = events[i].b;
      int pStarting = 0, pEnding = 0, pLying = 0;

      while(i < events.size() && events[i].b == p && events[i].type == ENDING){
        pEnding++;
        i++;
      }
      while(i < events.size() && events[i].b == p && events[i].type == LYING){
        pLying++;
        i++;
      }
      while(i < events.size() && events[i].b == p && events[i].type == STARTING){
        pStarting++;
        i++;
      }
      np = pLying;
      nr -= pLying;
      nr -= pEnding;

      float c;
      sah(node->min, node->max, p, nl, nr, np, k,c);
      if(c < c_ && !(p <= node->min[k]) && !(p >= node->max[k])){ // New best cost
        c_ = c;
        p_ = p;
        k_ = k;
      }
      nl += pStarting;
      nl += pLying;
      np = 0;
    }
  }
}

bool terminate(float c, size_t n){
  return (c > COST_INTERSECT*n);
}

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node) {
  //printf("SUBDIVIDE\n");
  //!\todo generate children, compute split position, move objets to children and subdivide if needed.
  
  //printf("prev_c=%f size=%lu\n", prev_c, node->objects.size());
  //if(terminate(prev_c, node->objects.size())){
  if(node->depth >= tree->depthLimit || node->objects.size() <= 1){
    node->leaf = true;
    return ;
  }
  
  int d = (node->depth) % 3; // Dimension to split
  KdTreeNode *node_left = initNode(false, d, node->depth + 1);
  KdTreeNode *node_right = initNode(false, d, node->depth + 1);

  node_left->min = node->min;
  node_left->max = node->max;
  node_right->min = node->min;
  node_right->max = node->max;
  
  float p = 0;
  float axis;
  float c = 0;
  findPlane(scene, node, p, axis, c);

  if(c == INFINITY  || c > COST_INTERSECT*node->objects.size()){
    node->leaf = true;
    return ;
  }

  node->split = p;
  node->axis = axis;
  node_left->max[axis] = p;
  node_right->min[axis] = p;

  node->left = node_left;
  node->right = node_right;
  
  size_t check = 0;
  for(size_t i = 0; i < node->objects.size(); i++){
    if(scene->objects[node->objects[i]]->geom.type == SPHERE){
      if(intersectSphereAabb(scene->objects[node->objects[i]]->geom.sphere.center, scene->objects[node->objects[i]]->geom.sphere.radius,node_left->min, node_left->max)){
        node_left->objects.push_back(node->objects[i]);
        check++;
      }
      if(intersectSphereAabb(scene->objects[node->objects[i]]->geom.sphere.center, scene->objects[node->objects[i]]->geom.sphere.radius,node_right->min, node_right->max)){
        node_right->objects.push_back(node->objects[i]);
        check++;
      }
    }
    else if(scene->objects[node->objects[i]]->geom.type == TRIANGLE){
      if(intersectTriangleAabb(scene->objects[node->objects[i]]->geom.triangle.p1, scene->objects[node->objects[i]]->geom.triangle.p2, scene->objects[node->objects[i]]->geom.triangle.p3, scene->objects[node->objects[i]]->geom.triangle.normal,node_left->min, node_left->max)){
        node_left->objects.push_back(node->objects[i]);
        check++;
      }
      if(intersectTriangleAabb(scene->objects[node->objects[i]]->geom.triangle.p1, scene->objects[node->objects[i]]->geom.triangle.p2, scene->objects[node->objects[i]]->geom.triangle.p3, scene->objects[node->objects[i]]->geom.triangle.normal,node_right->min, node_right->max)){
        node_right->objects.push_back(node->objects[i]);
        check++;
      }
    }
    
    if(check < i || check == 0){
      /*printf("wrong check: node.x=%f node.y=%f node.z=%f\n", scene->objects[node->objects[i]]->geom.sphere.center.x, scene->objects[node->objects[i]]->geom.sphere.center.y,
          scene->objects[node->objects[i]]->geom.sphere.center.z);*/
      printf("wrong check: p1(%f, %f, %f) p2(%f, %f , %f) p3(%f, %f, %f)\n", scene->objects[node->objects[i]]->geom.triangle.p1.x, scene->objects[node->objects[i]]->geom.triangle.p1.y,
          scene->objects[node->objects[i]]->geom.triangle.p1.z, scene->objects[node->objects[i]]->geom.triangle.p2.x, scene->objects[node->objects[i]]->geom.triangle.p2.y, scene->objects[node->objects[i]]->geom.triangle.p2.z,
          scene->objects[node->objects[i]]->geom.triangle.p3.x, scene->objects[node->objects[i]]->geom.triangle.p3.y, scene->objects[node->objects[i]]->geom.triangle.p3.z);
      //printf("Wrong check\n");
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

    float t; /*signed distance to the splitting plane */

    StackNode near; //The child of node for half-space containing the origin of R
    StackNode far; //The "other" child of node

    while(!stack->empty()){
      currentNode = stack->top();
      stack->pop();
      while(!currentNode.node->leaf){
        float diff = currentNode.node->right->min[currentNode.node->axis] - ray->orig[currentNode.node->axis];
        t = diff / ray->dir[currentNode.node->axis];

        near.tmin = currentNode.tmin;
        near.tmax = currentNode.tmax;
        far.tmin = currentNode.tmin;
        far.tmax = currentNode.tmax;

        if(diff > 0.0f){
          near.node = currentNode.node->left;
          far.node = currentNode.node->right;
        }else{
          near.node = currentNode.node->right;
          far.node = currentNode.node->left;
        }

        if((t > currentNode.tmax) || (t < 0.0f)){
          currentNode = near;
        }
        else{
          if(t < currentNode.tmin){
            currentNode = far;
          }else{
            far.tmin = t;
            far.tmax = currentNode.tmax;
            stack->push(far);

            near.tmax = t;
            currentNode = near;
          }
        }
      }
      // Leaf found -> Find intersection with objects
      for(size_t i = 0; i < currentNode.node->objects.size(); i++){
        Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
        if(scene->objects[currentNode.node->objects[i]]->geom.type == SPHERE){
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
        }else if(scene->objects[currentNode.node->objects[i]]->geom.type == TRIANGLE){
          if(intersectTriangle(ray, temp, scene->objects[currentNode.node->objects[i]])){
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
        }
        free(temp);
      }
      if(hasIntersection){ // If we find intersection we return true
        return true;
      }
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

  //call vanilla intersection on non kdtree object, then traverse the tree to compute other intersections

  float dist;

  Ray *ray_backup = new Ray(); //Ray backup -> we'll use it to find plane intersections
  rayInit(ray_backup, ray->orig, ray->dir, ray->tmin, ray->tmax);

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
    dist = distance(ray->orig, intersection->position);
    ray_backup->tmax = dist;
  }
  for(size_t i = 0; i < tree->outOfTree.size(); i++){ // Iterate through plane objects to find intersection
    Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
    if(intersectPlane(ray_backup, temp, scene->objects[tree->outOfTree[i]])){
      float temp_dist = ray_backup->tmax;
      if(hasIntersection){
        if(temp_dist < dist){
          dist = temp_dist;
          *intersection = *temp;
          ray->tmax = dist;
        }
      }
      else{
        hasIntersection = true;
        *intersection = *temp;
        dist = temp_dist;
        ray->tmax = dist;
      }
    }
    free(temp);
  }
  delete ray_backup;
  return hasIntersection;
}
