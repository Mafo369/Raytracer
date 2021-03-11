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

struct s_kdtreeNode
{
  bool leaf;                //! is this node a leaf ?
  int axis;                 //! axis index of the split, if not leaf
  float split;              //!position of the split
  int depth;                //!depth in the tree
  std::vector<int> objects; //! index of objects, if leaf
  KdTreeNode *left;         //!ptr to left child
  KdTreeNode *right;        //! ptr to right child
  vec3 min;                 //! min pos of node bounding box
  vec3 max;                 //! max pos of node bounding box
};

KdTreeNode *initNode(bool l, int a, int d)
{
  KdTreeNode *ret = new KdTreeNode();
  ret->leaf = l;
  ret->axis = a;
  ret->depth = d;
  ret->left = NULL;
  ret->right = NULL;
  return ret;
}

typedef struct s_stackNode
{
  float tmin;
  float tmax;
  KdTreeNode *node;
} StackNode;

struct s_kdtree
{
  int depthLimit;
  size_t objLimit;
  KdTreeNode *root;

  std::vector<int> outOfTree;
  std::vector<int> inTree;
};

typedef struct event_t
{
  int s;
  float b;
  int k;
  int type;
} Event;

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node);

KdTree *initKdTree(Scene *scene)
{

  //!\todo compute scene bbox, store object in outOfTree or inTree depending on type
  KdTree *tree = new KdTree();

  for (size_t i = 0; i < scene->objects.size(); i++)
  {
    if (scene->objects[i]->geom.type == PLANE)
    {
      tree->outOfTree.push_back(i);
    }
    else if (scene->objects[i]->geom.type == SPHERE)
    {
      tree->inTree.push_back(i);
    }
    else if (scene->objects[i]->geom.type == TRIANGLE)
    {
      tree->inTree.push_back(i);
    }
  }

  tree->depthLimit = 40;
  tree->objLimit = scene->objects.size() * sizeof(Object);

  KdTreeNode *root = initNode(false, 0, 0);

  std::vector<float> x_vector;
  std::vector<float> y_vector;
  std::vector<float> z_vector;

  for (size_t i = 0; i < tree->inTree.size(); i++)
  {
    if (scene->objects[tree->inTree[i]]->geom.type == SPHERE)
    {
      float rad = scene->objects[tree->inTree[i]]->geom.sphere.radius;

      float x = scene->objects[tree->inTree[i]]->geom.sphere.center.x;
      float y = scene->objects[tree->inTree[i]]->geom.sphere.center.y;
      float z = scene->objects[tree->inTree[i]]->geom.sphere.center.z;

      x_vector.push_back(x + rad);
      x_vector.push_back(x - rad);
      y_vector.push_back(y + rad);
      y_vector.push_back(y - rad);
      z_vector.push_back(z + rad);
      z_vector.push_back(z - rad);
    }
    else if (scene->objects[tree->inTree[i]]->geom.type == TRIANGLE)
    {
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

      for (int j = 0; j < 3; j++)
      {
        x_vector.push_back(px[j]);
        y_vector.push_back(py[j]);
        z_vector.push_back(pz[j]);
      }
    }
  }
  float xmin = *std::min_element(x_vector.begin(), x_vector.end());
  float ymin = *std::min_element(y_vector.begin(), y_vector.end());
  float zmin = *std::min_element(z_vector.begin(), z_vector.end());
  float xmax = *std::max_element(x_vector.begin(), x_vector.end());
  float ymax = *std::max_element(y_vector.begin(), y_vector.end());
  float zmax = *std::max_element(z_vector.begin(), z_vector.end());

  root->min = vec3(xmin, ymin, zmin);
  root->max = vec3(xmax, ymax, zmax);

  for (size_t i = 0; i < tree->inTree.size(); i++)
  {
    root->objects.push_back(tree->inTree[i]);
  }

  tree->root = root;
  subdivide(scene, tree, tree->root);
  return tree;
}

bool intersectTriangleAabb(vec3 p1, vec3 p2, vec3 p3, vec3 normal, vec3 aabbMin, vec3 aabbMax)
{
  float minx = min(p1.x, min(p2.x, p3.x));
  float miny = min(p1.y, min(p2.y, p3.y));
  float minz = min(p1.z, min(p2.z, p3.z));
  float maxx = max(p1.x, max(p2.x, p3.x));
  float maxy = max(p1.y, max(p2.y, p3.y));
  float maxz = max(p1.z, max(p2.z, p3.z));

  vec3 min = vec3(minx, miny, minz);
  vec3 max = vec3(maxx, maxy, maxz);

  return !((min.x > aabbMax.x) || (max.x < aabbMin.x) || (min.y > aabbMax.y) || (max.y < aabbMin.y) || (min.z > aabbMax.z) || (max.z < aabbMin.z));
}

//from http://blog.nuclex-games.com/tutorials/collision-detection/static-sphere-vs-aabb/
bool intersectSphereAabb(vec3 sphereCenter, float sphereRadius, vec3 aabbMin, vec3 aabbMax)
{
  vec3 closestPointInAabb = min(max(sphereCenter, aabbMin), aabbMax);
  vec3 seg = closestPointInAabb - sphereCenter;
  float distanceSquared = dot(seg, seg);
  // The AABB and the sphere overlap if the closest point within the rectangle is
  // within the sphere's radius
  return distanceSquared < (sphereRadius * sphereRadius);
}

float surfaceArea(vec3 min, vec3 max)
{
  float dx = max.x - min.x;
  float dy = max.y - min.y;
  float dz = max.z - min.z;
  return 2 * dx * dy + 2 * dx * dz + 2 * dy * dz;
}

float p_VSub_V(vec3 min1, vec3 max1, vec3 min2, vec3 max2)
{
  return surfaceArea(min1, max1) / surfaceArea(min2, max2);
}

float lambda(int nl, int nr, float pl, float pr)
{
  if (nl == 0 || nr == 0)
    return 0.8f;
  return 1.f;
}

float cost(int nl, int nr, float pl, float pr)
{
  float Kt = COST_TRAVERSE;
  float Ki = COST_INTERSECT;
  float lambda_p = lambda(nl, nr, pl, pr);
  return (lambda_p * (Kt + Ki * (pl*nl + pr*nr)));
}

void splitBox(int d, vec3 min, vec3 max, float split, vec3 &min_vl, vec3 &max_vl, vec3 &min_vr, vec3 &max_vr)
{
  min_vl = min;
  max_vl = max;
  min_vr = min;
  max_vr = max;
  max_vl[d] = split;
  min_vr[d] = split;
}

void sah(vec3 min, vec3 max, float p, int nl, int nr, int np, int k, float &c)
{
  c = INFINITY;
  vec3 min_vl, max_vl;
  vec3 min_vr, max_vr;

  splitBox(k, min, max, p, min_vl, max_vl, min_vr, max_vr); // Split box in axis k with plane p

  float pl, pr; 
  pl = p_VSub_V(min_vl, max_vl, min, max);
  pr = p_VSub_V(min_vr, max_vr, min, max);
  float cpl, cpr;
  cpl = cost(nl + np, nr, pl, pr);
  cpr = cost(nl, nr + np, pl, pr);
  if (cpl < cpr)
  {
    c = cpl;
  }
  else
  {
    c = cpr;
  }
}

void clipSphereToBox(Scene *scene, int sphere, vec3 min, vec3 max, vec3 &minb, vec3 &maxb)
{
  float radius = scene->objects[sphere]->geom.sphere.radius;
  minb = vec3(scene->objects[sphere]->geom.sphere.center.x - radius, scene->objects[sphere]->geom.sphere.center.y - radius, scene->objects[sphere]->geom.sphere.center.z - radius);
  maxb = vec3(scene->objects[sphere]->geom.sphere.center.x + radius, scene->objects[sphere]->geom.sphere.center.y + radius, scene->objects[sphere]->geom.sphere.center.z + radius);

  for (int k = 0; k < 3; k++)
  {
    if (min[k] > minb[k])
    {
      minb[k] = min[k];
    }
    if (max[k] < maxb[k])
    {
      maxb[k] = max[k];
    }
  }
}

void clipTriangleToBox(Scene *scene, int triangle, vec3 min, vec3 max, vec3 &minb, vec3 &maxb)
{
  vec3 v0 = scene->objects[triangle]->geom.triangle.p1;
  vec3 v1 = scene->objects[triangle]->geom.triangle.p2;
  vec3 v2 = scene->objects[triangle]->geom.triangle.p3;

  const float minx = std::min(v0.x, std::min(v1.x, v2.x));
  const float maxx = std::max(v0.x, std::max(v1.x, v2.x));
  const float miny = std::min(v0.y, std::min(v1.y, v2.y));
  const float maxy = std::max(v0.y, std::max(v1.y, v2.y));
  const float minz = std::min(v0.z, std::min(v1.z, v2.z));
  const float maxz = std::max(v0.z, std::max(v1.z, v2.z));

  minb = vec3(minx, miny, minz);
  maxb = vec3(maxx, maxy, maxz);

  for (int k = 0; k < 3; k++)
  {
    if (min[k] > minb[k])
    {
      minb[k] = min[k];
    }
    if (max[k] < maxb[k])
    {
      maxb[k] = max[k];
    }
  }
}

bool isPlanar(vec3 min, vec3 max)
{
  float dx = max.x - min.x;
  float dy = max.y - min.y;
  float dz = max.z - min.z;
  return dx <= 0.001 || dy <= 0.001 || dz <= 0.001;
}

bool comp_events(Event i, Event j)
{
  return i.b < j.b;
}

//Incremental sweep to find p
void findPlane(Scene *scene, KdTreeNode *node, float &p_, float &k_, float &c_)
{
  c_ = INFINITY;
  p_ = 0;

  for (int k = 0; k < 3; k++)
  {
    std::vector<Event> events;
    for (size_t i = 0; i < node->objects.size(); i++)
    {
      vec3 minb, maxb;
      if (scene->objects[node->objects[i]]->geom.type == SPHERE)
      {
        clipSphereToBox(scene, node->objects[i], node->min, node->max, minb, maxb);
      }
      else
      {
        clipTriangleToBox(scene, node->objects[i], node->min, node->max, minb, maxb);
      }
      if (isPlanar(minb, maxb))
      {
        Event e;
        e.s = node->objects[i];
        e.b = minb[k];
        e.k = k;
        e.type = LYING;
        events.push_back(e);
      }
      else
      {
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
    for (size_t i = 0; i < events.size(); i++)
    {
      float p = events[i].b;
      int pStarting = 0, pEnding = 0, pLying = 0;

      while (i < events.size() && events[i].b == p && events[i].type == ENDING)
      {
        pEnding++;
        i++;
      }
      while (i < events.size() && events[i].b == p && events[i].type == LYING)
      {
        pLying++;
        i++;
      }
      while (i < events.size() && events[i].b == p && events[i].type == STARTING)
      {
        pStarting++;
        i++;
      }
      np = pLying;
      nr -= pLying;
      nr -= pEnding;

      float c;
      sah(node->min, node->max, p, nl, nr, np, k, c);
      if (c < c_ && !(p <= node->min[k]) && !(p >= node->max[k]))
      { // New best cost
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

void subdivide(Scene *scene, KdTree *tree, KdTreeNode *node)
{
  //!\todo generate children, compute split position, move objets to children and subdivide if needed.

  int d = (node->depth) % 3; // Dimension to split (random value to initialize nodes, axis will be chosen in findPlane)
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

  if (c == INFINITY || c > COST_INTERSECT * node->objects.size()) // if Terminate()
  {
    node->leaf = true;
    return;
  }

  node->split = p;
  node->axis = axis;
  node_left->max[axis] = p;
  node_right->min[axis] = p;

  node->left = node_left;
  node->right = node_right;

  for (size_t i = 0; i < node->objects.size(); i++)
  {
    if (scene->objects[node->objects[i]]->geom.type == SPHERE)
    {
      bool left = intersectSphereAabb(scene->objects[node->objects[i]]->geom.sphere.center, scene->objects[node->objects[i]]->geom.sphere.radius, node_left->min, node_left->max);
      bool right = intersectSphereAabb(scene->objects[node->objects[i]]->geom.sphere.center, scene->objects[node->objects[i]]->geom.sphere.radius, node_right->min, node_right->max);
      assert(left || right);
      if (left)
      {
        node_left->objects.push_back(node->objects[i]);
      }
      if (right)
      {
        node_right->objects.push_back(node->objects[i]);
      }
    }
    else if (scene->objects[node->objects[i]]->geom.type == TRIANGLE)
    {
      bool left = intersectTriangleAabb(scene->objects[node->objects[i]]->geom.triangle.p1, scene->objects[node->objects[i]]->geom.triangle.p2, scene->objects[node->objects[i]]->geom.triangle.p3, scene->objects[node->objects[i]]->geom.triangle.normal, node_left->min, node_left->max);
      bool right = intersectTriangleAabb(scene->objects[node->objects[i]]->geom.triangle.p1, scene->objects[node->objects[i]]->geom.triangle.p2, scene->objects[node->objects[i]]->geom.triangle.p3, scene->objects[node->objects[i]]->geom.triangle.normal, node_right->min, node_right->max);
      assert(left || right);
      if (left)
      {
        node_left->objects.push_back(node->objects[i]);
      }
      if (right)
      {
        node_right->objects.push_back(node->objects[i]);
      }
    }
  }
  node->objects.clear();

  subdivide(scene, tree, node_left);
  subdivide(scene, tree, node_right);
}

//Reference [6] : Vlastimil Havran. Heuristic Ray Shooting Algorithms
bool traverse(Scene *scene, KdTree *tree, std::stack<StackNode> *stack, StackNode currentNode, Ray *ray, Intersection *intersection)
{
  //! \todo traverse kdtree to find intersection

  float a, b;
  float t; /*signed distance to the splitting plane */

  a = currentNode.tmin;
  b = currentNode.tmax;

  StackNode near; //The child of node for half-space containing the origin of R
  StackNode far;  //The "other" child of node

  while (!stack->empty())
  {
    currentNode = stack->top();
    stack->pop();
    a = currentNode.tmin;
    b = currentNode.tmax;
    while (!currentNode.node->leaf)
    {
      float diff = currentNode.node->split - ray->orig[currentNode.node->axis];
      t = diff / ray->dir[currentNode.node->axis];

      if (diff > 0.0)
      {
        near.node = currentNode.node->left;
        far.node = currentNode.node->right;
      }
      else
      {
        near.node = currentNode.node->right;
        far.node = currentNode.node->left;
      }

      if ((t > b) || (t < 0.0))
      {
        currentNode = near;
      }
      else
      {
        if (t < a)
        {
          currentNode = far;
        }
        else
        {
          far.tmin = t;
          far.tmax = b;
          stack->push(far);

          currentNode = near;
          b = t;
        }
      }
    }
    ray->tmin = a;
    ray->tmax = b;
    bool hasIntersection = false;
    float dist;
    // Leaf found -> Find intersection with objects in node
    for (size_t i = 0; i < currentNode.node->objects.size(); i++)
    {
      Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
      if (scene->objects[currentNode.node->objects[i]]->geom.type == SPHERE)
      {
        if (intersectSphere(ray, temp, scene->objects[currentNode.node->objects[i]]))
        {
          float temp_dist = ray->tmax;
          if (hasIntersection)
          {
            if (temp_dist < dist)
            {
              dist = temp_dist;
              *intersection = *temp;
            }
          }
          else
          {
            hasIntersection = true;
            *intersection = *temp;
            dist = temp_dist;
          }
        }
      }
      else if (scene->objects[currentNode.node->objects[i]]->geom.type == TRIANGLE)
      {
        if (intersectTriangle(ray, temp, scene->objects[currentNode.node->objects[i]]))
        {
          float temp_dist = ray->tmax;
          if (hasIntersection)
          {
            if (temp_dist < dist)
            {
              dist = temp_dist;
              *intersection = *temp;
            }
          }
          else
          {
            hasIntersection = true;
            *intersection = *temp;
            dist = temp_dist;
          }
        }
      }
      free(temp);
    }
    if (hasIntersection)
    { // If we find intersection we return true
      ray->tmax = dist;
      return true;
    }
  }
  return false;
}

// from http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-box-intersection/
static bool intersectAabb(Ray *theRay, vec3 min, vec3 max)
{
  float tmin, tmax, tymin, tymax, tzmin, tzmax;
  vec3 bounds[2] = {min, max};
  tmin = (bounds[theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
  tmax = (bounds[1 - theRay->sign[0]].x - theRay->orig.x) * theRay->invdir.x;
  tymin = (bounds[theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
  tymax = (bounds[1 - theRay->sign[1]].y - theRay->orig.y) * theRay->invdir.y;
  if ((tmin > tymax) || (tymin > tmax))
    return false;
  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;
  tzmin = (bounds[theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
  tzmax = (bounds[1 - theRay->sign[2]].z - theRay->orig.z) * theRay->invdir.z;
  if ((tmin > tzmax) || (tzmin > tmax))
    return false;
  if (tzmin > tmin)
    tmin = tzmin;
  if (tzmax < tmax)
    tmax = tzmax;
  if (tmin > theRay->tmin)
    theRay->tmin = tmin;
  if (tmax < theRay->tmax)
    theRay->tmax = tmax;
  return true;
}

bool intersectKdTree(Scene *scene, KdTree *tree, Ray *ray, Intersection *intersection)
{
  bool hasIntersection = false;

  //call vanilla intersection on non kdtree object, then traverse the tree to compute other intersections

  float dist;

  Ray *ray_backup = new Ray(); //Ray backup -> we'll use it to find plane intersections
  rayInit(ray_backup, ray->orig, ray->dir, ray->tmin, ray->tmax);

  if (intersectAabb(ray, tree->root->min, tree->root->max))
  { // If ray hits biggest bbox we traverse tree to find sphere intersections
    std::stack<StackNode> stack;
    StackNode startNode;
    startNode.node = tree->root;
    startNode.tmin = ray->tmin;
    startNode.tmax = ray->tmax;
    stack.push(startNode);
    hasIntersection = traverse(scene, tree, &stack, startNode, ray, intersection);
  }
  if (hasIntersection)
  { // If object intersection, use tmax as distance reference
    ray_backup->tmax = ray->tmax;
    dist = ray->tmax;
  }
  for (size_t i = 0; i < tree->outOfTree.size(); i++)
  { // Iterate through plane objects to find intersection
    Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
    if (intersectPlane(ray_backup, temp, scene->objects[tree->outOfTree[i]]))
    {
      float temp_dist = ray_backup->tmax;
      if (hasIntersection)
      {
        if (temp_dist < dist)
        {
          dist = temp_dist;
          *intersection = *temp;
          ray->tmax = dist;
        }
      }
      else
      {
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
