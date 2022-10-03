#include "intersection.hpp"

Ray transform(Ray* ray, glm::mat4 matrix){
  Ray transformedRay;
  vec4 origin = matrix * vec4(ray->orig, 1);
  vec4 direction = matrix * vec4(ray->dir, 0);
  rayInit(&transformedRay, vec3(origin.x, origin.y, origin.z), vec3(direction.x, direction.y, direction.z), ray->tmin, ray->tmax, ray->depth); 
  return transformedRay;
}

bool intersectPlane(Ray *ray, Intersection *intersection, Object *obj)
{

  //! \todo : compute intersection of the ray and the plane object

  vec3 n = obj->geom.plane.normal;
  float dn = dot(ray->dir, n);

  if (dn == 0.0)
  {
    //Pas des solutions
    return false;
  }
  else
  {
    float t = -(dot(ray->orig, n) + obj->geom.plane.dist) / dn;
    if (t >= ray->tmin && t <= ray->tmax)
    {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(obj->mat);
      intersection->normal = n;
      intersection->u = abs(fmod(intersection->position.x, 1.0f));
      intersection->v = abs(fmod(intersection->position.z, 1.0f));
      //if(intersection->position.x < 0 && intersection->position.z > 0){
      //  intersection->u = 1+fmod(intersection->position.x / 6, 1.0f);
      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
      //}
      //else if(intersection->position.x < 0 && intersection->position.z < 0){
      //  intersection->u = abs(fmod(intersection->position.x / 6, 1.0f));
      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
      //}
      //else if(intersection->position.x > 0 && intersection->position.z < 0){
      //  intersection->u = abs(fmod(intersection->position.x / 6, 1.0f));
      //  intersection->v = 1+fmod(intersection->position.z / 6, 1.0f);
      //}
      //else
      //{
      //  intersection->u = 1+fmod(intersection->position.x / 6, 1.0f);
      //  intersection->v = abs(fmod(intersection->position.z / 6, 1.0f));
      //}
      ray->tmax = t;
      return true;
    }
  }

  return false;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *obj)
{
  
  Ray transformedRay = transform(ray, obj->invTransform);

  vec3 oc = transformedRay.orig - point3(0,0,0);
  float a = dot(transformedRay.dir, transformedRay.dir);
  float b = 2.f * dot(transformedRay.dir, oc);
  float c = dot(oc, oc) - 1;

  float delta = b * b - 4 * a * c;

  if (delta == 0)
  {
    //Une solution
    float t = -b / (2*a);
    if (t >= ray->tmin && t <= ray->tmax)
    {
      intersection->position = ray->orig + (t * ray->dir);
      intersection->mat = &(obj->mat);
      vec3 objectPoint = obj->invTransform * vec4(intersection->position, 1);
      vec3 objectNormal = objectPoint - vec3(0,0,0);
      glm::mat4 normalMatrix = glm::transpose(obj->invTransform);
      vec4 normal4 = normalMatrix * vec4(objectNormal, 1);
      vec3 normal = vec3(normal4.x, normal4.y, normal4.z);
      intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
      intersection->transform = obj->transform;
      intersection->normal = normalize(normal);

      float pi = M_PI;
      auto theta = acos(-normal.y);
      auto phi = atan2(-normal.z, normal.x) + pi;

      intersection->u = phi / (2*pi);
      intersection->v = theta / pi;

      ray->tmax = t;
      return true;
    }
  }
  else if (delta > 0)
  {

    //Deux solutions
    float t1 = (-b + sqrtf(delta)) / (2*a);
    float t2 = (-b - sqrtf(delta)) / (2*a);
    float t;
    if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax && t2 >= transformedRay.tmin && t2 <= transformedRay.tmax)
    {
      t = std::min(t1, t2);
    }
    else if (t1 >= transformedRay.tmin && t1 <= transformedRay.tmax)
    {
      t = t1;
    }
    else if (t2 >= transformedRay.tmin && t2 <= transformedRay.tmax)
    {
      t = t2;
    }
    else
    {
      return false;
    }
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = &(obj->mat);
    vec3 objectPoint = obj->invTransform * vec4(intersection->position, 1);
    vec3 objectNormal = objectPoint - vec3(0,0,0);
    glm::mat4 normalMatrix = glm::transpose(obj->invTransform);
    vec3 normal = normalMatrix * vec4(objectNormal, 1);
    intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
    intersection->transform = obj->transform;
    intersection->normal = normalize(normal);

    float pi = M_PI;
    auto theta = acos(-normal.y);
    auto phi = atan2(-normal.z, normal.x) + pi;
    intersection->u = phi / (2*pi);
    intersection->v = theta / pi;


    ray->tmax = t;
    return true;
  }
  return false;
}

//Moller-Trumbore Ray Triangle Intersection
bool intersectTriangle(Ray *ray, Intersection *intersection, Object *obj)
{
  vec3 v1v2 = obj->geom.triangle.p2 - obj->geom.triangle.p1;
  vec3 v1v3 = obj->geom.triangle.p3 - obj->geom.triangle.p1;

  vec3 cross_rayDir_v1v3 = cross(ray->dir, v1v3);

  float det = dot(v1v2, cross_rayDir_v1v3);

  if (det > -acne_eps && det < acne_eps)
    return false;

  float inv_det = 1.f / det;

  vec3 o_minus_p1 = ray->orig - obj->geom.triangle.p1;

  float u = dot(o_minus_p1, cross_rayDir_v1v3) * inv_det;

  if (u < 0.f || u > 1.f)
    return false;

  vec3 cross_oMinusp1_v1v2 = cross(o_minus_p1, v1v2);

  float v = dot(ray->dir, cross_oMinusp1_v1v2) * inv_det;

  if (v < 0.f || u + v > 1.f)
    return false;

  float t = dot(v1v3, cross_oMinusp1_v1v2) * inv_det;

  if (t >= ray->tmin && t <= ray->tmax)
  {
    intersection->position = ray->orig + (t * ray->dir);
    intersection->mat = &(obj->mat);
    intersection->isOutside = dot(ray->dir, obj->geom.triangle.normal) < 0;
    //intersection->normal = intersection->isOutside ? obj->geom.triangle.normal : -obj->geom.triangle.normal;
    auto uv = obj->geom.triangle.tex[1] * u + obj->geom.triangle.tex[2] * v + obj->geom.triangle.tex[0] * (1.f - u - v);
    intersection->u = uv.x;
    intersection->v = uv.y;
    vec3 normal = obj->geom.triangle.n2 * u + obj->geom.triangle.n3 * v + obj->geom.triangle.n1 * (1.f - u - v);

    intersection->normal = normal;
    ray->tmax = t;
    return true;
  }
  return false;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection)
{
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();

  //!\todo loop on each object of the scene to compute intersection

  float dist;

  for (size_t i = 0; i < objectCount; i++)
  {
    Intersection *temp = (Intersection *)malloc(sizeof(Intersection));
    if (scene->objects[i]->geom.type == PLANE)
    {
      if (intersectPlane(ray, temp, scene->objects[i]))
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
    else if (scene->objects[i]->geom.type == SPHERE)
    {
      if (intersectSphere(ray, temp, scene->objects[i]))
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
    else if (scene->objects[i]->geom.type == TRIANGLE)
    {
      if (intersectTriangle(ray, temp, scene->objects[i]))
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
    else if (scene->objects[i]->geom.type == CUBE)
    {
      if (intersectCube(ray, temp, scene->objects[i]))
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
  return hasIntersection;
}

void check_axis(float origin, float direction, float& tmin, float& tmax){
  float tmin_numerator = (-1 - origin);
  float tmax_numerator = (1 - origin);

  if(abs(direction) >= acne_eps){
    tmin = tmin_numerator / direction;
    tmax = tmax_numerator / direction;
  }
  else
  {
    tmin = tmin_numerator * INFINITY;
    tmax = tmax_numerator * INFINITY;
  }

  if(tmin > tmax) std::swap(tmin, tmax);
}

vec3 computeCubeNormal(vec3 position){
    auto absPoint = glm::abs(position);
    auto maxc = std::max(std::max(absPoint.x, absPoint.y), absPoint.z);
    if(maxc == absPoint.x){
      return normalize(vec3(position.x, 0, 0));
    }
    else if(maxc == absPoint.y){
      return normalize(vec3(0, position.y, 0));
    }
    return normalize(vec3(0, 0, position.z)); 
}

vec2 cube_uv_front(point3 point){
  float u = abs(fmod(point.x + 1, 2.0f) / 2.0f);
  float v = abs(fmod(point.y + 1, 2.0f) / 2.0f);
  return vec2(u, v);
}

vec2 cube_uv_back(point3 point){
  float u = abs(fmod(1 - point.x, 2.0f) / 2.0f);
  float v = abs(fmod(point.y + 1, 2.0f) / 2.0f);
  return vec2(u, v);
}

vec2 cube_uv_left(point3 point){
  float u = abs(fmod(point.z + 1, 2.0f) / 2.0f);
  float v = abs(fmod(point.y + 1, 2.0f) / 2.0f);
  return vec2(u, v);
}

vec2 cube_uv_right(point3 point){
  float u = abs(fmod(1 - point.z, 2.0f) / 2.0f);
  float v = abs(fmod(point.y + 1, 2.0f) / 2.0f);
  return vec2(u, v);
}

vec2 cube_uv_up(point3 point){
  float u = abs(fmod(point.x + 1, 2.0f) / 2.0f);
  float v = abs(fmod(1 - point.z, 2.0f) / 2.0f);
  return vec2(u, v);
}

vec2 cube_uv_down(point3 point){
  float u = abs(fmod(point.x + 1, 2.0f) / 2.0f);
  float v = abs(fmod(point.z + 1, 2.0f) / 2.0f);
  return vec2(u, v);
}

bool intersectCube(Ray* ray, Intersection* intersection, Object* obj){
    
    Ray transformedRay = transform(ray, obj->invTransform);

    float xtmin, xtmax;
    check_axis(transformedRay.orig.x, transformedRay.dir.x, xtmin, xtmax);
    float ytmin, ytmax;
    check_axis(transformedRay.orig.y, transformedRay.dir.y, ytmin, ytmax);
    float ztmin, ztmax;
    check_axis(transformedRay.orig.z, transformedRay.dir.z, ztmin, ztmax);

    float tmin = std::max(std::max(xtmin, ytmin), ztmin);
    float tmax = std::min(std::min(xtmax, ytmax), ztmax);

    if(tmin > tmax || tmin < ray->tmin || tmin > ray->tmax)
      return false;

    intersection->position = ray->orig + (tmin * ray->dir);
    intersection->mat = &(obj->mat);

    vec3 objectPoint = obj->invTransform * vec4(intersection->position, 1);
    vec3 objectNormal = computeCubeNormal(objectPoint);
    glm::mat4 normalMatrix = glm::transpose(obj->invTransform);
    vec3 normal = normalMatrix * vec4(objectNormal, 1);
    intersection->isOutside = dot(transformedRay.dir, objectNormal) < 0;
    intersection->transform = obj->transform;
    intersection->normal = normalize(normal);

    vec2 uv;
    int face;
    if(objectNormal.x < 0){
      uv = cube_uv_left(objectPoint);
      face = 0;
    }
    else if(objectNormal.x > 0){
      uv = cube_uv_right(objectPoint);
      face = 1;
    }
    else if(objectNormal.z > 0) {
      uv = cube_uv_front(objectPoint);
      face = 2;
    }
    else if(objectNormal.z < 0){
      uv = cube_uv_back(objectPoint);
      face = 3;
    }
    else if(objectNormal.y > 0) {
      uv = cube_uv_up(objectPoint);
      face = 4;
    }
    else{
      uv = cube_uv_down(objectPoint);
      face = 5;
    }
    intersection->u = uv.x;
    intersection->v = uv.y;
    intersection->face = face;


    ray->tmax = tmin;
    return true;
}
