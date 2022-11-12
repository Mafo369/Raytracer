#pragma once
#include "../Object.h"

class Cube : public Object {
  public:
    Cube(std::shared_ptr<Material> mat, Transform transform) : Object(mat, transform) {}
    bool intersect(Ray* ray, Intersection* intersection) const override;
    Intersection sample(const Intersection& inter, const point2& u) const override {return Intersection();}
    Intersection sample(const point2& u) const override {return Intersection();}
    float pdf(const Intersection& inter, const vec3& wi) const override {return 0.f;}
    
    void check_axis(float origin, float direction, float& tmin, float& tmax) const;
    vec3 computeCubeNormal(vec3 position) const;

    vec2 cube_uv_front(point3 point) const;
    vec2 cube_uv_back(point3 point) const;
    vec2 cube_uv_left(point3 point) const;
    vec2 cube_uv_right(point3 point) const;
    vec2 cube_uv_up(point3 point) const;
    vec2 cube_uv_down(point3 point) const;
};
