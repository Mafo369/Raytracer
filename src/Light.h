#pragma once

#include "defines.h"
#include "scene.h"
#include "kdtree.h"
#include "Object.h"

class Light {
  public:
    virtual ~Light() = default;
    Light() = default;
    vec3 getColor() { return m_color; } 
    vec3 getPosition() { return m_position; }
    virtual vec3 getDirection(point3 p) = 0;
    virtual float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) = 0;
    virtual std::vector<vec3> getSamples() { return m_samples; }
    virtual color3 sample_Li(const Intersection& inter, const point2& u, vec3* wi, float* pdf) const = 0;
    virtual float pdf_Li(const Intersection& it, const vec3& wi) const = 0;
    bool isAmbient() { return m_ambient; }
    void setAmbient(bool ambient) { m_ambient = ambient; }
    bool is_shadowed(vec3 lightPosition, vec3 normal, vec3 point, Scene* scene, KdTree* tree);
    virtual float getSize() { return m_size; }
    virtual bool isDirectional() { return false; }
  protected:
    float m_size;
    bool m_ambient = false;
    color3 m_color;
    vec3 m_position;
    std::vector<vec3> m_samples;
};

#include <random>

class PointLight : public Light {
  public:
    PointLight(vec3 position , color3 color, float size=0.0 );
    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override;
    vec3 getDirection(point3 p) override;
    vec3 getLightPoint(point3 p, int c, float r);
    color3 sample_Li(const Intersection& inter, const point2& u, vec3* wi, float* pdf) const override { return color3(0); }
    float pdf_Li(const Intersection& it, const vec3& wi) const override { return 0.f; }
    ~PointLight();
  private:
    int m_shadowMin;
    int m_shadowMax;
};

class AmbientLight : public Light {
  public:
    AmbientLight(vec3 position , color3 color);
    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override;
    vec3 getDirection(point3 p) override;
    ~AmbientLight();
    color3 sample_Li(const Intersection& inter, const point2& u, vec3* wi, float* pdf) const override { return color3(0); }
    float pdf_Li(const Intersection& it, const vec3& wi) const override { return 0.f; }
  private:
};

class DirectLight : public Light {
  public:
    DirectLight(vec3 direction , color3 color);
    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override;
    vec3 getDirection(point3 p) override;
    bool isDirectional() override { return true; }
    ~DirectLight();
    color3 sample_Li(const Intersection& inter, const point2& u, vec3* wi, float* pdf) const override { return color3(0); }
    float pdf_Li(const Intersection& it, const vec3& wi) const override { return 0.f; }
  private:
};

class AreaLight : public Light {
  public:
    AreaLight(vec3 corner, vec3 full_uvec, int usteps, vec3 full_vvec, int vsteps, vec3 intensity);
    ~AreaLight();
    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override;
    point3 pointOnLight(float u, float v);
    vec3 getDirection(point3 p) override;
    void setup(Scene* scene);
    color3 sample_Li(const Intersection& inter, const point2& u, vec3* wi, float* pdf) const override { return color3(0); }
    float pdf_Li(const Intersection& it, const vec3& wi) const override { return 0.f; }
  private:
    float intensity;
    vec3 m_corner;
    vec3 uvec;
    vec3 vvec;
    int m_usteps;
    int m_vsteps;
    int nbSamples;


    Object* m_t1;
    Object* m_t2;
};

class ShapeLight : public Light {
  public:
    ShapeLight(vec3 position , color3 color, Object* obj) : Light() { m_shape = obj; m_position = position; m_color = color; }
    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override { return 0; }
    vec3 getDirection(point3 p) override { return vec3(0); };
    vec3 getLightPoint(point3 p, int c, float r) {return vec3(0);}
    color3 L(const Intersection& inter, const vec3& w) const {
      return dot(inter.normal, w) > 0.f ? m_color : color3(0);
    }
    color3 sample_Li(const Intersection& inter, const point2& u, vec3* wi, float* pdf) const override;
    float pdf_Li(const Intersection& it, const vec3& wi) const override;
    ~ShapeLight() {}
  private:
    Object* m_shape;
    const float m_area = 1.f;
};
