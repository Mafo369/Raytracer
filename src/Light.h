#pragma once

#include "defines.h"
#include "scene.h"
#include "kdtree.h"

class Light {
  public:
    virtual ~Light() = default;
    Light() = default;
    vec3 getColor() { return m_color; } 
    vec3 getPosition() { return m_position; }
    virtual vec3 getDirection(point3 p) = 0;
    virtual float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) = 0;
    virtual std::vector<vec3> getSamples() { return m_samples; }
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
  private:
};

class DirectLight : public Light {
  public:
    DirectLight(vec3 direction , color3 color);
    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override;
    vec3 getDirection(point3 p) override;
    bool isDirectional() override { return true; }
    ~DirectLight();
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

//
//#include "raytracer.h"
//
//class AreaLight : public Light {
//  public:
//    // AreaLight Interface
//    AreaLight() = default;
//    virtual color3 L(const Intersection &intr, const vec3 &w) const = 0;
//};
//
//class DiffuseAreaLight : public AreaLight {
//  public:
//    // DiffuseAreaLight Public Methods
//    DiffuseAreaLight(const color3 &Lemit, int nSamples,
//                     Object *shape,
//                     bool twoSided);
//    color3 L(const Intersection &intr, const vec3 &w) const {
//        return (twoSided || dot(intr.normal, w) > 0) ? Lemit : color3(0.f);
//    }
//    color3 Power() const;
//    color3 Sample_Li(const Intersection &ref, const point2 &u, vec3 *wo,
//                       float *pdf) const;
//    float Pdf_Li(const Intersection &, const vec3 &, Scene*, KdTree*) const;
//    color3 Sample_Le(const point2 &u1, const point2 &u2, float time,
//                       Ray *ray, vec3 *nLight, float *pdfPos,
//                       float *pdfDir) const;
//    void Pdf_Le(const Ray &, const vec3 &, float *pdfPos,
//                float *pdfDir) const;
//
//  protected:
//    // DiffuseAreaLight Protected Data
//    const color3 Lemit;
//    Object* shape;
//    // Added after book publication: by default, DiffuseAreaLights still
//    // only emit in the hemimsphere around the surface normal.  However,
//    // this behavior can now be overridden to give emission on both sides.
//    const bool twoSided;
//    const float area;
//};
//
