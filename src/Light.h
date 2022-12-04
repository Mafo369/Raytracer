#pragma once

#include "defines.h"
#include "scene.h"
#include "kdtree.h"
#include "Object.h"
#include <stb_image.h>
#include "sampling/sampling.h"

class Light {
  public:
    virtual ~Light() = default;
    Light() = default;
    vec3 getColor() { return m_color; } 
    vec3 getPosition() { return m_position; }
    virtual vec3 getDirection(point3 p) = 0;
    virtual float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) = 0;
    virtual std::vector<vec3> getSamples() { return m_samples; }
    virtual color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const = 0;
    virtual float pdf_Li(const Intersection& it, const vec3& wi) const = 0;
    virtual color3 Le(Ray* ray) const {
      return color3(0);
    }
    bool isAmbient() { return m_ambient; }
    void setAmbient(bool ambient) { m_ambient = ambient; }
    bool is_shadowed(vec3 lightPosition, vec3 normal, vec3 point, Scene* scene, KdTree* tree);
    virtual float getSize() { return m_size; }
    virtual bool isDirectional() { return false; }
    color3 m_color;
  protected:
    float m_size;
    bool m_ambient = false;
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
    color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const override { return color3(0); }
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
    color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const override { return color3(0); }
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
    color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const override { return color3(0); }
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
    color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const override { return color3(0); }
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
    color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const override;
    float pdf_Li(const Intersection& it, const vec3& wi) const override;
    ~ShapeLight() {}
  private:
    Object* m_shape;
    const float m_area = 1.f;
};

#include "Material.h"

class IBL : public Light {
  public:
    IBL(const std::string& filename){
      int n;
      m_pixels = stbi_loadf(filename.c_str(), &m_width, &m_height, &n, 0);
      float filter = (float)1 / std::max(m_width, m_height);
      std::unique_ptr<float[]> img(new float[m_width * m_height]);
      for (int v = 0; v < m_height-3; v+=3) {
          float vp = (float)v / (float)m_height;
          float sinTheta = std::sin(Pi * float(v + .5f) / float(m_height));
          for (int u = 0; u < m_width-3; u+=3) {
              float up = (float)u / (float)m_width;
              img[u + v * m_width] = luminance(color3(m_pixels[u + v * m_width], m_pixels[u + v * m_width + 1], m_pixels[u + v * m_width + 2]));
              img[u + v * m_width] *= sinTheta;
          }
      }
      std::cout << "Channels HDR: " << n << std::endl;
      std::cout << "Resolution: " << m_width << "x" << m_height << std::endl;
      m_distribution.reset(new Distribution2D(m_pixels, m_width, m_height));
    }

    float intensityAt(vec3 point, Scene* scene, KdTree* tree, vec3 view, Intersection* intersection) override { return 0; }
    vec3 getDirection(point3 p) override { return vec3(0); };
    vec3 getLightPoint(point3 p, int c, float r) {return vec3(0);}

    ~IBL() {
      stbi_image_free(m_pixels);
    }

    vec3 Le(Ray* ray) const override {
          vec3 dir = m_transform.getInvTransform() * ray->dir;
          double theta = std::acos(dir.y);
          double phi = std::atan2(dir.z, dir.x);
          if(phi<0)phi += 2*M_PI;

          int i = phi/(2*M_PI) * m_width;
          int j = theta/M_PI * m_height;

          int index = 3*i + 3*m_width*j;

          return vec3(m_pixels[index], m_pixels[index+1], m_pixels[index+2]);
    };

    color3 sample_Li(Scene* scene, KdTree* tree, const Intersection& inter, const point2& u, vec3* wi, float* pdf, bool* visibility) const override;
    float pdf_Li(const Intersection& it, const vec3& wi) const override;
    
    Transform m_transform;

    std::unique_ptr<Distribution2D> m_distribution;

    int m_width;
    int m_height;
    float* m_pixels;
};
