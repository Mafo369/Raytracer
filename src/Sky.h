#pragma once

#include "defines.h"
#include "ray.h"
#include <stb_image.h>
#include "sampling/sampling.h"
#include <iostream>

class Sky {
  public:
    virtual vec3 getRadiance(const Ray& ray) const = 0;
};

class UniformSky : public Sky {
  public:
    UniformSky(const vec3& color) : m_color(color) {};
    vec3 getRadiance(const Ray& ray) const {
      return m_color;
    }
    color3 m_color;
};

//class IBL : public Sky {
//  public:
//    IBL(const std::string& filename){
//      int n;
//      m_pixels = stbi_loadf(filename.c_str(), &m_width, &m_height, &n, 0);
//      std::cout << "Channels HDR: " << n << std::endl;
//      std::cout << "Resolution: " << m_width << "x" << m_height << std::endl;
//    }
//
//    ~IBL() {
//      stbi_image_free(m_pixels);
//    }
//
//    vec3 getRadiance(const Ray& ray) const override {
//          vec3 dir = m_transform.getInvTransform() * ray.dir;
//          double theta = std::acos(dir.y);
//          double phi = std::atan2(dir.z, dir.x);
//          if(phi<0)phi += 2*M_PI;
//
//          int i = phi/(2*M_PI) * m_width;
//          int j = theta/M_PI * m_height;
//
//          int index = 3*i + 3*m_width*j;
//
//          return vec3(m_pixels[index], m_pixels[index+1], m_pixels[index+2]);
//    };
//    
//    Transform m_transform;
//
//    int m_width;
//    int m_height;
//    float* m_pixels;
//};
