#pragma once

#include <memory>
#include <glm/glm.hpp>
#include "defines.h"
#include "image.h"
#include <iostream>
#include <algorithm>

class texture {
    public:
      texture() = default;
      virtual ~texture() = default;
      virtual color3 value(float u, float v) const = 0;
      virtual color3 value(float u, float v, int face) const = 0;

      virtual color3 value(float u, float v, const vec3 duv[2]) const {
          vec3 uvw = vec3(u, v, 0);
          vec3 uv = m_transform.transformTo(uvw);
          color3 texColor = value(uv.x, uv.y);
          if(dot(duv[0], duv[0]) + dot(duv[1], duv[1]) == 0) return texColor;

          //std::cout << "DERIVATIVE" << std::endl;
          //std::cout <<"uv color: " << glm::to_string(texColor) << std::endl;
          //std::cout << glm::to_string(uv) << std::endl;
          for(int i = 1; i < 32; i++){
            float x=0, y=0, fx=0.05f, fy=1.0f/3.0f;
            for ( int ix=i; ix>0; ix/=2 ) { x+=fx*(ix%2); fx/=2; }   // Halton sequence (base 2)
            for ( int iy=i; iy>0; iy/=3 ) { y+=fy*(iy%3); fy/=3; }   // Halton sequence (base 3)
            if(true){
              float r = sqrtf(x)*0.5f;
              x = r*sinf(y*(float)M_PI*2.0);
              y = r*cosf(y*(float)M_PI*2.0);
            }else
            {
              if ( x > 0.5f ) x-=1;
              if ( y > 0.5f ) y-=1;
            }
            vec3 new_uv = uv + x * duv[0] + y * duv[1]; 
            auto tdColor = value(new_uv.x, new_uv.y);
            //std::cout << glm::to_string(new_uv) << std::endl;
            //std::cout << "x: " << x << " y: " << y << std::endl;
            //std::cout << glm::to_string(duv[0]) << std::endl;
            //std::cout << glm::to_string(duv[1]) << std::endl;
            //std::cout << glm::to_string(tdColor) << std::endl;
            texColor += tdColor;
          }
          texColor = texColor / 32.f;
          return texColor;
      }

      Transform m_transform;
};

class solid_color : public texture {
    public:
      solid_color(color3 c) : color_value(c) {}
      ~solid_color() {}

      solid_color(float red, float green, float blue)
        : solid_color(color3(red,green,blue)) {}

      color3 value(float u, float v) const override {
          return color_value;
      }

      color3 value(float u, float v, int face) const override {
          return color_value;
      }

    private:
        color3 color_value;
};

class checker_texture : public texture {
  public:

      checker_texture(std::shared_ptr<texture> _even, std::shared_ptr<texture> _odd)
          : odd(_odd), even(_even) {}

      checker_texture(color3 c1, color3 c2)
          : odd(std::make_shared<solid_color>(c2)), even(std::make_shared<solid_color>(c1)) {}

      ~checker_texture() {}

      color3 value(float u, float v) const override {
          point3 u1;
          u1.x = u - (int) u;
          u1.y = v - (int) v;
          if(u1.x < 0.0)
            u1.x += 1.0;
          if(u1.y < 0.0)
            u1.y += 1.0;

          if(u1.x <= 0.5){
            if(u1.y <= 0.5)
              return even->value(u1.x, u1.y);
            else
              return odd->value(u1.x, u1.y);
          }else{
            if(u1.y <= 0.5)
              return odd->value(u1.x, u1.y);
            else
              return even->value(u1.x, u1.y);
          }
      }

      color3 value(float u, float v, int face) const override {
        return value(u, v);
      }

  public:
      std::shared_ptr<texture> odd;
      std::shared_ptr<texture> even;
};

class image_texture : public texture {
  public:
    image_texture(const char* filename) {
      m_image = loadPng(filename);
    }
    
    ~image_texture(){
      delete m_image;
    }

    color3 value(float u, float v) const override {
      int u2 = floor(std::clamp(u,0.f,1.f) * (m_image->width-1));
      int v2 = floor(std::clamp(v,0.f,1.f) * (m_image->height-1));
      return *getPixelPtr(m_image, u2, v2);
    }

    color3 value(float u, float v, int face) const override {
      return value(u, v);
    }

  private:
    RenderImage* m_image;
};

typedef struct s_FaceInfo{
  color3 ul;
  color3 ur;
  color3 bl;
  color3 br;
  color3 main;
}FaceInfo;

class AlignCheck : public texture {
  public:
    AlignCheck(FaceInfo left, FaceInfo right, FaceInfo front, FaceInfo back, FaceInfo up, FaceInfo down){
      m_faces.push_back(left);
      m_faces.push_back(right);
      m_faces.push_back(front);
      m_faces.push_back(back);
      m_faces.push_back(up);
      m_faces.push_back(down);
    }

    ~AlignCheck() {}

    color3 value(float u, float v) const override {
      std::cerr << 
        "Warning: This texture has multiple faces. Selecting default one..." << std::endl;
      return m_faces[0].main;
    }

    color3 value(float u, float v, int face) const override {
      if( v > 0.8 ){
        if( u < 0.2 ) return m_faces[face].ul;
        if( u > 0.8 ) return m_faces[face].ur;
      }
      else if ( v < 0.2 )
      {
        if ( u < 0.2 ) return m_faces[face].bl;
        if ( u > 0.8 ) return m_faces[face].br;
      }
      return m_faces[face].main;
    }

  private:
    std::vector<FaceInfo> m_faces;
};

class CubeMapTexture : public texture {
  public:
    CubeMapTexture(texture* left, texture* right, texture* front, texture* back, texture* up, texture* down){
      m_faces.push_back(left);
      m_faces.push_back(right);
      m_faces.push_back(front);
      m_faces.push_back(back);
      m_faces.push_back(up);
      m_faces.push_back(down);
    }

    ~CubeMapTexture() {}

    color3 value(float u, float v) const override {
      std::cerr << 
        "Warning: This texture has multiple faces. Selecting default one..." << std::endl;
      return m_faces[0]->value(u, v);
    }

    color3 value(float u, float v, int face) const override {
      return m_faces[face]->value(u, v);
    }

  private:
    std::vector<texture*> m_faces;
};
