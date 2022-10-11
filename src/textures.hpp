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
      virtual color3 value(float u, float v, int face = -1) const = 0;

      Transform m_transform;
};

class solid_color : public texture {
    public:
      solid_color(color3 c) : color_value(c) {}
      ~solid_color() {}

      solid_color(float red, float green, float blue)
        : solid_color(color3(red,green,blue)) {}

      color3 value(float u, float v, int face = -1) const override {
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

      color3 value(float u, float v, int face = -1) const override {
          point3 u1;
          u1.x = u - (int) u;
          u1.y = v - (int) v;
          if(u1.x < 0.0)
            u1.x += 1.0;
          if(u1.y < 0.0)
            u1.y += 1.0;
          if(u1.z < 0.0)
            u1.z += 1.0;
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
          //int u2 = floor(u * width);
          //int v2 = floor(v * height);
          //if ( (u2 + v2) % 2 == 0)
          //  return even->value(u, v);
          //else
          //  return odd->value(u, v);
      }

  public:
      std::shared_ptr<texture> odd;
      std::shared_ptr<texture> even;
      int width = 200;
      int height = 200;
};

class image_texture : public texture {
  public:
    image_texture(const char* filename) {
      m_image = loadPng(filename);
    }
    
    ~image_texture(){
      delete m_image;
    }

    color3 value(float u, float v, int face = -1) const override {
      int u2 = floor(u * (m_image->width-1));
      int v2 = floor(v * (m_image->height-1));
      return *getPixelPtr(m_image, u2, v2);
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

    color3 value(float u, float v, int face) const override {
      return m_faces[face]->value(u, v);
    }

  private:
    std::vector<texture*> m_faces;
};
