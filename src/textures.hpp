#pragma once

#include <memory>
#include <glm/glm.hpp>
#include "defines.h"
#include "image.h"
#include <iostream>

class texture {
    public:
        virtual color3 value(float u, float v) const = 0;
};

class solid_color : public texture {
    public:
        solid_color() {}
        solid_color(color3 c) : color_value(c) {}

        solid_color(float red, float green, float blue)
          : solid_color(color3(red,green,blue)) {}

        virtual color3 value(float u, float v) const override {
            return color_value;
        }

    private:
        color3 color_value;
};

class checker_texture : public texture {
    public:
        checker_texture() {}

        checker_texture(std::shared_ptr<texture> _even, std::shared_ptr<texture> _odd)
            : odd(_odd), even(_even) {}

        checker_texture(color3 c1, color3 c2)
            : odd(std::make_shared<solid_color>(c2)), even(std::make_shared<solid_color>(c1)) {}

        virtual color3 value(float u, float v) const override {
            int u2 = floor(u * width);
            int v2 = floor(v * height);
            if ( (u2 + v2) % 2 == 0)
              return odd->value(u, v);
            else
              return even->value(u, v);
        }

    public:
        std::shared_ptr<texture> odd;
        std::shared_ptr<texture> even;
        int width = 2;
        int height = 2;
};

class image_texture : public texture {
  public:
    image_texture() {}

    image_texture(const char* filename) {
      m_image = loadPng(filename);
    }

    virtual color3 value(float u, float v) const override {
      int u2 = floor((1.f-u) * (m_image->width-1));
      int v2 = floor((1.f-v) * (m_image->height-1));
      ////std::cout << u2 << std::endl;
      //std::cout << u2 << " " << v2 << std::endl;
      color3 color = *getPixelPtr(m_image, u2, v2);
      return color;
    }

  private:
    RenderImage* m_image;
};
