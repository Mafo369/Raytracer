#pragma once
#include "../image.h"
#include "../scene.h"

class ImageRenderer {
public:
  virtual void init(std::string& imageName, glm::vec2 resolution, Scene* scene) = 0;

  virtual void destroy() = 0;

  virtual void run() = 0;

protected:
  Scene* m_scene = nullptr;
  
  std::string m_imageName = "";
};
