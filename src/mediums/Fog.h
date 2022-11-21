#include "../Medium.h"

class Fog : public Medium {
public:
  Fog(float density, bool isUniform, float pUniform=0.25, int phase=0, 
      bool uniformSampleRay=true) : m_beta(density), m_isUniform(isUniform),
      m_pUniform(pUniform), m_phase(phase), m_uniformSampleRay(uniformSampleRay) {}
  ~Fog() {}

  float tr(const Ray& ray, float yFloor) const override;
  color3 sample(const Ray &ray, Scene* scene, KdTree* tree, float yFloor) const override;
private:
  float int_exponential(float y0, float ysol, float beta, float s, float uy) const;

  float m_beta;
  bool m_isUniform;
  float m_pUniform;
  int m_phase;
  bool m_uniformSampleRay;
};
