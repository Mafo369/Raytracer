#pragma once

#define TP11
#define TP12
#define TP13
#define TP14
#define TP15
#define TP21
#define TP22
#define TP23
#define TP23TR
#define TP24
#define TP31
#define TP32
#define IMP
#undef  TP23TR
//#define SAMPLEGLOSSY

#include <stdbool.h>
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_access.hpp>
#include <glm/gtc/matrix_transform.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/matrix_transform_2d.hpp>
using namespace glm;

typedef vec3 point3;
typedef vec2 point2;
typedef vec3 color3;

const float acne_eps = 1e-3;

namespace glm {
// not really needed, but it makes it easier to follow the book...
template <int N, typename T, qualifier P> T length_sq(const vec<N, T, P> &x) { return dot(x, x); }
} // namespace glm

// Utility Functions
inline float degrees_to_radians(float degrees) {
  return degrees * M_PI / 180.0f;
}


// Inspired by: https://graphics.cs.utah.edu/courses/cs6620/fall2019/?f=code&prj=7&file=scene.h
class Transform {
  public:
    Transform() : m_translation(0,0,0) { m_transform = mat4(1.f);
                                         m_invTransform = inverse(m_transform); }
  
    mat3 const& getTransform() const { return m_transform; }
    vec3 const& getTranslation() const { return m_translation; }
    mat3 const& getInvTransform() const { return m_invTransform; }

    vec3 transformTo( vec3 const& p ) const { return m_invTransform * (p - m_translation); }
    vec3 transformFrom( vec3 const& p ) const { return m_transform * p + m_translation; }

    vec3 vectorTransformTo( vec3 const& dir ) const { return transposeMult(m_transform, dir); }
    vec3 vectorTransformFrom( vec3 const& dir ) const { return transposeMult(m_invTransform, dir); }

    void translate(vec3 const& p) { m_translation += p; }
    void rotate(vec3 const& axis, float degrees) { 
      mat3 m = glm::mat3(0.f);
      auto theta = degrees_to_radians(degrees);
      float c = (float) cos(theta); 
      if(c == 1){
        m = glm::mat4(1.f);
        return;
      }
      float s = (float) sin(theta);
      float t = 1 - c;
      float tx = t * axis.x;
      float ty = t * axis.y;
      float tz = t * axis.z;
      float txy = tx * axis.y;
      float txz = tx * axis.z;
      float tyz = ty * axis.z;
      float sx = s * axis.x;
      float sy = s * axis.y;
      float sz = s * axis.z;
      
      m[0][0] = tx * axis.x + c;
      m[0][1] = txy + sz;
      m[0][2] = txz - sy;

      m[1][0] = txy - sz;
      m[1][1] = ty * axis.y + c;
      m[1][2] = tyz + sx;

      m[2][0] = txz + sy;
      m[2][1] = tyz - sx;
      m[2][2] = tz * axis.z + c;
      transform(m);
    }
    void scale(float sx, float sy, float sz) {
      auto m = mat3(0.f);
      m[0][0] = sx; m[1][1] = sy; m[2][2] = sz;
      transform(m);
    }

    void transform(mat3 const& m) { 
      m_transform = m * m_transform;
      m_translation = m * m_translation;
      m_invTransform = glm::inverse(m_transform);
    }

  private:

    // Multiplies the given vector with the transpose of the given matrix
    static vec3 transposeMult( mat3 const &m, vec3 const &dir )
    {
        vec3 d;
        d.x = dot(column(m, 0), dir);
        d.y = dot(column(m, 1), dir);
        d.z = dot(column(m, 2), dir);
        return d;
    }

    glm::mat3 m_transform;
    vec3 m_translation;
    mutable glm::mat3 m_invTransform;
};
