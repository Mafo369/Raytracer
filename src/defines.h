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
#include <glm/gtx/compatibility.hpp>
using namespace glm;

#define COMBINE_BRDFS_WITH_FRESNEL 1

typedef vec3 point3;
typedef vec2 point2;
typedef vec3 color3;

const float acne_eps = 1e-4;

namespace glm {
// not really needed, but it makes it easier to follow the book...
template <int N, typename T, qualifier P> T length_sq(const vec<N, T, P> &x) { return dot(x, x); }
} // namespace glm

inline bool isBlack(const vec3& color){
  return color.r + color.g + color.b == 0.;
}

// Utility Functions
inline float degrees_to_radians(float degrees) {
  return degrees * M_PI / 180.0f;
}

inline float sqr(float x) { return x*x; }

#include <random>

static thread_local std::default_random_engine engine;
static thread_local std::uniform_real_distribution<float> uniform01 {0, 1};

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

inline void
CoordinateSystem(const vec3 &v1, vec3 *v2, vec3 *v3) {
    if (std::abs(v1.x) > std::abs(v1.y))
        *v2 = vec3(-v1.z, 0, v1.x) /
              std::sqrt(v1.x * v1.x + v1.z * v1.z);
    else
        *v2 = vec3(0, v1.z, -v1.y) /
              std::sqrt(v1.y * v1.y + v1.z * v1.z);
    *v3 = cross(v1, *v2);
}

inline vec3 SphericalDirection(float sinTheta, float cosTheta, float phi, 
                               const vec3& x, const vec3& y, const vec3& z) {
  return sinTheta * cos(phi) * x + sinTheta * sin(phi) * y + cosTheta * z;
}

#include <vulkan/vulkan.h>
#include <unordered_map>
#include <vector>
#include <iostream>

static std::unordered_map<VkResult, std::string> ErrorDescriptions = {
    {VK_SUCCESS, "Command successfully completed"},
    {VK_NOT_READY, "A fence or query has not yet completed"},
    {VK_TIMEOUT, "A wait operation has not completed in the specified time"},
    {VK_EVENT_SET, "An event is signaled"},
    {VK_EVENT_RESET, "An event is unsignaled"},
    {VK_INCOMPLETE, "A return array was too small for the result"},
    {VK_SUBOPTIMAL_KHR, "A swapchain no longer matches the surface properties exactly, but can still be used to present to the surface successfully."},
    {VK_THREAD_IDLE_KHR, "A deferred operation is not complete but there is currently no work for this thread to do at the time of this call."},
    {VK_THREAD_DONE_KHR, "A deferred operation is not complete but there is no work remaining to assign to additional threads."},
    {VK_OPERATION_DEFERRED_KHR, "A deferred operation was requested and at least some of the work was deferred."},
    {VK_OPERATION_NOT_DEFERRED_KHR, "A deferred operation was requested and no operations were deferred."},
    {VK_PIPELINE_COMPILE_REQUIRED_EXT, "A requested pipeline creation would have required compilation, but the application requested compilation to not be performed."},
    {VK_ERROR_OUT_OF_HOST_MEMORY, "A host memory allocation has failed."},
    {VK_ERROR_OUT_OF_DEVICE_MEMORY, "A device memory allocation has failed."},
    {VK_ERROR_INITIALIZATION_FAILED, "Initialization of an object could not be completed for implementation-specific reasons."},
    {VK_ERROR_DEVICE_LOST, "The logical or physical device has been lost. See Lost Device"},
    {VK_ERROR_MEMORY_MAP_FAILED, "Mapping of a memory object has failed."},
    {VK_ERROR_LAYER_NOT_PRESENT, "A requested layer is not present or could not be loaded."},
    {VK_ERROR_EXTENSION_NOT_PRESENT, "A requested extension is not supported."},
    {VK_ERROR_FEATURE_NOT_PRESENT, "A requested feature is not supported."},
    {VK_ERROR_INCOMPATIBLE_DRIVER, "The requested version of Vulkan is not supported by the driver or is otherwise incompatible for implementation-specific reasons."},
    {VK_ERROR_TOO_MANY_OBJECTS, "Too many objects of the type have already been created."},
    {VK_ERROR_FORMAT_NOT_SUPPORTED, "A requested format is not supported on this device."},
    {VK_ERROR_FRAGMENTED_POOL, "A pool allocation has failed due to fragmentation of the pool’s memory. This must only be returned if no attempt to allocate host or device memory was made to accommodate the new allocation. This should be returned in preference to VK_ERROR_OUT_OF_POOL_MEMORY, but only if the implementation is certain that the pool allocation failure was due to fragmentation."},
    {VK_ERROR_SURFACE_LOST_KHR, "A surface is no longer available."},
    {VK_ERROR_NATIVE_WINDOW_IN_USE_KHR, "The requested window is already in use by Vulkan or another API in a manner which prevents it from being used again."},
    {VK_ERROR_OUT_OF_DATE_KHR, "A surface has changed in such a way that it is no longer compatible with the swapchain, and further presentation requests using the swapchain will fail. Applications must query the new surface properties and recreate their swapchain if they wish to continue presenting to the surface."},
    {VK_ERROR_INCOMPATIBLE_DISPLAY_KHR, "The display used by a swapchain does not use the same presentable image layout, or is incompatible in a way that prevents sharing an image."},
    {VK_ERROR_INVALID_SHADER_NV, "One or more shaders failed to compile or link. More details are reported back to the application via VK_EXT_debug_report if enabled."},
    {VK_ERROR_OUT_OF_POOL_MEMORY, "A pool memory allocation has failed. This must only be returned if no attempt to allocate host or device memory was made to accommodate the new allocation. If the failure was definitely due to fragmentation of the pool, VK_ERROR_FRAGMENTED_POOL should be returned instead."},
    {VK_ERROR_INVALID_EXTERNAL_HANDLE, "An external handle is not a valid handle of the specified type."},
    {VK_ERROR_FRAGMENTATION, "A descriptor pool creation has failed due to fragmentation."},
    {VK_ERROR_INVALID_DEVICE_ADDRESS_EXT, "A buffer creation failed because the requested address is not available."},
    {VK_ERROR_INVALID_OPAQUE_CAPTURE_ADDRESS, "A buffer creation or memory allocation failed because the requested address is not available. A shader group handle assignment failed because the requested shader group handle information is no longer valid."},
    {VK_ERROR_FULL_SCREEN_EXCLUSIVE_MODE_LOST_EXT, "An operation on a swapchain created with VK_FULL_SCREEN_EXCLUSIVE_APPLICATION_CONTROLLED_EXT failed as it did not have exlusive full-screen access. This may occur due to implementation-dependent reasons, outside of the application’s control."},
    {VK_ERROR_UNKNOWN, "An unknown error has occurred; either the application has provided invalid input, or an implementation failure has occurred."}
};


#define VK_CHECK(x)  VK_CHECK_CALL(x)

#define VK_CHECK_CALL(x) VulkanCheckErrorStatus(x, __FILE__, __LINE__)

static bool VulkanCheckErrorStatus(VkResult x, const char* file, int line)
{
    if(x != VK_SUCCESS)
    {
        std::cout << "\033[1;31;49m **Vulkan Function Call Error** Description : \033[0m" << ErrorDescriptions[x] << " \033[2;90;49m [at Line : " << line << " in File : " << file << "\033[0m]" << std::endl;
        return true;
    }
    else return false;
}

#define VK_LOG_SUCCESS(msg) std::cout << "\033[1;32m[VULKAN]\033[1;32m - SUCCESS : " << (msg) <<  " \033[0m\n";

#define VK_LOG(...) std::cout , "\033[1;32m[VULKAN]\033[1;33m - LOG : " , __VA_ARGS__ , " \033[0m" , std::endl

template <typename T>
static std::ostream& operator,(std::ostream& out, const T& t) {
  out << t;
  return out;
}

//overloaded version to handle all those special std::endl and others...
static std::ostream& operator,(std::ostream& out, std::ostream&(*f)(std::ostream&)) {
  out << f;
  return out;
}
