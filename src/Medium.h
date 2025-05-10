
#include "Object.h"
#include "defines.h"
#include "ray.h"
#include "scene.h"

class Medium
{
  public:
    virtual ~Medium() {};
    virtual float tr( const Ray& ray, float yFloor ) const                                  = 0;
    virtual color3 sample( const Ray& ray, Scene* scene, KdTree* tree, float yFloor ) const = 0;
};
