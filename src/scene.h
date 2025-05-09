#pragma once

#include "Material.h"
#include "defines.h"
#include "tiny_obj_loader.h"

// #include "Object.h"
#include "textures.hpp"

class Camera;
class Light;
class Object;
class Medium;
class Sky;

typedef std::vector<Object*> Objects;
typedef std::vector<Light*> Lights;

class Scene
{
  public:
    Lights lights;    //! the scene have several lights
    Lights envLights; //! the scene have several lights
    Objects objects;  //! the scene have several objects
    std::vector<Material> m_Materials;
    Camera* cam; //! the scene have one camera
    // color3 skyColor; //! the sky color, could be extended to a sky function ;)
    // texture* m_skyTexture = nullptr;
    Sky* sky;
    int depth      = 8;
    Medium* medium = nullptr;
    float ysol     = -12.f;

    Material& CreateMaterial( MaterialModel materialModel, MatType materialType = DIFFUSE ) {
        size_t id = m_Materials.size();
        return m_Materials.emplace_back( id, materialModel, materialType );
    }

    Material& GetMaterial( int index ) { return m_Materials[index]; }
};

//! create a new sphere structure
Object* initSphere( int materialIndex, Transform transform );
Object* initCube( int materialIndex, Transform transform );
Object* initPlane( vec3 normal, float d, int materialIndex );
Object* initTriangle( vec3 p1,
                      vec3 p2,
                      vec3 p3,
                      vec3 n,
                      vec2 t[3],
                      Transform transform,
                      int materialIndex );
Object* initSmoothTriangle( vec3 p1,
                            vec3 p2,
                            vec3 p3,
                            vec3 n,
                            vec2 t[3],
                            vec3 n1,
                            vec3 n2,
                            vec3 n3,
                            Transform transform,
                            int materialIndex );

//! release memory for the object obj
void freeObject( Object* obj );

//! init a new light at position with a give color (no special unit here for the moment)
Light* initPointLight( point3 position, color3 color, float size = 0 );
Light* initDirectLight( vec3 direction, color3 color );

Light* initAmbientLight( color3 color );

//! release memory for the light
void freeLight( Light* );

// allocate the momery for the scene
Scene* initScene();
void freeScene( Scene* scene );

void setCameraFOV( Scene* scene,
                   point3 position,
                   vec3 at,
                   vec3 up,
                   float fov,
                   float width,
                   float height,
                   float aperture      = 0.01,
                   float dist_to_focus = 1 );
void setSimpleCamera( Scene* scene,
                      point3 position,
                      vec3 at,
                      vec3 up,
                      float fov,
                      float w,
                      float h,
                      float radius,
                      float distFocus );

//! take ownership of obj freeScene will free obj) ... typically use addObject(scene, initPlane()
void addObject( Scene* scene, Object* obj );

//! take ownership of light : freeScene will free light) ... typically use addObject(scene,
//! initLight()
void addLight( Scene* scene, Light* light );

void setSkyColor( Scene* scene, color3 c );

void readObjToTriangleMesh( const char* file,
                            tinyobj::attrib_t& attrib,
                            std::vector<tinyobj::shape_t>& shapes,
                            std::vector<tinyobj::material_t>& materials );
