#define TINYOBJLOADER_IMPLEMENTATION

#include <iostream>
#include <string>

#include "Medium.h"
#include "Object.h"
#include "scene.h"
#include <algorithm>
#include <string.h>

//#include <glm/gtx/norm.hpp>

#include "Camera.h"
#include "Light.h"
#include "shapes/cube.h"
#include "shapes/plane.h"
#include "shapes/sphere.h"
#include "shapes/triangle.h"

#include "Material.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>
glm::vec3 extractScale( const glm::mat4& m ) {
    // length2 returns length squared i.e. v·v
    // no square root involved
    return glm::vec3( glm::length2( glm::vec3( m[0] ) ),
                      glm::length2( glm::vec3( m[1] ) ),
                      glm::length2( glm::vec3( m[2] ) ) );
}

Object* initSphere( std::shared_ptr<Material> mat, Transform transform ) {
    auto ret                = new Sphere( mat, transform );
    ret->geom.type          = SPHERE;
    ret->geom.sphere.center = transform.transformFrom( vec3( 0 ) );
    glm::vec3 scalesSq      = extractScale( transform.getTransform() );
    float const maxScaleSq  = *std::max_element(
        &scalesSq[0], &scalesSq[0] + scalesSq.length() ); // length gives the dimension here i.e. 3
    // one sqrt when you know the largest of the three
    float const largestScale = std::sqrt( maxScaleSq );
    ret->geom.sphere.radius  = 1.f * largestScale;
    return ret;
}

Object* initPlane( vec3 normal, float d, std::shared_ptr<Material> mat ) {
    // auto ret = new Plane(mat, Transform(1.f));
    // ret->geom.type = PLANE;
    // ret->geom.plane.normal = normalize(normal);
    // ret->geom.plane.dist = d;
    // return ret;
    return nullptr;
}

Object* initCube( std::shared_ptr<Material> mat, Transform transform ) {
    auto ret           = new Cube( mat, transform );
    ret->geom.type     = CUBE;
    ret->transform     = transform;
    ret->geom.cube.min = transform.transformFrom( vec3( -1 ) );
    ret->geom.cube.max = transform.transformFrom( vec3( 1 ) );
    return ret;
}

Object* initTriangle( vec3 p1,
                      vec3 p2,
                      vec3 p3,
                      vec3 n,
                      vec2 t[3],
                      Transform transform,
                      std::shared_ptr<Material> mat ) {
    auto ret                  = new Triangle( mat, transform );
    ret->geom.type            = TRIANGLE;
    ret->geom.triangle.p1     = p1;
    ret->geom.triangle.p2     = p2;
    ret->geom.triangle.p3     = p3;
    ret->geom.triangle.normal = n;
    ret->geom.triangle.tex[0] = t[0];
    ret->geom.triangle.tex[1] = t[1];
    ret->geom.triangle.tex[2] = t[2];
    return ret;
}

Object* initSmoothTriangle( vec3 p1,
                            vec3 p2,
                            vec3 p3,
                            vec3 n,
                            vec2 t[3],
                            vec3 n1,
                            vec3 n2,
                            vec3 n3,
                            Transform transform,
                            std::shared_ptr<Material> mat ) {
    auto ret                  = new Triangle( mat, transform );
    ret->geom.type            = TRIANGLE;
    ret->geom.triangle.p1     = p1;
    ret->geom.triangle.p2     = p2;
    ret->geom.triangle.p3     = p3;
    ret->geom.triangle.n1     = n1;
    ret->geom.triangle.n2     = n2;
    ret->geom.triangle.n3     = n3;
    ret->geom.triangle.normal = n;
    ret->geom.triangle.tex[0] = t[0];
    ret->geom.triangle.tex[1] = t[1];
    ret->geom.triangle.tex[2] = t[2];
    return ret;
}

void freeObject( Object* obj ) {
    free( obj );
}

Light* initPointLight( point3 position, color3 color, float size ) {
    Light* light = new PointLight( position, color, size );
    return light;
}

Light* initAmbientLight( color3 color ) {
    Light* light = new AmbientLight( vec3( 0, 0, 0 ), color );
    return light;
}

Light* initDirectLight( vec3 direction, color3 color ) {
    Light* light = new DirectLight( direction, color );
    return light;
}

void freeLight( Light* light ) {
    delete light;
}

Scene* initScene() {
    return new Scene;
}

void freeScene( Scene* scene ) {
    std::for_each( scene->objects.begin(), scene->objects.end(), freeObject );
    std::for_each( scene->lights.begin(), scene->lights.end(), freeLight );
    if ( scene->medium != nullptr ) delete scene->medium;
    delete scene;
}

void setCameraFOV( Scene* scene,
                   point3 position,
                   point3 at,
                   vec3 up,
                   float fov,
                   float width,
                   float height,
                   float aperture,
                   float dist_to_focus ) {
    scene->cam = new CameraFOV( position, at, up, fov, width, height, aperture, dist_to_focus );
}

void setSimpleCamera( Scene* scene,
                      point3 position,
                      point3 at,
                      vec3 up,
                      float fov,
                      float w,
                      float h,
                      float radius,
                      float focusDist ) {
    scene->cam = new SimpleCamera( position, at, up, fov, ( w / h ), w, h, radius, focusDist );
}

void addObject( Scene* scene, Object* obj ) {
    scene->objects.push_back( obj );
}

void addLight( Scene* scene, Light* light ) {
    scene->lights.push_back( light );
}

void setSkyColor( Scene* scene, color3 c ) {
    // scene->skyColor = c;
}

// https://github.com/tinyobjloader/tinyobjloader
void readObjToTriangleMesh( const char* file,
                            tinyobj::attrib_t& attrib,
                            std::vector<tinyobj::shape_t>& shapes,
                            std::vector<tinyobj::material_t>& materials ) {
    std::string inputfile = file;
    std::string err;

    bool ret =
        tinyobj::LoadObj( &attrib, &shapes, &materials, &err, inputfile.c_str(), "../assets/" );

    if ( !err.empty() ) { std::cout << err << std::endl; }

    if ( !ret ) exit( 1 );

    std::cout << "Model vertices: " << attrib.vertices.size() << std::endl;
}
