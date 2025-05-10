#include "example_scenes.h"
#include "Light.h"
#include "Material.h"
#include "Object.h"
#include "Sky.h"
#include "defines.h"
#include "image.h"
#include "materials/Blinn.h"
#include "materials/CookTorrance.h"
#include "ray.h"
#include "scene.h"
#include "shapes/plane.h"
#include "textures.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mediums/Fog.h"

#include <glm/gtc/matrix_transform.hpp>

void addObjectsFromFile( const char* filename,
                         Scene* scene,
                         Material& default_mat,
                         Transform& transform ) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    readObjToTriangleMesh( filename, attrib, shapes, materials );

    for ( size_t s = 0; s < shapes.size(); s++ ) {
        size_t index_offset = 0;
        for ( size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++ ) {
            int fv = shapes[s].mesh.num_face_vertices[f];

            std::vector<point3> vector;
            std::vector<point3> normals;
            std::vector<vec2> texture;
            std::vector<vec2> textures( fv );

            // Loop over vertices in the face.
            for ( int v = 0; v < fv; v++ ) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx   = attrib.vertices[3 * idx.vertex_index + 0];
                tinyobj::real_t vy   = attrib.vertices[3 * idx.vertex_index + 1];
                tinyobj::real_t vz   = attrib.vertices[3 * idx.vertex_index + 2];
                if ( !attrib.normals.empty() ) {
                    tinyobj::real_t nx = attrib.normals[3 * idx.normal_index + 0];
                    tinyobj::real_t ny = attrib.normals[3 * idx.normal_index + 1];
                    tinyobj::real_t nz = attrib.normals[3 * idx.normal_index + 2];
                    vec3 n             = point3( nx, ny, nz );
                    normals.push_back( n );
                }
                if ( !attrib.texcoords.empty() ) {
                    tinyobj::real_t tx = attrib.texcoords[2 * idx.texcoord_index + 0];
                    tinyobj::real_t ty = attrib.texcoords[2 * idx.texcoord_index + 1];
                    textures[v]        = vec2( tx, ty );
                }
                // Optional: vertex colors
                // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
                // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
                // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
                vec3 v0 = point3( vx, vy, vz );
                vector.push_back( v0 );
            }

            index_offset += fv;

            vec3 v0 = vector[0];
            vec3 v1 = vector[1];
            vec3 v2 = vector[2];
            vec3 n;
            if ( normals.empty() ) {
                vec3 v2v0 = v2 - v0;
                vec3 v1v0 = v1 - v0;
                n         = normalize( cross( v1v0, v2v0 ) );
                normals   = { n, n, n };
            }
            else { n = normalize( normals[0] + normals[1] + normals[2] ); }

            if ( !materials.empty() ) {
                int matIndex    = shapes[s].mesh.material_ids[f];
                auto& mat       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
                mat.m_IOR       = materials[matIndex].ior;
                mat.m_roughness = 1.f;
                mat.m_metalness = 0.0;
                mat.m_albedo.r  = materials[matIndex].diffuse[0];
                mat.m_albedo.g  = materials[matIndex].diffuse[1];
                mat.m_albedo.b  = materials[matIndex].diffuse[2];

                if ( materials[matIndex].emission[0] != 0 || materials[matIndex].emission[1] != 0 ||
                     materials[matIndex].emission[2] != 0 ) {
                    mat.m_emission = color3( materials[matIndex].emission[0],
                                             materials[matIndex].emission[1],
                                             materials[matIndex].emission[2] );
                    auto lightObj  = initSmoothTriangle( v0,
                                                        v1,
                                                        v2,
                                                        n,
                                                        textures.data(),
                                                        normals[0],
                                                        normals[1],
                                                        normals[2],
                                                        transform,
                                                        mat.m_UID );
                    auto light = new ShapeLight( vec3( 0 ), color3( mat.m_emission ), lightObj );
                    addLight( scene, light );
                }
                // mat.m_texture = new image_texture(materials[matIndex].diffuse_texname.c_str());
                if ( materials[matIndex].specular[0] == 1 && materials[matIndex].specular[1] == 1 &&
                     materials[matIndex].specular[2] == 1 ) {
                    mat.m_MatType   = SPECULAR;
                    mat.m_metalness = 1.0;
                    mat.m_roughness = 0.002;
                    // mat->m_specularColor.r = materials[matIndex].diffuse[0];
                    // mat->m_specularColor.g = materials[matIndex].diffuse[1];
                    // mat->m_specularColor.b = materials[matIndex].diffuse[2];
                }
                else {
                    // mat->m_specularColor.r = materials[matIndex].specular[0];
                    // mat->m_specularColor.g = materials[matIndex].specular[1];
                    // mat->m_specularColor.b = materials[matIndex].specular[2];
                }
                addObject( scene,
                           initSmoothTriangle( v0,
                                               v1,
                                               v2,
                                               n,
                                               textures.data(),
                                               normals[0],
                                               normals[1],
                                               normals[2],
                                               transform,
                                               mat.m_UID ) );
            }
            else {
                addObject( scene,
                           initSmoothTriangle( v0,
                                               v1,
                                               v2,
                                               n,
                                               textures.data(),
                                               normals[0],
                                               normals[1],
                                               normals[2],
                                               transform,
                                               default_mat.m_UID ) );
            }
        }
    }
}

Scene* initScene9() {
    Scene* scene = initScene();
    auto from    = point3( -0.23, 2.585, 5.3 );
    auto at      = vec3( -0.23, 2.585, -2.8 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 0.1, 0 ),
                     60,
                     float( WIDTH ),
                     float( HEIGHT ),
                     0.001,
                     glm::distance( from, at ) );
    // setSkyColor( scene, color3( 1.f ) );
    // auto sky   = new IBL( "/home/mafo/dev/Raytracer/assets/brown_photostudio_06_4k.hdr" );

    // addLight(scene, sky);
    // scene->envLights.push_back(sky);

    auto& mat       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat.m_roughness = 0.004;
    mat.m_metalness = 0.0000;
    mat.m_albedo    = color3( 0.f, 1.f, 0.95f );
    Transform t;

    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/cornell-box.obj", scene, mat, t );

    Transform tMonkey;
    tMonkey.scale( 0.58, 0.58, 0.58 );
    tMonkey.rotate( vec3( -1, 0, 0 ), 35 );
    tMonkey.rotate( vec3( 0, 1, 0 ), 20 );
    tMonkey.translate( vec3( -1.1, 3.43, -3.6 ) );
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/Suzanne.obj", scene, mat, tMonkey );

    Transform tPot;
    tPot.scale( 0.05, 0.05, 0.05 );
    tPot.rotate( vec3( -1, 0, 0 ), 90 );
    tPot.translate( vec3( -1.1, -0.15, -1 ) );
    // addObjectsFromFile( "../assets/teapot.obj", scene, mat, tPot );

    auto& mat2       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, TRANSPARENT );
    mat2.m_IOR       = 1.5;
    mat2.m_roughness = 0.01;
    mat2.m_albedo    = color3( 1.f );
    Transform t1;
    t1.scale( 0.7, 0.7, 0.7 );
    t1.translate( vec3( 0.6, 2, -1.8 ) );
    addObject( scene, initSphere( mat2.m_UID, t1 ) );

    // auto& mat4 = scene->CreateMaterial(MaterialModel::COOK_TORRANCE);
    // auto lightIntensity4 =
    //    1.23457 * 4.f * M_PI / ( 4. * M_PI * sqr(0.9) * M_PI );
    // lightIntensity4 = 4;
    // mat4.m_albedo = color3( 1.f );
    // mat4.m_emission = color3(lightIntensity4);
    // Transform t3;
    // t3.scale(0.9,0.9,0.9);
    // t3.translate(vec3(-0.25,5,-3.));
    // addObject(scene, initSphere(mat4, t3));

    // auto obj3   = initSphere( mat4, t3 );
    // auto light4 = new ShapeLight( point3( 3.75, 0, 0 ), color3( lightIntensity4 ), obj3 );
    // addLight(scene, light4);

    // addLight(scene, initPointLight(point3(-0.23, 5, -3), color3(1, 1, 1)));
    ////auto light = new AreaLight(vec3(-0.88,5,-3.57), vec3(0.64,0,0), 4, vec3(0,0,-1), 4,
    /// color3(1,1,1)); /addLight(scene, light);
    return scene;
}

Scene* initScene10() {
    Scene* scene = initScene();
    auto from    = point3( 0., -60, 12 );
    auto at      = vec3( 0, 0, 12 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 0.0, 1 ),
                     30,
                     float( WIDTH ),
                     float( HEIGHT ),
                     0.001,
                     glm::distance( from, at ) );
    setSkyColor( scene, color3( 0., 0., 0. ) );

    auto& mat           = scene->CreateMaterial( MaterialModel::BLINN );
    mat.m_albedo        = color3( 1 );
    mat.m_specularColor = color3( 0 );

    Transform modelMatrix;
    modelMatrix.translate( vec3( 0, 0, 12 ) );
    modelMatrix.translate( vec3( 0, 0, -12 ) );
    modelMatrix.scale( 32, 32, 32 );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix1;
    modelMatrix1.scale( 32, 32, 32 );
    modelMatrix1.rotate( vec3( 1, 0, 0 ), ( 180.f ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    ret            = new Plane( mat.m_UID, modelMatrix1 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix2;
    modelMatrix2.scale( 32, 32, 32 );
    modelMatrix2.rotate( vec3( 1, 0, 0 ), ( 90.f ) );
    modelMatrix2.translate( vec3( 0, 0, 12 ) );
    modelMatrix2.translate( vec3( 0, 20, 0 ) );
    ret            = new Plane( mat.m_UID, modelMatrix2 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1           = scene->CreateMaterial( MaterialModel::BLINN );
    mat1.m_albedo        = color3( 1, 0.5, 0.5 );
    mat1.m_specularColor = color3( 0 );

    Transform modelMatrix3;
    modelMatrix3.scale( 32, 32, 32 );
    modelMatrix3.rotate( vec3( 0, 1, 0 ), ( 90.f ) );
    modelMatrix3.translate( vec3( 0, 0, 12 ) );
    modelMatrix3.translate( vec3( -15, 0, 0 ) );
    ret            = new Plane( mat1.m_UID, modelMatrix3 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat2           = scene->CreateMaterial( MaterialModel::BLINN );
    mat2.m_specularColor = color3( 0 );
    mat2.m_albedo        = color3( 0.5, 0.5, 1.0 );

    Transform modelMatrix4;
    modelMatrix4.scale( 32, 32, 32 );
    modelMatrix4.rotate( vec3( 0, 1, 0 ), ( -90.f ) );
    modelMatrix4.translate( vec3( 0, 0, 12 ) );
    modelMatrix4.translate( vec3( 15, 0, 0 ) );
    ret            = new Plane( mat2.m_UID, modelMatrix4 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat3           = scene->CreateMaterial( MaterialModel::BLINN );
    mat3.m_albedo        = color3( 0.8, 0.2, 0.2 );
    mat3.m_specularColor = color3( 0.7f );
    mat3.m_shininess     = 20;
    mat3.m_reflection    = vec3( 0.7 );

    auto& mat4           = scene->CreateMaterial( MaterialModel::BLINN );
    mat4.m_albedo        = color3( 0.1, 0.1, 0.9 );
    mat4.m_specularColor = color3( 0.9f, 0.9, 1.0 ) * 0.8f;
    mat4.m_IOR           = 1.52;
    mat4.m_shininess     = 10;
    mat4.m_refraction    = vec3( 0.8 );
    mat4.m_absorption    = color3( 0.01, 0.001, 0.0001 );

    Transform modelMatrix5;
    modelMatrix5.scale( 4, 4, 4 );
    modelMatrix5.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix5.translate( vec3( -8, -6, 4 ) );
    addObject( scene, initSphere( mat4.m_UID, modelMatrix5 ) );

    Transform modelMatrix6;
    modelMatrix6.scale( 0.8, 0.8, 0.8 );
    modelMatrix6.rotate( vec3( 0, 0, 1 ), ( -30.f ) );
    modelMatrix6.translate( vec3( 2, 5, 0 ) );
    addObjectsFromFile( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/teapot.obj",
                        scene,
                        mat3,
                        modelMatrix6 );

    Transform modelMatrix7;
    modelMatrix7.scale( 0.3, 0.3, 0.3 );
    modelMatrix7.rotate( vec3( 0, 0, 1 ), ( -60.f ) );
    modelMatrix7.translate( vec3( 5, -6, 0 ) );
    addObjectsFromFile( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/teapot.obj",
                        scene,
                        mat3,
                        modelMatrix7 );

    addLight( scene, initPointLight( point3( 0, 0, 22 ), color3( 0.5f ) ) );
    addLight( scene, initAmbientLight( color3( 0.1 ) ) );
    return scene;
}

Scene* initScene11() {
    Scene* scene = initScene();
    auto from    = point3( 0.0, -70.0, 15.0 );
    auto at      = vec3( 2.0, 0.0, 3.0 );
    // setCameraFOV(scene, from, at, vec3(0, 0, 1), 30.f, float(WIDTH), float(HEIGHT), 0.01,
    // distance(from, at));
    setSimpleCamera(
        scene, from, at, vec3( 0.f, 0.f, 1 ), 30.f, float( WIDTH ), float( HEIGHT ), 0., 1 );
    // setSkyColor( scene, color3( 1., 1., 1 ) );

    // image_texture* sky  = new image_texture( "../assets/clouds.png" );
    // scene->m_skyTexture = sky;
    // Transform texSky;
    // texSky.scale( 1, .4, 1 );
    // texSky.translate( vec3( 0, -0.1, 0 ) );
    // scene->m_skyTexture->m_transform = texSky;

    auto& mat     = scene->CreateMaterial( MaterialModel::BLINN );
    mat.m_albedo  = color3( 1, 1, 1 );
    mat.m_texture = new checker_texture( color3( 0.2, 0.2, 0.2 ), color3( 0.6, 0.6, 0.6 ) );
    Transform texT;
    texT.scale( 0.01, 0.01, 0.01 );
    mat.m_texture->m_transform = texT;
    mat.m_specularColor        = color3( 0 );
    mat.m_reflection           = vec3( 0.5 );

    Transform modelMatrix;
    modelMatrix.scale( 1000, 1000, 1000 );
    modelMatrix.rotate( vec3( 0, 0, 1 ), 30 );

    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );
    // addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

    auto& mat1 = scene->CreateMaterial( MaterialModel::BLINN );
    mat1.m_texture =
        new image_texture( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/bricks.png" );
    mat1.m_specularColor = color3( 0.3 );
    mat1.m_shininess     = 10;

    Transform modelMatrix1;
    modelMatrix1.scale( 0.8, 0.8, 0.8 );
    modelMatrix1.rotate( vec3( 0, 0, 1 ), -50.f );
    modelMatrix1.translate( vec3( 2, -5, 0 ) );
    addObjectsFromFile( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/teapot.obj",
                        scene,
                        mat1,
                        modelMatrix1 );

    auto& mat2     = scene->CreateMaterial( MaterialModel::BLINN );
    mat2.m_texture = new checker_texture( color3( 0.7, 0, 0 ), color3( 0.3, 0, 0 ) );
    Transform texTransform;
    texTransform.scale( 0.25, 0.4, 1 );
    mat2.m_texture->m_transform = texTransform;
    mat2.m_specularColor        = color3( 0.8 );
    mat2.m_shininess            = 100;
    mat2.m_reflection           = vec3( 0.5f );

    Transform modelMatrix2;
    modelMatrix2.scale( 6, 6, 6 );
    modelMatrix2.translate( vec3( 15, 2, 6 ) );

    addObject( scene, initSphere( mat2.m_UID, modelMatrix2 ) );

    auto& mat3           = scene->CreateMaterial( MaterialModel::BLINN );
    mat3.m_albedo        = color3( 0 );
    mat3.m_specularColor = color3( 0.8 );
    mat3.m_shininess     = 100;
    mat3.m_refraction    = vec3( 1.0 );
    mat3.m_IOR           = 1.52;

    Transform modelMatrix3;
    modelMatrix3.scale( 5, 5, 5 );
    modelMatrix3.translate( vec3( -8, -16, 5 ) );

    addObject( scene, initSphere( mat3.m_UID, modelMatrix3 ) );

    addLight( scene, initAmbientLight( color3( 0.2 ) ) );
    addLight( scene, initDirectLight( vec3( -1, 0.2, -1 ), color3( 0.6 ) ) );
    addLight( scene, initDirectLight( vec3( 1, 0.3, -1 ), color3( 0.4 ) ) );

    return scene;
}

Scene* initScene12() {
    Scene* scene = initScene();
    auto from    = point3( 0.0, -70.0, 25.0 );
    auto at      = vec3( -2.0, 0.0, 3.0 );
    // setCameraFOV(scene, from, at, vec3(0, 0, 1), 30.f, float(WIDTH), float(HEIGHT), 0.01,
    // distance(from, at));
    setSimpleCamera(
        scene, from, at, vec3( 0.f, 0.f, 1 ), 25.f, float( WIDTH ), float( HEIGHT ), 1.5, 70 );
    // setSkyColor( scene, color3( 1., 1., 1 ) );

    // image_texture* sky  = new image_texture( "../assets/clouds.png" );
    // scene->m_skyTexture = sky;
    // Transform texSky;
    // texSky.scale( 1, .4, 1 );
    // texSky.translate( vec3( 0, -0.1, 0 ) );
    // scene->m_skyTexture->m_transform = texSky;

    auto& mat           = scene->CreateMaterial( MaterialModel::BLINN );
    mat.m_albedo        = color3( 0.3 );
    mat.m_specularColor = color3( 0.1 );
    mat.m_shininess     = 50;

    Transform modelMatrix;
    modelMatrix.scale( 500, 500, 500 );

    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );
    // addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

    auto& mat1           = scene->CreateMaterial( MaterialModel::BLINN );
    mat1.m_texture       = new image_texture( "../assets/bricks.png" );
    mat1.m_specularColor = color3( 0.3 );
    mat1.m_shininess     = 10;

    Transform modelMatrix1;
    modelMatrix1.scale( 0.7, 0.7, 0.7 );
    modelMatrix1.rotate( vec3( 0, 0, 1 ), -50.f );
    modelMatrix1.translate( vec3( 0, 0, 0 ) );
    addObjectsFromFile( "../assets/teapot.obj", scene, mat1, modelMatrix1 );

    auto& mat2     = scene->CreateMaterial( MaterialModel::BLINN );
    mat2.m_texture = new checker_texture( color3( 0.7, 0, 0 ), color3( 0.3, 0, 0 ) );
    Transform texTransform;
    texTransform.scale( 0.25, 0.4, 1 );
    mat2.m_texture->m_transform = texTransform;
    mat2.m_specularColor        = color3( 0.8 );
    mat2.m_shininess            = 100;
    mat2.m_reflection           = vec3( 0.5f );

    Transform modelMatrix3;
    modelMatrix3.scale( 5, 5, 5 );
    modelMatrix3.translate( vec3( 35, 70, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix3 ) );

    Transform modelMatrix4;
    modelMatrix4.scale( 5, 5, 5 );
    modelMatrix4.translate( vec3( 30, 60, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix4 ) );

    Transform modelMatrix5;
    modelMatrix5.scale( 5, 5, 5 );
    modelMatrix5.translate( vec3( 25, 50, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix5 ) );

    Transform modelMatrix6;
    modelMatrix6.scale( 5, 5, 5 );
    modelMatrix6.translate( vec3( 20, 40, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix6 ) );

    Transform modelMatrix7;
    modelMatrix7.scale( 5, 5, 5 );
    modelMatrix7.translate( vec3( 15, 30, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix7 ) );

    Transform modelMatrix8;
    modelMatrix8.scale( 5, 5, 5 );
    modelMatrix8.translate( vec3( 10, 20, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix8 ) );

    Transform modelMatrix9;
    modelMatrix9.scale( 5, 5, 5 );
    modelMatrix9.translate( vec3( 5, 10, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix9 ) );

    Transform modelMatrix10;
    modelMatrix10.scale( 5, 5, 5 );
    modelMatrix10.translate( vec3( -5, -10, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix10 ) );

    Transform modelMatrix11;
    modelMatrix11.scale( 5, 5, 5 );
    modelMatrix11.translate( vec3( -10, -20, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix11 ) );

    Transform modelMatrix12;
    modelMatrix12.scale( 5, 5, 5 );
    modelMatrix12.translate( vec3( -15, -30, 5 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix12 ) );

    addLight( scene, initAmbientLight( color3( 0.2 ) ) );
    addLight( scene, initDirectLight( vec3( -1, 0.2, -1 ), color3( 0.6 ) ) );
    addLight( scene, initDirectLight( vec3( 1, 0.3, -1 ), color3( 0.4 ) ) );

    return scene;
}

Scene* initScene13() {
    Scene* scene = initScene();
    auto from    = point3( 42, -42.0, 15.0 );
    auto at      = vec3( 6.0, 0.0, -3.0 );
    // setCameraFOV(scene, from, at, vec3(0, 0, 1), 30.f, float(WIDTH), float(HEIGHT), 0.01,
    // distance(from, at));
    setSimpleCamera(
        scene, from, at, vec3( 0.f, 0.f, 1 ), 40.f, float( WIDTH ), float( HEIGHT ), 0, 1 );
    // setSkyColor( scene, color3( 0., 0., 0 ) );

    // image_texture* sky  = new image_texture( "../assets/clouds.png" );
    // scene->m_skyTexture = sky;
    // Transform texSky;
    // texSky.scale( 1, .4, 1 );
    // texSky.translate( vec3( 0, -0.1, 0 ) );
    // scene->m_skyTexture->m_transform = texSky;

    auto& mat     = scene->CreateMaterial( MaterialModel::BLINN );
    mat.m_albedo  = color3( 1, 1, 1 );
    mat.m_texture = new checker_texture( color3( 0.5, 0.5, 0.7 ), color3( 1 ) );
    Transform texCheck;
    texCheck.scale( 0.003, 0.003, 0.003 );
    mat.m_texture->m_transform = texCheck;
    mat.m_specularColor        = color3( 0 );

    Transform modelMatrix;
    modelMatrix.scale( 1000, 1000, 1000 );

    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );
    // addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

    auto& matWall     = scene->CreateMaterial( MaterialModel::BLINN );
    matWall.m_texture = new image_texture( "../assets/bricks.png" );
    Transform texCheckWall;
    texCheckWall.scale( 0.05, 0.05, 0.05 );
    matWall.m_texture->m_transform = texCheckWall;
    matWall.m_specularColor        = color3( 0.3 );
    matWall.m_shininess            = 10;

    Transform modelMatrixWall;
    modelMatrixWall.scale( 100, 100, 100 );
    modelMatrixWall.rotate( vec3( 1, 0, 0 ), 90 );
    modelMatrixWall.translate( vec3( 0, 10, 0 ) );

    auto retWall       = new Plane( matWall.m_UID, modelMatrixWall );
    retWall->geom.type = PLANE;
    addObject( scene, retWall );
    // addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

    auto& mat1           = scene->CreateMaterial( MaterialModel::BLINN );
    mat1.m_albedo        = color3( 0.8, 0.2, 0.2 );
    mat1.m_specularColor = color3( 0.8 );
    mat1.m_shininess     = 100;
    mat1.m_reflection    = vec3( 0.3 );

    Transform modelMatrix1;
    modelMatrix1.scale( 0.5, 0.5, 0.5 );
    modelMatrix1.rotate( vec3( 0, 0, 1 ), -50.f );
    modelMatrix1.translate( vec3( 13, -21, 0 ) );
    addObjectsFromFile( "../assets/teapot.obj", scene, mat1, modelMatrix1 );

    auto& mat2           = scene->CreateMaterial( MaterialModel::BLINN );
    mat2.m_albedo        = color3( 0 );
    mat2.m_specularColor = color3( 0.8 );
    mat2.m_shininess     = 100;
    mat2.m_reflection    = vec3( 0.8f );

    Transform modelMatrix3;
    modelMatrix3.scale( 6, 6, 6 );
    modelMatrix3.translate( vec3( -28, 0, 6 ) );
    addObject( scene, initSphere( mat2.m_UID, modelMatrix3 ) );

    auto& mat3           = scene->CreateMaterial( MaterialModel::BLINN );
    mat3.m_albedo        = color3( 0 );
    mat3.m_specularColor = color3( 0.8 );
    mat3.m_shininess     = 50;
    mat3.m_reflection    = vec3( 0.8f );
    // mat3.m_reflectionGloss = 0.05;

    Transform modelMatrix4;
    modelMatrix4.scale( 6, 6, 6 );
    modelMatrix4.translate( vec3( -14, 0, 6 ) );
    addObject( scene, initSphere( mat3.m_UID, modelMatrix4 ) );

    auto& mat4           = scene->CreateMaterial( MaterialModel::BLINN );
    mat4.m_albedo        = color3( 0 );
    mat4.m_specularColor = color3( 0.8 );
    mat4.m_shininess     = 20;
    mat4.m_reflection    = vec3( 0.8f );
    // mat4.m_reflectionGloss = 0.1;

    Transform modelMatrix5;
    modelMatrix5.scale( 6, 6, 6 );
    modelMatrix5.translate( vec3( 0, 0, 6 ) );
    addObject( scene, initSphere( mat4.m_UID, modelMatrix5 ) );

    auto& mat5           = scene->CreateMaterial( MaterialModel::BLINN );
    mat5.m_albedo        = color3( 0 );
    mat5.m_specularColor = color3( 0.8 );
    mat5.m_shininess     = 10;
    mat5.m_reflection    = vec3( 0.8f );
    // mat5.m_reflectionGloss = 0.2;

    Transform modelMatrix6;
    modelMatrix6.scale( 6, 6, 6 );
    modelMatrix6.translate( vec3( 14, 0, 6 ) );
    addObject( scene, initSphere( mat5.m_UID, modelMatrix6 ) );

    auto& mat6           = scene->CreateMaterial( MaterialModel::BLINN );
    mat6.m_albedo        = color3( 0 );
    mat6.m_specularColor = color3( 0.8 );
    mat6.m_shininess     = 50;
    mat6.m_refraction    = vec3( 0.8 );
    mat6.m_IOR           = 2.0;
    // mat6.m_refractionGloss = 0.05;

    Transform modelMatrix7;
    modelMatrix7.scale( 2, 2, 2 );
    modelMatrix7.translate( vec3( 22, -30, 2 ) );
    addObject( scene, initSphere( mat6.m_UID, modelMatrix7 ) );

    addLight( scene, initAmbientLight( color3( 0.2 ) ) );

    addLight( scene, initPointLight( vec3( -50, -100, 50 ), vec3( 1, 1, 1 ), 5 ) );

    return scene;
}

Scene* initScene14() {
    Scene* scene = initScene();
    auto from    = point3( 0., -60, 16 );
    auto at      = vec3( 0, 0, 11 );
    setSimpleCamera(
        scene, from, at, vec3( 0, 0.0, 1 ), 30, float( WIDTH ), float( HEIGHT ), 0, 1 );
    setSkyColor( scene, color3( 0., 0., 0. ) );

    auto& mat           = scene->CreateMaterial( MaterialModel::BLINN );
    mat.m_albedo        = color3( 1 );
    mat.m_specularColor = color3( 0 );

    Transform modelMatrix;
    modelMatrix.scale( 32, 32, 32 );
    modelMatrix.translate( vec3( 0, 0, 12 ) );
    modelMatrix.translate( vec3( 0, 0, -12 ) );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix1;
    modelMatrix1.scale( 32, 32, 32 );
    modelMatrix1.rotate( vec3( 1, 0, 0 ), ( 180.f ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    ret            = new Plane( mat.m_UID, modelMatrix1 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix2;
    modelMatrix2.scale( 32, 32, 32 );
    modelMatrix2.rotate( vec3( 1, 0, 0 ), ( 90.f ) );
    modelMatrix2.translate( vec3( 0, 0, 12 ) );
    modelMatrix2.translate( vec3( 0, 20, 0 ) );
    ret            = new Plane( mat.m_UID, modelMatrix2 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1           = scene->CreateMaterial( MaterialModel::BLINN );
    mat1.m_albedo        = color3( 1, 0.2, 0.2 );
    mat1.m_specularColor = color3( 0 );

    Transform modelMatrix3;
    modelMatrix3.scale( 32, 32, 32 );
    modelMatrix3.rotate( vec3( 0, 1, 0 ), ( 90.f ) );
    modelMatrix3.translate( vec3( 0, 0, 12 ) );
    modelMatrix3.translate( vec3( -15, 0, 0 ) );
    ret            = new Plane( mat1.m_UID, modelMatrix3 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat2           = scene->CreateMaterial( MaterialModel::BLINN );
    mat2.m_specularColor = color3( 0 );
    mat2.m_albedo        = color3( 0.2, 0.2, 1.0 );

    Transform modelMatrix4;
    modelMatrix4.scale( 32, 32, 32 );
    modelMatrix4.rotate( vec3( 0, 1, 0 ), ( -90.f ) );
    modelMatrix4.translate( vec3( 0, 0, 12 ) );
    modelMatrix4.translate( vec3( 15, 0, 0 ) );
    ret            = new Plane( mat2.m_UID, modelMatrix4 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat3           = scene->CreateMaterial( MaterialModel::BLINN );
    mat3.m_albedo        = color3( 0.8, 0.2, 0.2 );
    mat3.m_specularColor = color3( 0.7f );
    mat3.m_shininess     = 20;

    auto& mat4           = scene->CreateMaterial( MaterialModel::BLINN );
    mat4.m_albedo        = color3( 0.2, 0.8, 0.2 );
    mat4.m_specularColor = color3( 0.9f, 0.9, 1.0 ) * 0.8f;
    mat4.m_shininess     = 10;

    Transform modelMatrix5;
    modelMatrix5.scale( 4.5, 4.5, 4.5 );
    modelMatrix5.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix5.translate( vec3( 8, -6, 4.5 ) );
    addObject( scene, initSphere( mat4.m_UID, modelMatrix5 ) );

    Transform modelMatrix6;
    modelMatrix6.scale( 0.75, 0.75, 0.75 );
    modelMatrix6.rotate( vec3( 0, 0, 1 ), ( 30.f ) );
    modelMatrix6.translate( vec3( -4.5, 5, 0 ) );
    addObjectsFromFile( "../assets/teapot.obj", scene, mat3, modelMatrix6 );

    // auto& matE = scene->CreateMaterial(MaterialModel::BLINN);
    // matE->m_emission = color3(0.5);
    // matE->m_diffuseColor = color3(0.5f);

    // Transform modemMatrixE;
    // modemMatrixE.scale(5.f, 5.f, 5.f);
    // modemMatrixE.translate(vec3(0,0,19));
    // addObject(scene, initSphere(matE, modemMatrixE));

    addLight( scene, initPointLight( point3( 0, 0, 19 ), color3( 0.5f ), 5.0 ) );
    return scene;
}

Scene* initScene15() {
    Scene* scene = initScene();
    auto from    = point3( 0., -60, 16 );
    auto at      = vec3( 0, 0, 11 );
    setSimpleCamera(
        scene, from, at, vec3( 0, 0.0, 1 ), 30, float( WIDTH ), float( HEIGHT ), 0, 1 );
    setSkyColor( scene, color3( 0., 0., 0. ) );

    scene->medium       = new Fog( 0.04, false );
    auto& mat           = scene->CreateMaterial( MaterialModel::BLINN );
    mat.m_albedo        = color3( 0.3 );
    mat.m_specularColor = color3( 0 );

    Transform modelMatrix;
    modelMatrix.scale( 32, 32, 32 );
    modelMatrix.translate( vec3( 0, 0, 12 ) );
    modelMatrix.translate( vec3( 0, 0, -12 ) );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix1;
    modelMatrix1.scale( 32, 32, 32 );
    modelMatrix1.rotate( vec3( 1, 0, 0 ), ( 180.f ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    ret            = new Plane( mat.m_UID, modelMatrix1 );
    ret->geom.type = PLANE;
    // addObject(scene, ret);

    auto& matBehind = scene->CreateMaterial( MaterialModel::BLINN );
    // matBehind->m_diffuseColor = color3(1.0, 1., 1.0);
    auto walltext = new image_texture( "../assets/walltext.png" );
    Transform walltextT;
    // walltextT.scale(0.5, 0.5, 0.5);
    walltextT.translate( vec3( 0, -.19, 0 ) );
    walltext->m_transform = walltextT;
    matBehind.m_texture   = walltext;
    Transform modelMatrix2;
    modelMatrix2.scale( 32, 32, 32 );
    // modelMatrix2.scale(21,21,21);
    modelMatrix2.rotate( vec3( 1, 0, 0 ), ( 90.f ) );
    modelMatrix2.translate( vec3( 0, 0, 12 ) );
    modelMatrix2.translate( vec3( 0, 20, 0 ) );
    ret            = new Plane( mat.m_UID, modelMatrix2 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1           = scene->CreateMaterial( MaterialModel::BLINN );
    mat1.m_albedo        = color3( 1, 0.2, 0.2 );
    mat1.m_specularColor = color3( 0 );

    Transform modelMatrix3;
    modelMatrix3.scale( 32, 32, 32 );
    modelMatrix3.rotate( vec3( 0, 1, 0 ), ( 90.f ) );
    modelMatrix3.translate( vec3( 0, 0, 12 ) );
    modelMatrix3.translate( vec3( -15, 0, 0 ) );
    ret            = new Plane( mat1.m_UID, modelMatrix3 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat2           = scene->CreateMaterial( MaterialModel::BLINN );
    mat2.m_specularColor = color3( 0 );
    mat2.m_albedo        = color3( 0.2, 0.2, 1.0 );

    Transform modelMatrix4;
    modelMatrix4.scale( 32, 32, 32 );
    modelMatrix4.rotate( vec3( 0, 1, 0 ), ( -90.f ) );
    modelMatrix4.translate( vec3( 0, 0, 12 ) );
    modelMatrix4.translate( vec3( 15, 0, 0 ) );
    ret            = new Plane( mat2.m_UID, modelMatrix4 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat3           = scene->CreateMaterial( MaterialModel::BLINN );
    mat3.m_albedo        = color3( 0.8, 0.2, 0.2 );
    mat3.m_specularColor = color3( 0 );

    auto& mat4 = scene->CreateMaterial( MaterialModel::BLINN );
    // mat4.m_diffuseColor = color3(0.2, 0.8, 0.2);
    // mat4.m_diffuseColor = color3(1);
    // mat4.m_specularColor = color3(1);
    // mat4.m_roughness = 0.1;
    mat4.m_refraction = color3( 1 );
    // mat4.m_reflection = color3(1);
    mat4.m_IOR = 1.5;

    // auto& mat3 = std::make_shared<Blinn>(true);
    // mat3.m_diffuseColor = color3(0.8, 0.2, 0.2);
    // mat3.m_specularColor = color3(0.7f);
    // mat3.m_IOR = 1.15;
    // mat3.m_roughness = 0.1;

    // auto& mat4 = scene->CreateMaterial(MaterialModel::BLINN);
    // mat4.m_diffuseColor = color3(0.0, 1., 0.0);
    // mat4.m_specularColor = color3(0.9f, 0.9, 1.0) * 0.8f;
    // mat4.m_IOR = 1.8;
    // mat4.m_roughness = 0.1;

    // auto& mat4 = scene->CreateMaterial(MaterialModel::BLINN);
    ////mat4.m_diffuseColor = color3(0.0, 1., 0.0);
    // auto balltex =  new image_texture("../assets/6ball.png");
    ////Transform balltexT;
    ////balltexT.rotate(vec3(0,0,1), 60);
    ////balltex->m_transform = balltexT;
    // mat4.m_texture = balltex;
    // mat4.m_specularColor = color3(0.9f, 0.9, 1.0) * 0.8f;
    // mat4.m_IOR = 4.0;
    // mat4.m_roughness = 0.04;

    Transform modelMatrix5;
    modelMatrix5.scale( 4.5, 4.5, 4.5 );
    modelMatrix5.rotate( vec3( 0, 0, 1 ), ( 150.f ) );
    modelMatrix5.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix5.translate( vec3( 0, -2, 13.5 ) );
    addObject( scene, initSphere( mat4.m_UID, modelMatrix5 ) );

    Transform modelMatrix7;
    modelMatrix7.scale( 4.5, 4.5, 4.5 );
    modelMatrix7.rotate( vec3( 0, 0, 1 ), ( 150.f ) );
    modelMatrix7.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix7.translate( vec3( -8, -6, 8.5 ) );
    addObject( scene, initSphere( mat3.m_UID, modelMatrix7 ) );

    Transform modelMatrix6;
    modelMatrix6.scale( 0.75, 0.75, 0.75 );
    modelMatrix6.rotate( vec3( 0, 0, 1 ), ( 30.f ) );
    modelMatrix6.translate( vec3( -4.5, 5, 0 ) );
    // addObjectsFromFile("../assets/teapot.obj", scene, mat3, modelMatrix6);

    auto lightSize      = 4.f;
    auto lightIntensity = 99000.f * 4.f * M_PI / ( 4. * M_PI * lightSize * lightSize * M_PI );
    auto& matE          = scene->CreateMaterial( MaterialModel::BLINN );
    matE.m_emission     = color3( lightIntensity );
    matE.m_albedo       = color3( 1.f );

    Transform modemMatrixE;
    // modemMatrixE.scale(lightSize, lightSize, lightSize);
    // modemMatrixE.translate(vec3(-2,0,24.5));
    addObject( scene, initSphere( matE.m_UID, modemMatrixE ) );

    // addLight(scene, initPointLight(point3(0, 0, 17), color3(0.5f), lightSize));
    auto obj   = initSphere( matE.m_UID, modemMatrixE );
    auto light = new ShapeLight( point3( 0, 0, 17 ), color3( 2000 ), obj );
    // addLight(scene, light);
    return scene;
}

Scene* initScene16() {
    Scene* scene = initScene();
    auto from    = point3( 0., -60, 16 );
    auto at      = vec3( 0, 0, 11 );
    setSimpleCamera(
        scene, from, at, vec3( 0, 0.0, 1 ), 30, float( WIDTH ), float( HEIGHT ), 0, 1 );

    scene->sky = new UniformSky( vec3( 0 ) );

    auto& mat       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat.m_albedo    = color3( 1 );
    mat.m_metalness = 0.0;
    mat.m_roughness = 1.0;

    Transform modelMatrix;
    modelMatrix.scale( 32, 32, 32 );
    modelMatrix.translate( vec3( 0, 0, 12 ) );
    modelMatrix.translate( vec3( 0, 0, -12 ) );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix1;
    modelMatrix1.scale( 32, 32, 32 );
    modelMatrix1.rotate( vec3( 1, 0, 0 ), ( 180.f ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    ret            = new Plane( mat.m_UID, modelMatrix1 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix2;
    modelMatrix2.scale( 32, 32, 32 );
    modelMatrix2.rotate( vec3( 1, 0, 0 ), ( 90.f ) );
    modelMatrix2.translate( vec3( 0, 0, 12 ) );
    modelMatrix2.translate( vec3( 0, 20, 0 ) );
    ret            = new Plane( mat.m_UID, modelMatrix2 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat1.m_albedo    = color3( 1, 0.2, 0.2 );
    mat1.m_metalness = 0.0;
    mat1.m_roughness = 1.0;

    Transform modelMatrix3;
    modelMatrix3.scale( 32, 32, 32 );
    modelMatrix3.rotate( vec3( 0, 1, 0 ), ( 90.f ) );
    modelMatrix3.translate( vec3( 0, 0, 12 ) );
    modelMatrix3.translate( vec3( -15, 0, 0 ) );
    ret            = new Plane( mat1.m_UID, modelMatrix3 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat2       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat2.m_albedo    = color3( 0.8, 0.2, 0.2 );
    mat2.m_metalness = 0.0;
    mat2.m_roughness = 1.0;

    Transform modelMatrix4;
    modelMatrix4.scale( 32, 32, 32 );
    modelMatrix4.rotate( vec3( 0, 1, 0 ), ( -90.f ) );
    modelMatrix4.translate( vec3( 0, 0, 12 ) );
    modelMatrix4.translate( vec3( 15, 0, 0 ) );
    ret            = new Plane( mat2.m_UID, modelMatrix4 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat3       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat3.m_albedo    = color3( 0.8, 0.2, 0.2 );
    mat3.m_metalness = 0.0;
    mat3.m_roughness = 1.0;

    auto& mat4       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat4.m_albedo    = color3( 0.0, 1., 0.0 );
    mat4.m_metalness = 0.0;
    mat4.m_roughness = 0.4;

    Transform modelMatrix5;
    modelMatrix5.scale( 4.5, 4.5, 4.5 );
    modelMatrix5.rotate( vec3( 0, 0, 1 ), ( 150.f ) );
    modelMatrix5.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix5.translate( vec3( 8, -6, 4.5 ) );
    addObject( scene, initSphere( mat4.m_UID, modelMatrix5 ) );

    Transform modelMatrix6;
    modelMatrix6.scale( 0.75, 0.75, 0.75 );
    modelMatrix6.rotate( vec3( 0, 0, 1 ), ( 30.f ) );
    modelMatrix6.translate( vec3( -4.5, 5, 0 ) );
    // addObjectsFromFile("../assets/teapot.obj", scene, mat3, modelMatrix6);

    // auto& matE = scene->CreateMaterial(MaterialModel::COOK_TORRANCE);
    // matE->m_emission = color3(0.5);
    // matE->m_diffuseColor = color3(0.5f);

    // Transform modemMatrixE;
    // modemMatrixE.scale(5.f, 5.f, 5.f);
    // modemMatrixE.translate(vec3(0,0,19));
    // addObject(scene, initSphere(matE, modemMatrixE));

    // addLight( scene, initPointLight( point3( 0, 0, 19 ), color3( 0.5f ), 5.0 ) );
    //
    auto lightSize = 3.f;
    auto lightIntensity =
        /*99000.f*/ 10000 * 4.f * M_PI / ( 4. * M_PI * lightSize * lightSize * M_PI );
    lightIntensity  = 10.0;
    auto& matE      = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    matE.m_emission = color3( lightIntensity );
    matE.m_albedo   = color3( 1.f );

    Transform modemMatrixE;
    modemMatrixE.scale( lightSize, lightSize, lightSize );
    modemMatrixE.translate( vec3( 0, 0, 19 ) );
    addObject( scene, initSphere( matE.m_UID, modemMatrixE ) );

    // addLight(scene, initPointLight(point3(0, 0, 17), color3(0.5f), lightSize));
    auto obj   = initSphere( matE.m_UID, modemMatrixE );
    auto light = new ShapeLight( point3( 0, 0, 19 ), color3( lightIntensity ), obj );
    addLight( scene, light );
    return scene;
}

Scene* initScene17() {
    Scene* scene = initScene();
    auto from    = point3( 0., -60, 16 );
    auto at      = vec3( 0, 0, 11 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 0.0, 1 ),
                     30,
                     float( WIDTH ),
                     float( HEIGHT ),
                     0.001,
                     distance( from, at ) );

    scene->sky = new UniformSky( vec3( 0 ) );

    auto& mat       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat.m_albedo    = color3( 0.8 );
    mat.m_metalness = 0.0;
    mat.m_roughness = 1.0;

    Transform modelMatrix;
    modelMatrix.scale( 32, 32, 32 );
    modelMatrix.translate( vec3( 0, 0, 12 ) );
    modelMatrix.translate( vec3( 0, 0, -12 ) );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    Transform modelMatrix1;
    modelMatrix1.scale( 32, 32, 32 );
    modelMatrix1.rotate( vec3( 1, 0, 0 ), ( 180.f ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    modelMatrix1.translate( vec3( 0, 0, 12 ) );
    ret            = new Plane( mat.m_UID, modelMatrix1 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& matBehind = scene->CreateMaterial( MaterialModel::BLINN );
    auto walltext =
        new image_texture( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/walltext.png" );
    Transform walltextT;
    walltextT.translate( vec3( 0, -.19, 0 ) );
    walltext->m_transform = walltextT;
    matBehind.m_texture   = walltext;
    Transform modelMatrix2;
    modelMatrix2.scale( 32, 32, 32 );
    // modelMatrix2.scale(21,21,21);
    modelMatrix2.rotate( vec3( 1, 0, 0 ), ( 90.f ) );
    modelMatrix2.translate( vec3( 0, 0, 12 ) );
    modelMatrix2.translate( vec3( 0, 20, 0 ) );
    ret            = new Plane( mat.m_UID, modelMatrix2 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat1.m_albedo    = color3( 1, 0.2, 0.2 );
    mat1.m_metalness = 0.0;
    mat1.m_roughness = 1.0;

    Transform modelMatrix3;
    modelMatrix3.scale( 32, 32, 32 );
    modelMatrix3.rotate( vec3( 0, 1, 0 ), ( 90.f ) );
    modelMatrix3.translate( vec3( 0, 0, 12 ) );
    modelMatrix3.translate( vec3( -15, 0, 0 ) );
    ret            = new Plane( mat1.m_UID, modelMatrix3 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat2       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat2.m_albedo    = color3( 0.2, 0.2, 1.0 );
    mat2.m_metalness = 0.0;
    mat2.m_roughness = 1.0;

    Transform modelMatrix4;
    modelMatrix4.scale( 32, 32, 32 );
    modelMatrix4.rotate( vec3( 0, 1, 0 ), ( -90.f ) );
    modelMatrix4.translate( vec3( 0, 0, 12 ) );
    modelMatrix4.translate( vec3( 15, 0, 0 ) );
    ret            = new Plane( mat2.m_UID, modelMatrix4 );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat3       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat3.m_albedo    = color3( 0.6, 0.2, 0.2 );
    mat3.m_roughness = 0.5;
    mat3.m_IOR       = 1.01;

    auto& mat4    = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, TRANSPARENT );
    mat4.m_albedo = color3( 1 );
    mat4.m_IOR    = 1.5;

    Transform modelMatrix5;
    modelMatrix5.scale( 4.5, 4.5, 4.5 );
    modelMatrix5.rotate( vec3( 0, 0, 1 ), ( 150.f ) );
    modelMatrix5.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix5.translate( vec3( 0, -2, 13.5 ) );
    addObject( scene, initSphere( mat4.m_UID, modelMatrix5 ) );

    Transform modelMatrix7;
    modelMatrix7.scale( 4.5, 4.5, 4.5 );
    modelMatrix7.rotate( vec3( 0, 0, 1 ), ( 150.f ) );
    modelMatrix7.rotate( vec3( 0, 1, 0 ), ( 30.f ) );
    modelMatrix7.translate( vec3( -8, -6, 8.5 ) );
    addObject( scene, initSphere( mat3.m_UID, modelMatrix7 ) );

    Transform modelMatrix6;
    modelMatrix6.scale( 0.75, 0.75, 0.75 );
    modelMatrix6.rotate( vec3( 0, 0, 1 ), ( 30.f ) );
    modelMatrix6.translate( vec3( -4.5, 5, 0 ) );
    // addObjectsFromFile("../assets/teapot.obj", scene, mat3, modelMatrix6);

    auto lightSize = 0.9f;
    auto lightIntensity =
        /*99000.f*/ 10000 * 4.f * M_PI / ( 4. * M_PI * lightSize * lightSize * M_PI );
    lightIntensity  = 4.0;
    auto& matE      = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    matE.m_emission = color3( lightIntensity );
    matE.m_albedo   = color3( 1.f );

    Transform modemMatrixE;
    modemMatrixE.scale( lightSize, lightSize, lightSize );
    modemMatrixE.translate( vec3( 0, 0, 5 ) );
    addObject( scene, initSphere( matE.m_UID, modemMatrixE ) );

    // addLight(scene, initPointLight(point3(0, 0, 17), color3(0.5f), lightSize));
    auto obj   = initSphere( matE.m_UID, modemMatrixE );
    auto light = new ShapeLight( point3( 0, 0, 5 ), color3( lightIntensity ), obj );
    addLight( scene, light );
    return scene;
}

Scene* initScene18() {
    Scene* scene = initScene();
    auto from    = point3( 1.9166, 0.4598, 1.1936 );
    auto at      = vec3( 0.7520, 0.33266, 0.4188 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 1, 0 ),
                     80,
                     (float)WIDTH,
                     (float)HEIGHT,
                     0.04,
                     glm::length( at - from ) * 1.5 );

    // auto sky =
    //     new IBL( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/rainforest_trail_4k.hdr" );
    auto sky =
        new IBL( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/brown_photostudio_06_4k.hdr" );

    addLight( scene, sky );
    scene->envLights.push_back( sky );

    // auto& matLight            = scene->CreateMaterial(MaterialModel::COOK_TORRANCE);
    // auto lightIntensity1 = 100;
    // matLight->m_albedo       = color3( 1.f );
    // matLight->m_emission     = color3( lightIntensity1 );
    // Transform tLight;
    // tLight.scale( 0.2, 0.2, 0.2 );
    // tLight.translate( vec3( -1.25, 1.5, 0 ) );
    // addObject( scene, initSphere( matLight, tLight ) );

    // auto obj0   = initSphere( matLight, tLight );
    // auto lightLight = new ShapeLight( point3( 0, 0, 0 ), color3( lightIntensity1 ), obj0 );
    // addLight( scene, lightLight );

    auto& mat       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat.m_albedo    = color3( 0.3 );
    mat.m_metalness = 0.1;
    mat.m_roughness = 0.64;

    Transform modelMatrix;
    modelMatrix.scale( 5, 1.3, 1 );
    modelMatrix.rotate( vec3( -1, 0, 0 ), 90 );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat1.m_roughness = 0.004;
    mat1.m_metalness = 1.0;
    mat1.m_albedo    = color3( 0.8 );
    Transform t0;
    t0.scale( 0.5, 0.5, 0.5 );
    t0.translate( vec3( 0, 0.5, 0 ) );
    addObject( scene, initSphere( mat1.m_UID, t0 ) );

    auto& mat2       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, TRANSPARENT );
    mat2.m_IOR       = 1.5;
    mat2.m_roughness = 0.01;
    mat2.m_albedo    = color3( 1.f );
    Transform t1;
    t1.scale( 0.5, 0.5, 0.5 );
    t1.translate( vec3( -1.2, 0.5, 0 ) );
    addObject( scene, initSphere( mat2.m_UID, t1 ) );

    auto& mat3       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat3.m_roughness = 0.004;
    mat3.m_metalness = 0.0000;
    mat3.m_albedo    = color3( 0.f, 1.f, 0.95f );
    // mat3.m_albedo = color3(0.5);
    Transform t2;
    t2.scale( 0.5, 0.5, 0.5 );
    t2.translate( vec3( 1.2, 0.5, 0 ) );
    addObject( scene, initSphere( mat3.m_UID, t2 ) );

    auto& mat4       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat4.m_roughness = 0.1;
    mat4.m_metalness = 0.01;
    mat4.m_albedo    = color3( 0, 0.f, 0 );
    // mat4.m_albedo = color3(0.5);
    mat4.m_texture = new image_texture(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/pf-s125-ake7011-a.png" );
    Transform t3;
    t3.scale( 0.5, 0.5, 0.5 );
    t3.rotate( vec3( 1, 0, 0 ), 90 );
    t3.translate( vec3( 2.4, 0.5, 0 ) );
    addObject( scene, initSphere( mat4.m_UID, t3 ) );

    auto& mat5       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat5.m_roughness = 0.2;
    mat5.m_metalness = 0.9;
    mat5.m_albedo    = color3( 1.f );
    // mat5.m_albedo = color3(0.5);
    Transform t4;
    t4.scale( 0.5, 0.5, 0.5 );
    t4.translate( vec3( -2.4, 0.5, 0 ) );
    addObject( scene, initSphere( mat5.m_UID, t4 ) );

    auto lightSize      = 0.01f;
    auto lightIntensity = 0.012 * 4.f * M_PI / ( 4. * M_PI * lightSize * lightSize * M_PI );
    auto& matE          = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    matE.m_emission     = color3( lightIntensity );
    matE.m_albedo       = color3( 1.f );

    Transform modemMatrixE;
    modemMatrixE.scale( lightSize, lightSize, lightSize );
    modemMatrixE.translate( vec3( 0, 5.0, 5 ) );
    // addObject( scene, initSphere( matE, modemMatrixE ) );

    return scene;
}

Scene* initScene19() {
    Scene* scene = initScene();
    auto from    = point3( -3.5, -2, 7.2 );
    auto at      = vec3( 3.0, -2.4, 2.6 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 1, 0 ),
                     80,
                     (float)WIDTH,
                     (float)HEIGHT,
                     0.04,
                     glm::length( at - from ) );

    setSkyColor( scene, color3( 1.f ) );
    auto sky =
        new IBL( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/brown_photostudio_06_4k.hdr" );

    addLight( scene, sky );
    scene->envLights.push_back( sky );

    // auto left  = new image_texture( "../assets/Standard-Cube-Map1/StandardCubeMap/left.png" );
    // auto right = new image_texture( "../assets/Standard-Cube-Map1/StandardCubeMap/right.png" );
    // auto front = new image_texture( "../assets/Standard-Cube-Map1/StandardCubeMap/front.png" );
    // auto back  = new image_texture( "../assets/Standard-Cube-Map1/StandardCubeMap/back.png" );
    // auto up    = new image_texture( "../assets/Standard-Cube-Map1/StandardCubeMap/top.png" );
    // auto down  = new image_texture( "../assets/Standard-Cube-Map1/StandardCubeMap/bottom.png" );

    // auto& matSky       = scene->CreateMaterial(MaterialModel::COOK_TORRANCE);
    // matSky->m_texture = new CubeMapTexture( left, right, front, back, up, down );
    // matSky->m_albedo  = color3( 0, 0, 0 );
    // Transform skyboxT;
    // skyboxT.scale( 50, 50, 50 );
    // addObject( scene, initCube( matSky, skyboxT ) );

    float metalness = 0.01;
    float y         = -6.4;
    for ( int i = 0; i < 7; i++ ) {
        float x         = -2;
        float roughness = 0.01;
        for ( int j = 0; j < 7; j++ ) {
            auto& mat3       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
            mat3.m_roughness = roughness;
            mat3.m_metalness = metalness;
            mat3.m_albedo    = color3( 0.9f, 0, 0 );
            Transform t2;
            t2.scale( 0.5, 0.5, 0.5 );
            t2.translate( vec3( x, y, 0 ) );
            addObject( scene, initSphere( mat3.m_UID, t2 ) );
            roughness += 0.14;
            x += 1.2;
        }
        metalness += 0.14;
        y += 1.2;
        std::cout << "rough: " << roughness << std::endl;
    }
    std::cout << "metal: " << metalness << std::endl;

    metalness = 0.01;
    y         = -6.4;
    for ( int i = 0; i < 7; i++ ) {
        float x         = 8.6;
        float roughness = 0.01;
        for ( int j = 0; j < 7; j++ ) {
            auto& mat3       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
            mat3.m_roughness = roughness;
            mat3.m_metalness = metalness;
            mat3.m_albedo    = color3( 0.8, 0.1, 0.f );
            Transform t2;
            t2.scale( 0.5, 0.5, 0.5 );
            t2.translate( vec3( 6, y, x ) );
            addObject( scene, initSphere( mat3.m_UID, t2 ) );
            roughness += 0.14;
            x -= 1.2;
        }
        metalness += 0.14;
        y += 1.2;
        std::cout << "rough: " << roughness << std::endl;
    }
    std::cout << "metal: " << metalness << std::endl;

    auto lightSize      = 0.01f;
    auto lightIntensity = 0.09 * 4.f * M_PI / ( 4. * M_PI * lightSize * lightSize * M_PI );
    auto& matE          = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    matE.m_emission     = color3( lightIntensity );
    matE.m_albedo       = color3( 1.f );

    Transform modemMatrixE;
    modemMatrixE.scale( lightSize, lightSize, lightSize );
    modemMatrixE.translate( vec3( -8, -1.5, 15 ) );
    // addObject( scene, initSphere( matE, modemMatrixE ) );

    return scene;
}

Scene* initScene20() {
    Scene* scene = initScene();
    auto from    = point3( 0, 0.9, 3 );
    auto at      = vec3( 0, 0.5, 0 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 1, 0 ),
                     90,
                     (float)WIDTH,
                     (float)HEIGHT,
                     0.02,
                     glm::length( at - from ) );

    auto sky =
        new IBL( "C:/Users/marco/Documents/GitRepos/Raytracer/assets/rainforest_trail_4k.hdr" );
    // auto sky = new IBL( "/home/mafo/dev/Raytracer/assets/spaichingen_hill_4k.hdr" );
    Transform skyT;
    //// skyT.rotate(vec3(0,1,0), 80);
    // sky->m_transform = skyT;
    //  scene->sky       = sky;
    addLight( scene, sky );
    scene->envLights.push_back( sky );

    auto& mat       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat.m_albedo    = color3( 0.6 );
    mat.m_metalness = 0.0;
    mat.m_roughness = 0.8;
    // auto& mat             = scene->CreateMaterial(MaterialModel::COOK_TORRANCE, TRANSPARENT);
    // mat.m_albedo  = color3( 0.6 );
    // mat.m_IOR = 1.1;
    // mat.m_metalness = 0.0;
    // mat.m_roughness = 0.8;

    Transform modelMatrix;
    modelMatrix.scale( 5, 1.3, 1 );
    modelMatrix.rotate( vec3( -1, 0, 0 ), 90 );
    auto ret       = new Plane( mat.m_UID, modelMatrix );
    ret->geom.type = PLANE;
    addObject( scene, ret );

    auto& mat1       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat1.m_roughness = 0.01;
    mat1.m_metalness = 1.0;
    mat1.m_albedo    = color3( 1.f );
    Transform t0;
    t0.scale( 0.75, 0.75, 0.75 );
    t0.translate( vec3( 0, 0., 0 ) );
    // addObjectsFromFile("../assets/bunny.obj", scene, mat1, t0);

    auto& mat2       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, TRANSPARENT );
    mat2.m_IOR       = 1.5;
    mat2.m_roughness = 0.01;
    mat2.m_albedo    = color3( 1.f );
    Transform t1;
    t1.scale( 0.75, 0.75, 0.75 );
    t1.translate( vec3( 1.2, 0., 0 ) );
    // addObjectsFromFile("../assets/bunny.obj", scene, mat2, t1);

    auto& mat3       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat3.m_roughness = 0.0001;
    mat3.m_metalness = 0.0000;
    mat3.m_albedo    = color3( 0.f, 1.f, 0.95f );
    Transform t2;
    t2.scale( 0.75, 0.75, 0.75 );
    t2.translate( vec3( -1.2, 0., 0 ) );
    // addObjectsFromFile("../assets/bunny.obj", scene, mat3, t2);

    auto& mat4       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat4.m_roughness = 0.5;
    mat4.m_metalness = 0.1;
    mat4.m_albedo    = color3( 0, 1.f, 0 );
    Transform t3;
    t3.scale( 0.75, 0.75, 0.75 );
    t3.translate( vec3( -2.4, 0., 0 ) );
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/bunny.obj", scene, mat4, t3 );

    auto& mat5       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, SPECULAR );
    mat5.m_roughness = 0.2;
    mat5.m_metalness = 0.9;
    mat5.m_albedo    = color3( 1.f );
    Transform t4;
    t4.scale( 0.75, 0.75, 0.75 );
    t4.translate( vec3( 2.4, 0., 0 ) );
    // addObjectsFromFile("../assets/bunny.obj", scene, mat5, t4);

    return scene;
}

Scene* initScene21() {
    Scene* scene = initScene();
    auto from    = point3( 0, 2, 15 );
    auto at      = vec3( 0, -2, 2.15 );
    setSimpleCamera( scene,
                     from,
                     at,
                     vec3( 0, 1, 0 ),
                     28,
                     (float)WIDTH,
                     (float)HEIGHT,
                     0.001,
                     glm::length( at - from ) );

    // auto sky = new IBL("/home/mafo/dev/Raytracer/assets/rainforest_trail_4k.hdr");
    // auto sky = new IBL("/home/mafo/dev/Raytracer/assets/spaichingen_hill_4k.hdr");
    auto sky = new UniformSky( vec3( 0 ) );
    Transform skyT;
    // skyT.rotate(vec3(0,1,0), 80);
    // sky->m_transform = skyT;
    scene->sky = sky;

    auto lightIntensity = 800 * 4.f * M_PI / ( 4. * M_PI * sqr( 0.5 ) * M_PI );
    auto& mat           = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat.m_albedo        = color3( 1.f );
    mat.m_emission      = color3( lightIntensity );

    Transform modelMatrix;
    modelMatrix.scale( 0.5, 0.5, 0.5 );
    modelMatrix.translate( vec3( 10, 10, 4 ) );
    // addObject(scene, initSphere(mat, modelMatrix));

    auto obj   = initSphere( mat.m_UID, modelMatrix );
    auto light = new ShapeLight( point3( 10, 10, 4 ), color3( lightIntensity ), obj );
    // addLight(scene, light);

    auto& mat1           = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    auto lightIntensity1 = 100;
    mat1.m_albedo        = color3( 1.f );
    mat1.m_emission      = color3( lightIntensity1 );
    Transform t0;
    t0.scale( 0.1, 0.1, 0.1 );
    t0.translate( vec3( -1.25, 0., 0 ) );
    addObject( scene, initSphere( mat1.m_UID, t0 ) );

    auto obj0   = initSphere( mat1.m_UID, t0 );
    auto light0 = new ShapeLight( point3( -1.25, 0, 0 ), color3( lightIntensity1 ), obj0 );
    addLight( scene, light0 );

    auto& mat2           = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    auto lightIntensity2 = 901.803;
    mat2.m_albedo        = color3( 1.f );
    mat2.m_emission      = color3( lightIntensity2 );
    Transform t1;
    t1.scale( 0.03333, 0.03333, 0.03333 );
    t1.translate( vec3( -3.75, 0., 0 ) );
    addObject( scene, initSphere( mat2.m_UID, t1 ) );

    auto obj1   = initSphere( mat2.m_UID, t1 );
    auto light1 = new ShapeLight( point3( -3.75, 0, 0 ), color3( lightIntensity2 ), obj1 );
    addLight( scene, light1 );

    auto& mat3           = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    auto lightIntensity3 = 11.1111;
    mat3.m_albedo        = color3( 1.f );
    mat3.m_emission      = color3( lightIntensity3 );
    Transform t2;
    t2.scale( 0.3, 0.3, 0.3 );
    t2.translate( vec3( 1.25, 0., 0 ) );
    addObject( scene, initSphere( mat3.m_UID, t2 ) );

    auto obj2   = initSphere( mat3.m_UID, t2 );
    auto light2 = new ShapeLight( point3( 1.25, 0, 0 ), color3( lightIntensity3 ), obj2 );
    addLight( scene, light2 );

    auto& mat4           = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    auto lightIntensity4 = 1.23457 * 4.f * M_PI / ( 4. * M_PI * sqr( 0.9 ) * M_PI );
    lightIntensity4      = 4;
    lightIntensity4      = 1.23457;
    mat4.m_albedo        = color3( 1.f );
    mat4.m_emission      = color3( lightIntensity4 );
    Transform t3;
    t3.scale( 0.9, 0.9, 0.9 );
    t3.translate( vec3( 3.75, 0., 0 ) );
    addObject( scene, initSphere( mat4.m_UID, t3 ) );

    auto obj3   = initSphere( mat4.m_UID, t3 );
    auto light4 = new ShapeLight( point3( 3.75, 0, 0 ), color3( lightIntensity4 ), obj3 );
    addLight( scene, light4 );

    auto& mat5       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, DIFFUSE );
    mat5.m_roughness = 0.005;
    mat5.m_metalness = 0.0;
    mat5.m_albedo    = color3( 0.07, 0.09, 0.13 );
    Transform t4;
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/veach_mi/plate1.obj", scene, mat5, t4 );

    auto& mat6       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, DIFFUSE );
    mat6.m_roughness = 0.02;
    mat6.m_metalness = 0.0;
    mat6.m_albedo    = color3( 0.07, 0.09, 0.13 );
    Transform t5;
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/veach_mi/plate2.obj", scene, mat6, t5 );

    auto& mat7       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE, DIFFUSE );
    mat7.m_roughness = 0.05;
    mat7.m_metalness = 0.0;
    mat7.m_albedo    = color3( 0.07, 0.09, 0.13 );
    Transform t6;
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/veach_mi/plate3.obj", scene, mat7, t6 );

    auto& mat8       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat8.m_roughness = 0.1;
    mat8.m_metalness = 0.0;
    mat8.m_albedo    = color3( 0.07, 0.09, 0.13 );
    Transform t7;
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/veach_mi/plate4.obj", scene, mat8, t7 );

    auto& mat9       = scene->CreateMaterial( MaterialModel::COOK_TORRANCE );
    mat9.m_albedo    = color3( 0.4, 0.4, 0.4 );
    mat9.m_metalness = 0.0;
    mat9.m_roughness = 1.f;
    Transform t8;
    addObjectsFromFile(
        "C:/Users/marco/Documents/GitRepos/Raytracer/assets/veach_mi/floor.obj", scene, mat9, t8 );

    return scene;
}

Scene* parseScene( int sceneId ) {
    Scene* scene = NULL;
    switch ( sceneId ) {
    case 9:
        scene = initScene9();
        break;
    case 10:
        scene = initScene10();
        break;
    case 11:
        scene = initScene11();
        break;
    case 12:
        scene = initScene12();
        break;
    case 13:
        scene = initScene13();
        break;
    case 14:
        scene = initScene14();
        break;
    case 15:
        scene = initScene15();
        break;
    case 16:
        scene = initScene16();
        break;
    case 17:
        scene = initScene17();
        break;
    case 18:
        scene = initScene18();
        break;
    case 19:
        scene = initScene19();
        break;
    case 20:
        scene = initScene20();
        break;
    case 21:
        scene = initScene21();
        break;

    default:
        scene = initScene19();
        break;
    }
    return scene;
}
