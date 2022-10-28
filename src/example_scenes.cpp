#include "example_scenes.h"

#include "defines.h"
#include "image.h"
#include "ray.h"
#include "raytracer.h"
#include "scene.h"
#include "textures.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Light.h"
#include "Object.h"
#include "materials/Blinn.h"
#include "shapes/plane.h"
#include "Camera.h"

#include <glm/gtc/matrix_transform.hpp>

//CookTorrance* mat_lib[] = {
//    /* bunny glass */
//    {1.05, 2.2, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, TRANSPARENT },
//
//    /* specular black phenolic */
//    {1.072, 0.0588, {1.0, 0.824, 0.945}, {0.002, 0.002, 0.003}},
//
//    /* specular blue phenolic */
//    {1.1051, 0.0568, {0.005, 0.013, 0.032}, {1.0, 0.748, 0.718}},
//
//    /* specular green phenolic */
//    {1.1051, 0.0567, {0.006, 0.026, 0.022}, {1.0, 0.739, 0.721}},
//
//    /* specular white phenolic */
//    {1.1022, 0.0579, {0.286, 0.235, 0.128}, {1.0, 0.766, 0.762}},
//
//    /* marron plastic */
//    {1.0893, 0.0604, {0.202, 0.035, 0.033}, {1.0, 0.857, 0.866}},
//
//    /* purple paint */
//    {1.1382, 0.0886, {0.301, 0.034, 0.039}, {1.0, 0.992, 0.98}},
//
//    /* red specular plastic */
//    {1.0771, 0.0589, {0.26, 0.036, 0.014}, {1.0, 0.852, 1.172}},
//
//    /* green acrylic */
//    {1.1481, 0.0625, {0.016, 0.073, 0.04}, {1.0, 1.056, 1.146}},
//
//    /* blue acrylic */
//    {1.1153, 0.068, {0.012, 0.036, 0.106}, {1.0, 0.965, 1.07}},
//
//    {1.51, 2.2, {1.0, 1.0, 1.0}, {1.0, .0, .0}, TRANSPARENT },
//};

//Material mat_lib[] = {
//    /* nickel */
//    {2.4449, 0.0681, {1.0, 0.882, 0.786}, {0.014, 0.012, 0.012}},
//
//    /* specular black phenolic */
//    {1.072, 0.0588, {1.0, 0.824, 0.945}, {0.002, 0.002, 0.003}},
//
//    /* specular blue phenolic */
//    {1.1051, 0.0568, {0.005, 0.013, 0.032}, {1.0, 0.748, 0.718}},
//
//    /* specular green phenolic */
//    {1.1051, 0.0567, {0.006, 0.026, 0.022}, {1.0, 0.739, 0.721}},
//
//    /* specular white phenolic */
//    {1.1022, 0.0579, {0.286, 0.235, 0.128}, {1.0, 0.766, 0.762}},
//
//    /* marron plastic */
//    {1.0893, 0.0604, {0.202, 0.035, 0.033}, {1.0, 0.857, 0.866}},
//
//    /* purple paint */
//    {1.1382, 0.0886, {0.301, 0.034, 0.039}, {1.0, 0.992, 0.98}},
//
//    /* red specular plastic */
//    {1.0771, 0.0589, {0.26, 0.036, 0.014}, {1.0, 0.852, 1.172}},
//
//    /* green acrylic */
//    {1.1481, 0.0625, {0.016, 0.073, 0.04}, {1.0, 1.056, 1.146}},
//
//    /* blue acrylic */
//    {1.1153, 0.068, {0.012, 0.036, 0.106}, {1.0, 0.965, 1.07}}};

Scene *initScene0() {
  Scene *scene = initScene();
  auto from = point3(1, 0.6, 1);
  auto at = vec3(0, 0.3, 0);
  setSimpleCamera(scene, from, at, vec3(0, 1, 0), 60,
            (float)WIDTH, (float)HEIGHT, 0., glm::length(at - from));
   image_texture* sky = new image_texture("../assets/clouds.png");
  scene->m_skyTexture = sky;
  Transform texSky;
  texSky.scale(1,.4,1);
  texSky.translate(vec3(0,-0.1,0));
  scene->m_skyTexture->m_transform = texSky;

 setSkyColor(scene, color3(0.52f, 0.8f, 0.9f));
  //auto mat = std::make_shared<CookTorrance>(false);
  //mat->m_IOR = 4.0;
  //mat->m_roughness = 0.05;
  //mat->m_specularColor = color3(1.0f);
  //mat->m_diffuseColor = color3(1.0f);
  //Transform t0;
  //t0.translate(vec3(0,0.25,0));
  //t0.scale(0.25,0.25,0.25);
  //addObject(scene, initSphere(mat, t0));

  //auto mat1 = std::make_shared<CookTorrance>(false);
  //mat1->m_IOR = 4.0;
  //mat1->m_roughness = 0.05;
  //mat1->m_specularColor = color3(1.f, 0.f, 0.f);
  //mat1->m_diffuseColor = color3(1.f, 0.f, 0.f);
  //Transform t1;
  //t1.translate(vec3(0.5,0.25,0));
  //t1.scale(0.25,0.25,0.25);
  //addObject(scene, initSphere(mat1, t1));

  auto mat2 = std::make_shared<CookTorrance>(false);
  mat2->m_IOR = 1.3;
  mat2->m_roughness = 0.5;
  mat2->m_transparent = true;
  mat2->m_specularColor = color3(1.f, 1.f, 1.f);
  mat2->m_diffuseColor = color3(1.f, 1.f, 1.f);
  Transform t2;
  t2.translate(vec3(0.3,1.25,0.3));
  t2.scale(0.25,0.25,0.25);
  addObject(scene, initSphere(mat2, t2));

  auto mat3 = std::make_shared<CookTorrance>(false);
  mat3->m_IOR = 4.0;
  mat3->m_roughness = 1.05;
  mat3->m_specularColor = color3(0.f, 0.f, 1.f);
  mat3->m_diffuseColor = color3(0.f, 0.f, 1.f);
  Transform t3;
  t3.translate(vec3(-6,1.0,-6.0));
  t3.scale(0.08,0.08,0.08);
  addObject(scene, initSphere(mat3, t3));

  auto mat4 = std::make_shared<CookTorrance>(false);
  mat4->m_IOR = 4.0;
  mat4->m_roughness = 0.05;
  mat4->m_specularColor = color3(0.6f);
  mat4->m_diffuseColor = color3(0.6f);

  Transform t4;
  t4.rotate(vec3(-1,0,0), 90);
  auto ret = new Plane(mat4, t4);
  ret->geom.type = PLANE;
  addObject(scene, ret);

  //addLight(scene, initPointLight(point3(1, 5, 1), color3(1, 1, 1)));
  addLight(scene, initDirectLight(vec3(0, -1, 0), color3(0.6)));
  //addLight(scene, initDirectLight(vec3(1, 0.3, -1), color3(0.4)));
  //addLight(scene, initPointLight(point3(4, 10, -2), color3(1, 1, 1)));
  //addLight(scene, initPointLight(point3(4, 5, 0), color3(3, 3, 3)));
  
  //auto light = new AreaLight(vec3(13,10,13), vec3(-6,0,0), 10, vec3(0,0,-6), 10, color3(1,1,1));
  //addLight(scene, light);

  //light = new AreaLight(vec3(1,7,-5), vec3(6,3,0), 10, vec3(0,3,6), 10, color3(1,1,1));
  //addLight(scene, light);

  return scene;
}

Scene *initScene1() {

  Scene *scene = initScene();
  //setCamera(scene, point3(-9, 0, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
  //          (float)WIDTH / (float)HEIGHT);
  //setSkyColor(scene, color3(0.2, 0.2, 0.7));

  //Material mat->
  //mat->m_IOR = 1.12;
  //mat->m_roughness = 2.0;
  //mat->m_specularColor = color3(0.4f);
  //mat->m_diffuseColor = color3(0.6f);

  //for (int i = 0; i < 10; ++i) {
  //  mat->m_diffuseColor = color3(0.301, 0.034, 0.039);
  //  mat->m_specularColor = color3(1.0, 0.992, 0.98);
  //  mat->m_IOR = 1.1382;
  //  mat->m_roughness = 0.0886;
  //  mat->m_roughness = ((float)10 - i) / (10 * 9.f);
  //  addObject(scene, initSphere(point3(0, 0, -1.5 + i / 9.f * 3.f), .15, mat->lib[i]));
  //}
  //for (int i = 0; i < 10; ++i) {
  //  mat->m_diffuseColor = color3(0.012, 0.036, 0.106);
  //  mat->m_specularColor = color3(1.0, 0.965, 1.07);
  //  mat->m_IOR = 1.1153;
  //  mat->m_roughness = 0.068;
  //  mat->m_roughness = ((float)i + 1) / 10.f;
  //  addObject(scene, initSphere(point3(0, 1, -1.5 + i / 9.f * 3.f), .15, mat->lib[i]));
  //}
  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(1.0, 0.882, 0.786);
  //mat->m_IOR = 2.4449;
  //mat->m_roughness = 0.0681;
  //addObject(scene, initSphere(point3(-3.f, 1.f, 0.f), 2., mat->lib[0]));
  //
  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(1.0, 0.882, 0.786);
  //mat->m_IOR = 2.4449;
  //mat->m_roughness = 0.0681;
  //addObject(scene, initSphere(point3(-12.f, 1.f, 0.f), 0.8, mat->lib[1]));

  //mat->m_diffuseColor = color3(0.016, 0.073, 0.04);
  //mat->m_specularColor = color3(1.0, 1.056, 1.146);
  //mat->m_IOR = 1.1481;
  //mat->m_roughness = 0.0625;
  //addObject(scene, initPlane(vec3(0, 1, 0), +1, mat->);

  //addLight(scene, initPointLight(point3(10, 10.7, 1), .5f * color3(3, 3, 3)));
  //addLight(scene, initPointLight(point3(8, 20, 3), .5f * color3(4, 4, 4)));
  //addLight(scene, initPointLight(point3(4, 30, -1), .5f * color3(5, 5, 5)));
  ////addLight(scene, initPointLight(point3(10, 10, 10), color3(1, 1, 1)));
  ////addLight(scene, initPointLight(point3(4, 10, -2), color3(1, 1, 1)));
  return scene;
}

Scene *initScene2() {
  Scene *scene = initScene();
  //setCamera(scene, point3(0.5, 3, 1), vec3(0, 0, 0.6), vec3(0, 0, 1), 60,
  //          (float)WIDTH / (float)HEIGHT);
  //setSkyColor(scene, color3(0.2, 0.2, 0.7));
  //Material mat->
  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(1.0, 0.882, 0.786);
  //mat->m_IOR = 2.4449;
  //mat->m_roughness = 0.0681;

  //mat->m_diffuseColor = color3(0.05, 0.05, 0.05);
  //mat->m_specularColor = color3(0.95);
  //mat->m_IOR = 1.1022;
  //mat->m_roughness = 0.0579;

  //addObject(scene, initPlane(vec3(0, 0, 1), 0, mat->);

  //mat->m_diffuseColor = color3(0.005, 0.013, 0.032);
  //mat->m_specularColor = color3(1.0, 0.748, 0.718);
  //for (int i = 0; i < 4; ++i) {
  //  mat->m_IOR = 1.1051 + (-0.1 + float(i) / 3.f * 0.4);
  //  for (int j = 0; j < 10; ++j) {
  //    mat->m_roughness = 0.0568 + (-0.1 + float(j) / 9.f * 0.3);
  //    addObject(scene, initSphere(point3(-1.5 + float(j) / 9.f * 3.f, 0,
  //                                       0.4 + float(i) * 0.4f),
  //                                .15, mat->);
  //  }
  //}

  //addLight(scene, initPointLight(point3(-20, 5, 10), color3(30, 30, 30)));
  //addLight(scene, initPointLight(point3(10, 10, 10), color3(30, 30, 30)));
  //addLight(scene, initPointLight(point3(50, -100, 10), color3(1, 0.7, 2)));
  return scene;
}

Scene *initScene3() {
  Scene *scene = initScene();
  //setCamera(scene, point3(4.5, .8, 4.5), vec3(0, 0.3, 0), vec3(0, 1, 0), 46.5,
  //          (float)WIDTH / (float)HEIGHT);
  //setSkyColor(scene, color3(0.2, 0.2, 0.7));
  //Material mat->
  //mat->m_diffuseColor = color3(0.301, 0.034, 0.039);
  //mat->m_specularColor = color3(1.0, 0.992, 0.98);
  //mat->m_IOR = 1.1382;
  //mat->m_roughness = 0.0886;

  //addLight(scene, initPointLight(point3(0, 1.7, 1), .5f * color3(3, 3, 3)));
  //addLight(scene, initPointLight(point3(3, 2, 3), .5f * color3(4, 4, 4)));
  //addLight(scene, initPointLight(point3(4, 3, -1), .5f * color3(5, 5, 5)));

  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(0.7, 0.882, 0.786);
  //mat->m_IOR = 6;
  //mat->m_roughness = 0.0181;
  //auto transform = glm::translate(glm::mat->(1.f), vec3(0,0.1,0)) * glm::scale(glm::mat->(1.f), vec3(0.3,0.3,0.3));
  //addObject(scene, initSphere(mat-> transform));

  //mat->m_diffuseColor = color3(0.26, 0.036, 0.014);
  //mat->m_specularColor = color3(1.0, 0.852, 1.172);
  //mat->m_IOR = 1.3771;
  //mat->m_roughness = 0.01589;
  //transform = glm::translate(glm::mat->(1.f), vec3(1,-.05,0)) * glm::scale(glm::mat->(1.f), vec3(0.15,0.15,0.15));
  //addObject(scene, initSphere(mat-> transform));

  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(0.7, 0.882, 0.786);
  //mat->m_IOR = 3;
  //mat->m_roughness = 0.00181;
  //transform = glm::translate(glm::mat->(1.f), vec3(3,0.05,2)) * glm::scale(glm::mat->(1.f), vec3(0.25,0.25,0.25));
  //addObject(scene, initSphere(mat-> transform));

  //mat->m_diffuseColor = color3(0.46, 0.136, 0.114);
  //mat->m_specularColor = color3(0.8, 0.852, 0.8172);
  //mat->m_IOR = 1.5771;
  //mat->m_roughness = 0.01589;
  //transform = glm::translate(glm::mat->(1.f), vec3(1.3,0.,2.6)) * glm::scale(glm::mat->(1.f), vec3(0.215,0.215,0.215));
  //addObject(scene, initSphere(mat-> transform));

  //mat->m_diffuseColor = color3(0.06, 0.26, 0.22);
  //mat->m_specularColor = color3(0.70, 0.739, 0.721);
  //mat->m_IOR = 1.3051;
  //mat->m_roughness = 0.567;
  //transform = glm::translate(glm::mat->(1.f), vec3(1.9,0.05,2.2)) * glm::scale(glm::mat->(1.f), vec3(0.25,0.25,0.25));
  //addObject(scene, initSphere(mat-> transform));

  //mat->m_diffuseColor = color3(0.012, 0.036, 0.406);
  //mat->m_specularColor = color3(1.0, 0.965, 1.07);
  //mat->m_IOR = 1.1153;
  //mat->m_roughness = 0.068;
  //mat->m_roughness = 0.18;
  //transform = glm::translate(glm::mat->(1.f), vec3(0,0,1)) * glm::scale(glm::mat->(1.f), vec3(0.2,0.2,0.2));
  //addObject(scene, initSphere(mat-> transform));

  //mat->m_diffuseColor = color3(.2, 0.4, .3);
  //mat->m_specularColor = color3(.2, 0.2, .2);
  //mat->m_IOR = 1.382;
  //mat->m_roughness = 0.05886;
  //addObject(scene, initPlane(vec3(0, 1, 0), 0.2, mat->);

  //mat->m_diffuseColor = color3(.5, 0.09, .07);
  //mat->m_specularColor = color3(.2, .2, .1);
  //mat->m_IOR = 1.8382;
  //mat->m_roughness = 0.886;
  //addObject(scene, initPlane(vec3(1, 0.0, -1.0), 2, mat->);

  //mat->m_diffuseColor = color3(0.1, 0.3, .05);
  //mat->m_specularColor = color3(.5, .5, .5);
  //mat->m_IOR = 1.9382;
  //mat->m_roughness = 0.0886;
  //addObject(scene, initPlane(vec3(0.3, -0.2, 1), 4, mat->);
  return scene;
}

Scene *initScene4() {
  Scene *scene = initScene();
  //setCamera(scene, point3(6, 4, 6), vec3(0, 1, 0), vec3(0, 1, 0), 90,
  //          (float)WIDTH / (float)HEIGHT);
  //setSkyColor(scene, color3(0.2, 0.2, 0.7));
  //Material mat->
  //mat->m_diffuseColor = color3(0.301, 0.034, 0.039);
  //mat->m_specularColor = color3(1.0, 0.992, 0.98);
  //mat->m_IOR = 1.1382;
  //mat->m_roughness = 0.0886;

  //addLight(scene, initPointLight(point3(10, 10.7, 1), .5f * color3(3, 3, 3)));
  //addLight(scene, initPointLight(point3(8, 20, 3), .5f * color3(4, 4, 4)));
  //addLight(scene, initPointLight(point3(4, 30, -1), .5f * color3(5, 5, 5)));

  //mat->m_diffuseColor = color3(.2, 0.4, .3);
  //mat->m_specularColor = color3(.2, 0.2, .2);
  //mat->m_IOR = 2.382;
  //mat->m_roughness = 0.005886;
  //addObject(scene, initPlane(vec3(0, 1, 0), 0.2, mat->);

  //mat->m_diffuseColor = color3(.5, 0.09, .07);
  //mat->m_specularColor = color3(.2, .2, .1);
  //mat->m_IOR = 2.8382;
  //mat->m_roughness = 0.00886;
  //addObject(scene, initPlane(vec3(1, 0.0, 0.0), 2, mat->);

  //mat->m_diffuseColor = color3(0.1, 0.3, .05);
  //mat->m_specularColor = color3(.5, .5, .5);
  //mat->m_IOR = 2.9382;
  //mat->m_roughness = 0.00886;
  //addObject(scene, initPlane(vec3(0, 0, 1), 4, mat->);

  //for (int i = 0; i < 600; i++) {
  //  addObject(scene,
  //            initSphere(point3(1 + rand() % 650 / 100.0, rand() % 650 / 100.0,
  //                              1 + rand() % 650 / 100.0),
  //                       .05 + rand() % 200 / 1000.0, mat->lib[rand() % 10]));
  //}
  return scene;
}
Scene *initScene5() {

  Scene *scene = initScene();
  //setCamera(scene, point3(3, 1, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
  //          (float)WIDTH / (float)HEIGHT);
  //setSkyColor(scene, color3(0.1f, 0.3f, 0.5f));
  //Material mat->
  //mat->m_diffuseColor = color3(0.301, 0.034, 0.039);
  //mat->m_specularColor = color3(1.0, 0.992, 0.98);
  //mat->m_IOR = 1.1382;
  //mat->m_roughness = 0.0886;
  ////addLight(scene, initPointLight(point3(0, 1.7, 1), .5f * color3(3, 3, 3)));
  ////addLight(scene, initPointLight(point3(3, 2, 3), .5f * color3(4, 4, 4)));
  ////addLight(scene, initPointLight(point3(4, 3, -1), .5f * color3(5, 5, 5)));
  ////addLight(scene, initPointLight(point3(1, 0, 1), .5f * color3(3, 3, 3)));

  //vec3 v0 = vec3(1,0,0);
  //vec3 v1 = vec3(0,1,0);
  //vec3 v2 = vec3(0,0,-1);

  //vec3 v1v0 = v1-v0;
  //vec3 v2v0 = v2-v0;
  //vec3 n = normalize(cross(v2v0, v1v0));

  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(0.7, 0.882, 0.786);
  //mat->m_IOR = 6;
  //mat->m_roughness = 0.0181;
  //vec2 textures[3];
  //addObject(scene, initTriangle(v0, v1, v2, n,textures,mat->lib[0]));

  //v1 = vec3(1,0,0);
  //v0 = vec3(0,1,0);
  //v2 = vec3(0,0,1);

  //v1v0 = v1-v0;
  //v2v0 = v2-v0;
  //n = normalize(cross(v2v0, v1v0));
  //mat->m_diffuseColor = color3(0.5, 0,0);
  //mat->m_specularColor = color3(0.5, 0, 0);
  //mat->m_IOR = 1.01;
  //mat->m_roughness = 2.2;
  //mat->mtype = DIFFUSE;
  //mat->m_texture = new image_texture("../assets/container2.png");
  //textures[0] = vec2(0,0);
  //textures[1] = vec2(1,0);
  //textures[2] = vec2(0,1);
  //addObject(scene, initTriangle(v0, v1, v2, n,textures,mat->);
  //v0 = vec3(1,-1, 1);
  //textures[0] = vec2(1,1);
  //textures[1] = vec2(1,0);
  //textures[2] = vec2(0,1);
  //addObject(scene, initTriangle(v0, v1, v2, n,textures,mat->);
  //
  //mat->m_texture = nullptr;
  //mat->m_diffuseColor = color3(0.016, 0.073, 0.04);
  //mat->m_specularColor = color3(1.0, 1.056, 1.146);
  //mat->m_IOR = 1.1481;
  //mat->m_roughness = 0.0625;
  //addObject(scene, initPlane(vec3(0, 1, 0), +1, mat->);

  //addObject(scene, initSphere(point3(-3, 1, 0), 2, mat->);

  //addLight(scene, initPointLight(point3(4, 2, 0), color3(1,1,1)));
  //addLight(scene, initPointLight(point3(10, 100, 100), color3(1, 1, 1)));
  //addLight(scene, initPointLight(point3(10, 10, -16), color3(1, 1, 1)));
  return scene;
}

void addObjectsFromFile(const char *filename, Scene *scene, std::shared_ptr<Material> default_mat, Transform& transform){
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  auto mat = std::make_shared<CookTorrance>(false);
  
  readObjToTriangleMesh(filename, attrib, shapes, materials);

  for(size_t s = 0; s < shapes.size(); s++){
    size_t index_offset = 0;
    for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size();f++){
      int fv = shapes[s].mesh.num_face_vertices[f];

      std::vector<point3> vector;
      std::vector<point3> normals;
      std::vector<vec2> texture;
      vec2 textures[fv];

      // Loop over vertices in the face.
      for (int v = 0; v < fv; v++) {
        // access to vertex
        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
        tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
        tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
        tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
        tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
        tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
        tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
        if(!attrib.texcoords.empty()){
          tinyobj::real_t tx = attrib.texcoords[2*idx.texcoord_index+0];
          tinyobj::real_t ty = attrib.texcoords[2*idx.texcoord_index+1];
          textures[v] = vec2(tx, ty);
        }
        // Optional: vertex colors
        // tinyobj::real_t red = attrib.colors[3*idx.vertex_index+0];
        // tinyobj::real_t green = attrib.colors[3*idx.vertex_index+1];
        // tinyobj::real_t blue = attrib.colors[3*idx.vertex_index+2];
        vec3 v0 = point3(vx, vy, vz);
        vec3 n = point3(nx, ny, nz);
        vector.push_back(v0);
        normals.push_back(n);
      }

      index_offset += fv;

      vec3 v0 = vector[0];
      vec3 v1 = vector[1];
      vec3 v2 = vector[2];
      //vec3 n = normalize(cross(v2v0, v1v0));

      vec3 n = normalize(normals[0] + normals[1] + normals[2]);

      if(!materials.empty()){
        int matIndex = shapes[s].mesh.material_ids[f];
        mat->m_IOR = materials[matIndex].ior;
        mat->m_roughness = materials[matIndex].ior;
        mat->m_diffuseColor.r = materials[matIndex].diffuse[0];
        mat->m_diffuseColor.g = materials[matIndex].diffuse[1];
        mat->m_diffuseColor.b = materials[matIndex].diffuse[2];
        //mat.m_texture = new image_texture(materials[matIndex].diffuse_texname.c_str());
        if(materials[matIndex].specular[0] == 0 && materials[matIndex].specular[1] == 0 && materials[matIndex].specular[1] == 0){
          mat->m_specularColor.r = materials[matIndex].diffuse[0];
          mat->m_specularColor.g = materials[matIndex].diffuse[1];
          mat->m_specularColor.b = materials[matIndex].diffuse[2];
        }
        else{
          mat->m_specularColor.r = materials[matIndex].specular[0];
          mat->m_specularColor.g = materials[matIndex].specular[1];
          mat->m_specularColor.b = materials[matIndex].specular[2];
        }
        addObject(scene, initSmoothTriangle(v0, v1, v2, n, textures, normals[0], normals[1], normals[2], transform, mat));
      }
      else{
        addObject(scene, initSmoothTriangle(v0, v1, v2, n, textures, normals[0], normals[1], normals[2], transform, default_mat));
      }
    }
  }
}

Scene *initScene6() {

  Scene *scene = initScene();
  //setCamera(scene, point3(-4.5, 2.5, 3.5), vec3(0, 0, 0), vec3(0, 0.5, 0), 60,
  //          (float)WIDTH / (float)HEIGHT);
  //setSkyColor(scene, color3(0.459f, 0.f, 0.878f));
  //
  ////addLight(scene, initPointLight(point3(100, 1, 1), color3(50, 50, 50)));
  ////addLight(scene, initPointLight(point3(1, 0.1, 0.5), color3(5, 5, 5)));
  ////addLight(scene, initPointLight(point3(0, 0, 0), color3(5, 5, 5)));
  ////addLight(scene, initPointLight(point3(-6, 3, 7), color3(10, 10, 10)));

  ////addLight(scene, initPointLight(point3(6, 3, -7), color3(1, 1, 1)));

  //addLight(scene, initPointLight(point3(-6, 4, 0), color3(1, 1, 1)));
  //addLight(scene, initPointLight(point3(4, 10, 10), color3(1, 1, 1)));

  //addLight(scene, initPointLight(point3(-5, 4, 5), color3(1, 1, 1)));
  //addLight(scene, initPointLight(point3(-5, 4, 4), color3(1, 1, 1)));

  ////addLight(scene, initPointLight(point3(20, 1, -20), color3(1, 1, 1)));

  //Material mat->
  //mat->m_IOR = 1.3f;
  //mat->m_diffuseColor = color3(.001f);
  //mat->m_specularColor = color3(1.0f);
  //mat->m_roughness = 1.f;

  //addObjectsFromFile("../assets/new_bunny.obj", scene, mat->lib[0]);

  //mat->m_diffuseColor = color3(0.7, 0, 0);
  //mat->m_specularColor = color3(0.95, 0, 0);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.5;
  //mat->mtype = DIFFUSE;
  //addObject(scene, initSphere(point3(4, 0., -2), 0.4, mat->);

  //mat->m_diffuseColor = color3(0., 0, 0.7);
  //mat->m_specularColor = color3(0., 0, 0.95);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.5;
  //mat->mtype = DIFFUSE;
  //addObject(scene, initSphere(point3(-0.5, 0., 0), 0.2, mat->);

  //mat->m_diffuseColor = color3(0.7, 0.5, 0);
  //mat->m_specularColor = color3(0.95, 0.6, 0);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.5;
  //mat->mtype = DIFFUSE;
  //addObject(scene, initSphere(point3(-0.9, 0., 3), 0.2, mat->);

  //addObject(scene, initSphere(point3(4, 0., 3), 0.2, mat->lib[2]));

  //mat->m_diffuseColor = color3(0.5);
  //mat->m_specularColor = color3(0.5);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.0681;
  //mat->mtype = DIFFUSE;
  //mat->m_texture = new checker_texture(color3(1.0, 1.0, 1.0), color3(0, 0, 0));

  //mat->m_diffuseColor = color3(0.5);
  //mat->m_specularColor = color3(0.5);
  //mat->m_IOR = 1.1;
  ////mat->m_roughness = 0.0681;
  //mat->m_roughness = 1.2;
  //mat->mtype = DIFFUSE;
  ////mat->m_texture = new checker_texture(color3(0.2, 0.3, 0.1), color3(0.9, 0.9, 0.9));
  //mat->m_texture = new image_texture("../assets/chessboardtexture.png");
  //addObject(scene, initPlane(vec3(0, 1, 0), 0.4, mat->);
  //mat->m_texture = nullptr;

  //mat->m_diffuseColor = color3(0, 0.8, 0);
  //mat->m_specularColor = color3(0.6, 0.6, 0);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.0681;
  ////addObject(scene, initPlane(vec3(0, 0, 1), +10, mat->);

  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(0.7, 0.882, 0.786);
  //mat->m_IOR = 6;
  //mat->m_roughness = 0.0181;
  //addObject(scene, initSphere(point3(0, 0., -4), .3, mat->);
  //mat->m_texture = nullptr;

  //mat->m_diffuseColor = color3(0.26, 0.036, 0.014);
  //mat->m_specularColor = color3(1.0, 0.852, 1.172);
  //mat->m_IOR = 1.3771;
  //mat->m_roughness = 0.01589;
  //addObject(scene, initSphere(point3(-1, -.05, 0), .15, mat->);

  //mat->m_diffuseColor = color3(0.014, 0.012, 0.012);
  //mat->m_specularColor = color3(0.7, 0.882, 0.786);
  //mat->m_IOR = 3;
  //mat->m_roughness = 0.00181;
  //addObject(scene, initSphere(point3(1, 0.05, 2), .25, mat->);
  //mat->m_texture = nullptr;

  //mat->m_diffuseColor = color3(0.46, 0.136, 0.114);
  //mat->m_specularColor = color3(0.8, 0.852, 0.8172);
  //mat->m_IOR = 1.5771;
  //mat->m_roughness = 0.01589;
  //addObject(scene, initSphere(point3(1.6, 0., 2.6), 0.215, mat->);

  //mat->m_diffuseColor = color3(0.06, 0.26, 0.22);
  //mat->m_specularColor = color3(0.70, 0.739, 0.721);
  //mat->m_IOR = 1.3051;
  //mat->m_roughness = 0.567;
  //mat->m_texture = new image_texture("../assets/earthmap.png");
  //addObject(scene, initSphere(point3(1.9, 0.05, 2.2), .25, mat->);

  //mat->m_diffuseColor = color3(0.012, 0.036, 0.406);
  //mat->m_specularColor = color3(1.0, 0.965, 1.07);
  //mat->m_IOR = 1.1153;
  //mat->m_roughness = 0.068;
  //mat->m_roughness = 0.18;
  //addObject(scene, initSphere(point3(-2, 0, 1.5), .20, mat->);
  
  return scene;
}

Scene *initScene7() {

  Scene *scene = initScene();
  //setCamera(scene, point3(8, 5, 7), vec3(0, 0, 0), vec3(0, 1, 0), 60,
  //          float(WIDTH) / float(HEIGHT));
  //setSkyColor(scene, color3(0.2, 0.8, 0.7)); 

  //Material mat->
  //mat->m_diffuseColor = color3(0.5);
  //mat->m_specularColor = color3(0.5);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.0681;
  //mat->mtype = DIFFUSE;
  //addObjectsFromFile("../assets/werewolf.obj", scene, mat->;
  
  //mat->m_diffuseColor = color3(.4, 0.8, .4);
  //mat->m_specularColor = color3(.4, 0.6, .2);
  //mat->m_IOR = 1.382;
  //mat->m_roughness = 0.05886;
  //addObject(scene, initPlane(vec3(100, 1, 0), +100, mat->);
  
  //Material mat->
  //mat->m_diffuseColor = color3(0.5);
  //mat->m_specularColor = color3(0.5);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.0681;
  //mat->mtype = DIFFUSE;
  //mat->m_texture = new checker_texture(color3(0.2, 0.3, 0.1), color3(0.9, 0.9, 0.9));
  //addObject(scene, initPlane(vec3(0, 1, 0), 1, mat->);

  //mat->m_texture = nullptr;
  //mat->m_diffuseColor = color3(1.0f, 0, 0);
  //mat->m_specularColor = color3(0.5f, 0.f, 0.f);
  //mat->m_roughness = 0.5f;
  //addObject(scene, initSphere(point3(-2, 0, -1), .25, mat->);

  //mat->m_diffuseColor = color3(0, 1.0f, 0);
  //mat->m_specularColor = color3(0.f, 1.f, 0.f);
  //mat->m_roughness = 0.5f;
  //addObject(scene, initSphere(point3(0, 0, -1), .25, mat->);

  //mat->m_diffuseColor = color3(0.5, 0 ,1.0f);
  //mat->m_specularColor = color3(0.5f, 0.f, 1.f);
  //mat->m_roughness = 0.5f;
  //addObject(scene, initSphere(point3(-1, 0, 0), .25, mat->);

  //mat->m_diffuseColor = color3(0.5, 1.0f, 0.5);
  //mat->m_specularColor = color3(0.5f, 1.f, 0.5f);
  //mat->m_roughness = 0.5f;
  //addObject(scene, initSphere(point3(-2, 0, 0), .25, mat->);

  addLight(scene, initPointLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initPointLight(point3(4, 10, -2), color3(1, 1, 1)));
  addLight(scene, initPointLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initPointLight(point3(0, 40, 0), color3(1, 1, 1)));
  return scene;
}

Scene *initScene8() {
  Scene *scene = initScene();
//  setCamera(scene, point3(8, 3, -4), vec3(0, 1, -1), vec3(0, 1, 0), 88,
//            (float)WIDTH / (float)HEIGHT, 0.05, glm::length(point3(8, 3, -4) - vec3(3, 2.5, -2)) / 2.f);
//  setSkyColor(scene, color3(0.4, 0.9, 0.9));
//
//  //addLight(scene, initPointLight(point3(52, 10, 52), color3(1, 1, 1)));
//  //addLight(scene, initPointLight(point3(52, 10, 16), color3(1, 1, 1)));
//  //addLight(scene, initPointLight(point3(16, 10, 52), color3(1, 1, 1)));
//
//  //addLight(scene, initPointLight(point3(0, 50, 0), color3(1, 1, 1)));
//  
//  glm::mat-> glassSphereT = glm::translate(glm::mat->(1.f), vec3(3,2.5,-2)) *glm::scale(glm::mat->(1.f), vec3(2,2,2));
//  addObject(scene, initSphere(mat->lib[10], glassSphereT));
//
//  //addLight(scene, initPointLight(point3(5, 30, 5),color3(1, 1, 1)));
//
// //addLight(scene, initPointLight(point3(-5, 5, 0), color3(1,1,1)));
// //addLight(scene, initPointLight(point3(0, 5, -5), color3(1,1,1)));
//  
//
// addLight(scene, initPointLight(point3(19, 19, 19), color3(1,1,1)));
// addLight(scene, initPointLight(point3(-19, 19, -19), color3(1,1,1)));
// addLight(scene, initPointLight(point3(-19, 19, 19), color3(1,1,1)));
// addLight(scene, initPointLight(point3(19, 19, -19), color3(1,1,1)));
// addLight(scene, initPointLight(point3(0, 1, -19), color3(1,1,1)));
// 
//  auto light = new AreaLight(vec3(5,4,-5), vec3(-4,0,0), 2, vec3(0,0,-4), 2, color3(1,1,1));
//  addLight(scene, light);
//  //light = new AreaLight(vec3(54,10,18), vec3(-4,0,0), 2, vec3(0,0,-4), 2, color3(1,1,1));
//  //addLight(scene, light);
//  //light = new AreaLight(vec3(18,10,54), vec3(-4,0,0), 2, vec3(0,0,-4), 2, color3(1,1,1));
//  //addLight(scene, light);
//  //light = new AreaLight(vec3(2,50,2), vec3(-4,0,0), 2, vec3(0,0,-4), 2, color3(1,1,1));
//  //addLight(scene, light);
//  //light = new AreaLight(vec3(7,30,7), vec3(-4,0,0), 2, vec3(0,0,-4), 2, color3(1,1,1));
//  //addLight(scene, light);
//
//  Material mat->;
//  mat->.m_diffuseColor = color3(0.f, 0.f, 0.5f);
//  mat->.m_specularColor = color3(0.f, 0.f, 0.7f);
//  mat->.m_roughness = 0.005;
//  mat->.m_IOR = 0.01;
//  mat->.mtype = DIFFUSE;
//
//  //addObject(scene, initSphere(point3(-5, 1.8, -5), 1.8f, mat->lib[8]));
//  //addObject(scene, initSphere(point3(18, 1.8, 26), 1.8f, mat->));
//  //addObject(scene, initSphere(point3(26, 1.8, 18), 1.8f, mat->lib[8]));
//  //
//  glm::mat-> modelMatrix = glm::translate(glm::mat->(1.f), vec3(-5,1.8,-5)) *glm::scale(glm::mat->(1.f), vec3(1.8,1.8,1.8));
//  addObject(scene, initSphere(mat->lib[8], modelMatrix));
//  //glm::mat-> transform1 = glm::translate(glm::mat->(1.f), vec3(3, 2, -7))*glm::scale(glm::mat->(1.f), vec3(1.05,1.05,1.05));
//  //addObject(scene, initSphere(mat->, transform1));
//  //glm::mat-> transform2 = glm::translate(glm::mat->(1.f), vec3(26, 1.8, 18))*glm::scale(glm::mat->(1.f), vec3(1.8,1.8,1.8));
//  //addObject(scene, initSphere(mat->lib[8], transform2));
//
//  //addLight(scene, initPointLight(point3(0, 3.5, 0), color3(1,1,1)));
//
//  auto left = new image_texture("../assets/negx.png");
//  auto right = new image_texture("../assets/posx.png");
//  auto front = new image_texture("../assets/negz.png");
//  auto back = new image_texture("../assets/posz.png");
//  auto up = new image_texture("../assets/posy.png");
//  auto down = new image_texture("../assets/negy.png");
//
//  mat->.m_texture = new CubeMapTexture(left, right, front, back, up, down);
//  mat->.m_diffuseColor = color3(0,0,0);
//  mat->.m_specularColor = color3(0,0,0);
//  glm::mat-> transform = glm::scale(glm::mat->(1.f), vec3(50,50,50));
//  addObject(scene, initCube(mat->, transform));
//  mat->.m_texture = nullptr;
//
//
//  Material mat->
//  mat->m_diffuseColor = color3(0.5);
//  mat->m_specularColor = color3(0.5);
//  mat->m_IOR = 1.1;
//  //mat->m_roughness = 0.0681;
//  mat->m_roughness = 1.2;
//  mat->mtype = DIFFUSE;
//  //mat->m_texture = new checker_texture(color3(0.2, 0.3, 0.1), color3(0.9, 0.9, 0.9));
//  mat->m_texture = new image_texture("../assets/chessboardtexture.png");
//  //mat->m_texture = new CubeMapTexture(color3(1,1,1), color3(1,0,0), color3(1,1,0), color3(0,1,0), color3(0,1,1));
//  //addObject(scene, initPlane(vec3(0, 1, 0), 0, mat->);
//
  return scene;
}

Scene *initScene9() {
  Scene *scene = initScene();
  //setCamera(scene, point3(-0.23, 2.585, 4.3), vec3(-0.23, 2.585, -2.8), vec3(0, 0.01, 0), 60,
  //          float(WIDTH) / float(HEIGHT), 0);
  //setSkyColor(scene, color3(0.2, 0.8, 0.7)); 

  //Material mat->
  //mat->m_diffuseColor = color3(0.5);
  //mat->m_specularColor = color3(0.5);
  //mat->m_IOR = 1.5;
  //mat->m_roughness = 0.0681;
  //mat->mtype = DIFFUSE;
  //addObjectsFromFile("../assets/cornell-box.obj", scene, mat->;

  //addLight(scene, initPointLight(point3(-0.23, 5, -3), color3(1, 1, 1)));
  ////auto light = new AreaLight(vec3(-0.88,5,-3.57), vec3(0.64,0,0), 4, vec3(0,0,-1), 4, color3(1,1,1));
  ////addLight(scene, light);
  return scene;
}

Scene *initScene10() {
  Scene *scene = initScene();
  auto from = point3(0., -60, 12);
  auto at = vec3(0, 0, 12);
  setCameraFOV(scene, from, at, vec3(0, 0.0, 1), 30,
            float(WIDTH), float(HEIGHT), 0.001, glm::distance(from, at));
  setSkyColor(scene, color3(0., 0., 0.)); 

  auto mat = std::make_shared<Blinn>();
  mat->m_diffuseColor = color3(1);
  mat->m_specularColor = color3(0);

  Transform modelMatrix;
  modelMatrix.translate(vec3(0,0,12));
  modelMatrix.translate(vec3(0,0,-12));
  modelMatrix.scale(32,32,32);
  auto ret = new Plane(mat, modelMatrix);
  ret->geom.type = PLANE;
  addObject(scene, ret);

  Transform modelMatrix1;
  modelMatrix1.scale(32,32,32);
  modelMatrix1.rotate(vec3(1,0,0), (180.f)); 
  modelMatrix1.translate(vec3(0,0,12));
  modelMatrix1.translate(vec3(0,0,12));
  ret = new Plane(mat, modelMatrix1);
  ret->geom.type = PLANE;
  addObject(scene, ret);

  Transform modelMatrix2;
  modelMatrix2.scale(32,32,32);
  modelMatrix2.rotate(vec3(1,0,0), (90.f));
  modelMatrix2.translate(vec3(0,0,12));
  modelMatrix2.translate(vec3(0,20,0));
  ret = new Plane(mat, modelMatrix2);
  ret->geom.type = PLANE;
  addObject(scene, ret);

  auto mat1 = std::make_shared<Blinn>();
  mat1->m_diffuseColor = color3(1, 0.5, 0.5);
  mat1->m_specularColor = color3(0);

  Transform modelMatrix3;
  modelMatrix3.scale(32,32,32);
  modelMatrix3.rotate(vec3(0,1,0), (90.f));
  modelMatrix3.translate(vec3(0,0,12));
  modelMatrix3.translate(vec3(-15,0,0));
  ret = new Plane(mat1, modelMatrix3);
  ret->geom.type = PLANE;
  addObject(scene, ret);

  auto mat2 = std::make_shared<Blinn>();
  mat2->m_specularColor = color3(0);
  mat2->m_diffuseColor = color3(0.5, 0.5, 1.0);

  Transform modelMatrix4;
  modelMatrix4.scale(32,32,32);
  modelMatrix4.rotate(vec3(0,1,0), (-90.f));
  modelMatrix4.translate(vec3(0,0,12));
  modelMatrix4.translate(vec3(15,0,0));
  ret = new Plane(mat2, modelMatrix4);
  ret->geom.type = PLANE;
  addObject(scene, ret);

  auto mat3 = std::make_shared<Blinn>();
  mat3->m_diffuseColor = color3(0.8, 0.2, 0.2);
  mat3->m_specularColor = color3(0.7f);
  mat3->m_shininess = 20;
  mat3->m_reflection = vec3(0.7);

  auto mat4 = std::make_shared<Blinn>();
  mat4->m_diffuseColor = color3(0.1, 0.1, 0.9);
  mat4->m_specularColor = color3(0.9f, 0.9, 1.0) * 0.8f;
  mat4->m_IOR = 1.52;
  mat4->m_shininess = 10;
  mat4->m_refraction = vec3(0.8);
  mat4->m_absorption = color3(0.01, 0.001, 0.0001);

  Transform modelMatrix5;
  modelMatrix5.scale(4, 4, 4);
  modelMatrix5.rotate(vec3(0,1,0), (30.f));
  modelMatrix5.translate(vec3(-8,-6,4));
  addObject(scene, initSphere(mat4, modelMatrix5));

  Transform modelMatrix6;
  modelMatrix6.scale(0.8, 0.8, 0.8);
  modelMatrix6.rotate(vec3(0,0,1), (-30.f)); 
  modelMatrix6.translate(vec3(2,5,0));
  addObjectsFromFile("../assets/teapot.obj", scene, mat3, modelMatrix6);

  Transform modelMatrix7;
  modelMatrix7.scale(0.3, 0.3, 0.3);
  modelMatrix7.rotate(vec3(0,0,1), (-60.f));
  modelMatrix7.translate(vec3(5,-6,0));
  addObjectsFromFile("../assets/teapot.obj", scene, mat3, modelMatrix7);


  addLight(scene, initPointLight(point3(0, 0, 22), color3(0.5f)));
  addLight(scene, initAmbientLight(color3(0.1)));
  return scene;
}
 
Scene *initScene11() {
  Scene *scene = initScene();
  auto from = point3(0.0, -70.0, 15.0);
  auto at = vec3(2.0, 0.0, 3.0);
  //setCameraFOV(scene, from, at, vec3(0, 0, 1), 30.f, float(WIDTH), float(HEIGHT), 0.01, distance(from, at));
  setSimpleCamera(scene, from, at, vec3(0.f, 0.f, 1), 30.f, float(WIDTH), float(HEIGHT), 0., 1 );
  setSkyColor(scene, color3(1.,1.,1)); 

  image_texture* sky = new image_texture("../assets/clouds.png");
  scene->m_skyTexture = sky;
  Transform texSky;
  texSky.scale(1,.4,1);
  texSky.translate(vec3(0,-0.1,0));
  scene->m_skyTexture->m_transform = texSky;

  auto mat = std::make_shared<Blinn>();
  mat->m_diffuseColor = color3(1,1,1);
  mat->m_texture = new checker_texture(color3(0.2, 0.2, 0.2), color3(0.6, 0.6, 0.6));
  Transform texT;
  texT.scale(0.01,0.01,0.01);
  mat->m_texture->m_transform = texT;
  mat->m_specularColor = color3(0);
  mat->m_reflection = vec3(0.5);

  Transform modelMatrix;
  modelMatrix.scale(1000, 1000, 1000);
  modelMatrix.rotate(vec3(0,0,1), 30);

  auto ret = new Plane(mat, modelMatrix);
  ret->geom.type = PLANE;
  addObject(scene, ret);
  //addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

  auto mat1 = std::make_shared<Blinn>();
  mat1->m_texture = new image_texture("../assets/bricks.png");
  mat1->m_specularColor = color3(0.3);
  mat1->m_shininess = 10;

  Transform modelMatrix1;
  modelMatrix1.scale(0.8, 0.8, 0.8);
  modelMatrix1.rotate(vec3(0,0,1), -50.f);
  modelMatrix1.translate(vec3(2,-5,0));
  addObjectsFromFile("../assets/teapot.obj", scene, mat1, modelMatrix1);

  auto mat2 = std::make_shared<Blinn>();
  mat2->m_texture = new checker_texture(color3(0.7,0,0), color3(0.3,0,0)); 
  Transform texTransform;
  texTransform.scale(0.25, 0.4,1);
  mat2->m_texture->m_transform = texTransform;
  mat2->m_specularColor = color3(0.8);
  mat2->m_shininess = 100;
  mat2->m_reflection = vec3(0.5f);

  Transform modelMatrix2;
  modelMatrix2.scale(6, 6, 6);
  modelMatrix2.translate(vec3(15, 2, 6));

  addObject(scene, initSphere(mat2, modelMatrix2));

  auto mat3 = std::make_shared<Blinn>();
  mat3->m_diffuseColor = color3(0);
  mat3->m_specularColor = color3(0.8);
  mat3->m_shininess = 100;
  mat3->m_refraction = vec3(1.0);
  mat3->m_IOR = 1.52;

  Transform modelMatrix3;
  modelMatrix3.scale(5, 5, 5);
  modelMatrix3.translate(vec3(-8, -16, 5));

  addObject(scene, initSphere(mat3, modelMatrix3));

  addLight(scene, initAmbientLight(color3(0.2)));
  addLight(scene, initDirectLight(vec3(-1, 0.2, -1), color3(0.6)));
  addLight(scene, initDirectLight(vec3(1, 0.3, -1), color3(0.4)));


  return scene;
}

Scene *initScene12() {
  Scene *scene = initScene();
  auto from = point3(0.0, -70.0, 25.0);
  auto at = vec3(-2.0, 0.0, 3.0);
  //setCameraFOV(scene, from, at, vec3(0, 0, 1), 30.f, float(WIDTH), float(HEIGHT), 0.01, distance(from, at));
  setSimpleCamera(scene, from, at, vec3(0.f, 0.f, 1), 25.f, float(WIDTH), float(HEIGHT), 1.5, 70 );
  setSkyColor(scene, color3(1.,1.,1)); 

  image_texture* sky = new image_texture("../assets/clouds.png");
  scene->m_skyTexture = sky;
  Transform texSky;
  texSky.scale(1,.4,1);
  texSky.translate(vec3(0,-0.1,0));
  scene->m_skyTexture->m_transform = texSky;

  auto mat = std::make_shared<Blinn>();
  mat->m_diffuseColor = color3(0.3);
  mat->m_specularColor = color3(0.1);
  mat->m_shininess = 50;

  Transform modelMatrix;
  modelMatrix.scale(500, 500, 500);

  auto ret = new Plane(mat, modelMatrix);
  ret->geom.type = PLANE;
  addObject(scene, ret);
  //addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

  auto mat1 = std::make_shared<Blinn>();
  mat1->m_texture = new image_texture("../assets/bricks.png");
  mat1->m_specularColor = color3(0.3);
  mat1->m_shininess = 10;

  Transform modelMatrix1;
  modelMatrix1.scale(0.7, 0.7, 0.7);
  modelMatrix1.rotate(vec3(0,0,1), -50.f);
  modelMatrix1.translate(vec3(0,0,0));
  addObjectsFromFile("../assets/teapot.obj", scene, mat1, modelMatrix1);

  auto mat2 = std::make_shared<Blinn>();
  mat2->m_texture = new checker_texture(color3(0.7,0,0), color3(0.3,0,0)); 
  Transform texTransform;
  texTransform.scale(0.25, 0.4,1);
  mat2->m_texture->m_transform = texTransform;
  mat2->m_specularColor = color3(0.8);
  mat2->m_shininess = 100;
  mat2->m_reflection = vec3(0.5f);

  Transform modelMatrix3;
  modelMatrix3.scale(5, 5, 5);
  modelMatrix3.translate(vec3(35, 70, 5));
  addObject(scene, initSphere(mat2, modelMatrix3));

  Transform modelMatrix4;
  modelMatrix4.scale(5, 5, 5);
  modelMatrix4.translate(vec3(30, 60, 5));
  addObject(scene, initSphere(mat2, modelMatrix4));

  Transform modelMatrix5;
  modelMatrix5.scale(5, 5, 5);
  modelMatrix5.translate(vec3(25, 50, 5));
  addObject(scene, initSphere(mat2, modelMatrix5));

  Transform modelMatrix6;
  modelMatrix6.scale(5, 5, 5);
  modelMatrix6.translate(vec3(20, 40, 5));
  addObject(scene, initSphere(mat2, modelMatrix6));

  Transform modelMatrix7;
  modelMatrix7.scale(5, 5, 5);
  modelMatrix7.translate(vec3(15, 30, 5));
  addObject(scene, initSphere(mat2, modelMatrix7));

  Transform modelMatrix8;
  modelMatrix8.scale(5, 5, 5);
  modelMatrix8.translate(vec3(10, 20, 5));
  addObject(scene, initSphere(mat2, modelMatrix8));

  Transform modelMatrix9;
  modelMatrix9.scale(5, 5, 5);
  modelMatrix9.translate(vec3(5, 10, 5));
  addObject(scene, initSphere(mat2, modelMatrix9));

  Transform modelMatrix10;
  modelMatrix10.scale(5, 5, 5);
  modelMatrix10.translate(vec3(-5, -10, 5));
  addObject(scene, initSphere(mat2, modelMatrix10));

  Transform modelMatrix11;
  modelMatrix11.scale(5, 5, 5);
  modelMatrix11.translate(vec3(-10, -20, 5));
  addObject(scene, initSphere(mat2, modelMatrix11));

  Transform modelMatrix12;
  modelMatrix12.scale(5, 5, 5);
  modelMatrix12.translate(vec3(-15, -30, 5));
  addObject(scene, initSphere(mat2, modelMatrix12));


  addLight(scene, initAmbientLight(color3(0.2)));
  addLight(scene, initDirectLight(vec3(-1, 0.2, -1), color3(0.6)));
  addLight(scene, initDirectLight(vec3(1, 0.3, -1), color3(0.4)));


  return scene;
}

Scene *initScene13() {
  Scene *scene = initScene();
  auto from = point3(42, -42.0, 15.0);
  auto at = vec3(6.0, 0.0, -3.0);
  //setCameraFOV(scene, from, at, vec3(0, 0, 1), 30.f, float(WIDTH), float(HEIGHT), 0.01, distance(from, at));
  setSimpleCamera(scene, from, at, vec3(0.f, 0.f, 1), 40.f, float(WIDTH), float(HEIGHT), 0, 1 );
  setSkyColor(scene, color3(1.,1.,1)); 

  image_texture* sky = new image_texture("../assets/clouds.png");
  scene->m_skyTexture = sky;
  Transform texSky;
  texSky.scale(1,.4,1);
  texSky.translate(vec3(0,-0.1,0));
  scene->m_skyTexture->m_transform = texSky;

  auto mat = std::make_shared<Blinn>();
  mat->m_diffuseColor = color3(1,1,1);
  mat->m_texture = new checker_texture(color3(0.5, 0.5, 0.7), color3(1));
  Transform texCheck;
  texCheck.scale(0.003, 0.003, 0.003);
  mat->m_texture->m_transform = texCheck;
  mat->m_specularColor = color3(0);

  Transform modelMatrix;
  modelMatrix.scale(1000, 1000, 1000);

  auto ret = new Plane(mat, modelMatrix);
  ret->geom.type = PLANE;
  addObject(scene, ret);
  //addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);
  
  auto matWall = std::make_shared<Blinn>();
  matWall->m_texture = new image_texture("../assets/bricks.png");
  Transform texCheckWall;
  texCheckWall.scale(0.05, 0.05, 0.05);
  matWall->m_texture->m_transform = texCheckWall;
  matWall->m_specularColor = color3(0.3);
  matWall->m_shininess = 10;

  Transform modelMatrixWall;
  modelMatrixWall.scale(100, 100, 100);
  modelMatrixWall.rotate(vec3(1,0,0), 90);
  modelMatrixWall.translate(vec3(0,10,0));

  auto retWall = new Plane(matWall, modelMatrixWall);
  retWall->geom.type = PLANE;
  addObject(scene, retWall);
  //addObjectsFromFile("../assets/plane.obj", scene, mat, modelMatrix);

  auto mat1 = std::make_shared<Blinn>();
  mat1->m_diffuseColor = color3(0.8, 0.2, 0.2);
  mat1->m_specularColor = color3(0.8);
  mat1->m_shininess = 100;
  mat1->m_reflection = vec3(0.3);

  Transform modelMatrix1;
  modelMatrix1.scale(0.5, 0.5, 0.5);
  modelMatrix1.rotate(vec3(0,0,1), -50.f);
  modelMatrix1.translate(vec3(13,-21,0));
  addObjectsFromFile("../assets/teapot.obj", scene, mat1, modelMatrix1);

  auto mat2 = std::make_shared<Blinn>();
  mat2->m_diffuseColor = color3(0);
  mat2->m_specularColor = color3(0.8);
  mat2->m_shininess = 100;
  mat2->m_reflection = vec3(0.8f);

  Transform modelMatrix3;
  modelMatrix3.scale(6, 6, 6);
  modelMatrix3.translate(vec3(-28, 0, 6));
  addObject(scene, initSphere(mat2, modelMatrix3));

  auto mat3 = std::make_shared<Blinn>();
  mat3->m_diffuseColor = color3(0);
  mat3->m_specularColor = color3(0.8);
  mat3->m_shininess = 50;
  mat3->m_reflection = vec3(0.8f);
  mat3->m_reflectionGloss = 0.05;

  Transform modelMatrix4;
  modelMatrix4.scale(6, 6, 6);
  modelMatrix4.translate(vec3(-14, 0, 6));
  addObject(scene, initSphere(mat3, modelMatrix4));

  auto mat4 = std::make_shared<Blinn>();
  mat4->m_diffuseColor = color3(0);
  mat4->m_specularColor = color3(0.8);
  mat4->m_shininess = 20;
  mat4->m_reflection = vec3(0.8f);
  mat4->m_reflectionGloss = 0.1;

  Transform modelMatrix5;
  modelMatrix5.scale(6, 6, 6);
  modelMatrix5.translate(vec3(0, 0, 6));
  addObject(scene, initSphere(mat4, modelMatrix5));

  auto mat5 = std::make_shared<Blinn>();
  mat5->m_diffuseColor = color3(0);
  mat5->m_specularColor = color3(0.8);
  mat5->m_shininess = 10;
  mat5->m_reflection = vec3(0.8f);
  mat5->m_reflectionGloss = 0.2;

  Transform modelMatrix6;
  modelMatrix6.scale(6, 6, 6);
  modelMatrix6.translate(vec3(14, 0, 6));
  addObject(scene, initSphere(mat5, modelMatrix6));

  auto mat6 = std::make_shared<Blinn>();
  mat6->m_diffuseColor = color3(0);
  mat6->m_specularColor = color3(0.8);
  mat6->m_shininess = 50;
  mat6->m_refraction = vec3(0.8);
  mat6->m_IOR = 2.0;
  mat6->m_refractionGloss = 0.05;

  Transform modelMatrix7;
  modelMatrix7.scale(2, 2, 2);
  modelMatrix7.translate(vec3(22, -30, 2));
  addObject(scene, initSphere(mat6, modelMatrix7));

  addLight(scene, initAmbientLight(color3(0.2)));

  addLight(scene, initPointLight(vec3(-50, -100, 50), vec3(1,1,1), 5));

  return scene;
}


Scene* parseScene(int sceneId){
  Scene *scene = NULL;
  switch (sceneId) {
  case 0:
    scene = initScene0();
    break;
  case 1:
    scene = initScene1();
    break;
  case 2:
    scene = initScene2();
    break;
  case 3:
    scene = initScene3();
    break;
  case 4:
    scene = initScene4();
    break;
  case 5:
    scene = initScene5();
    break;
  case 6:
    scene = initScene6();
    break;
  case 7:
    scene = initScene7();
    break;
  case 8:
    scene = initScene8();
    break;
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

  default:
    scene = initScene0();
    break;
  }
  return scene;
}

