#include "defines.h"
#include "image.h"
#include "ray.h"
#include "raytracer.h"
#include "scene.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <random>

#define WIDTH 600
#define HEIGHT 400

Material mat_lib[] = {
    /* nickel */
    {2.4449, 0.0681, {1.0, 0.882, 0.786}, {0.014, 0.012, 0.012}, LAMBERTIAN},

    /* specular black phenolic */
    {1.072, 0.0588, {1.0, 0.824, 0.945}, {0.002, 0.002, 0.003}, LAMBERTIAN},

    /* specular blue phenolic */
    {1.1051, 0.0568, {0.005, 0.013, 0.032}, {1.0, 0.748, 0.718}, LAMBERTIAN},

    /* specular green phenolic */
    {1.1051, 0.0567, {0.006, 0.026, 0.022}, {1.0, 0.739, 0.721}, LAMBERTIAN},

    /* specular white phenolic */
    {1.1022, 0.0579, {0.286, 0.235, 0.128}, {1.0, 0.766, 0.762}, LAMBERTIAN},

    /* marron plastic */
    {1.0893, 0.0604, {0.202, 0.035, 0.033}, {1.0, 0.857, 0.866}, LAMBERTIAN},

    /* purple paint */
    {1.1382, 0.0886, {0.301, 0.034, 0.039}, {1.0, 0.992, 0.98}, LAMBERTIAN},

    /* red specular plastic */
    {1.0771, 0.0589, {0.26, 0.036, 0.014}, {1.0, 0.852, 1.172}, LAMBERTIAN},

    /* green acrylic */
    {1.1481, 0.0625, {0.016, 0.073, 0.04}, {1.0, 1.056, 1.146}, LAMBERTIAN},

    /* blue acrylic */
    {1.1153, 0.068, {0.012, 0.036, 0.106}, {1.0, 0.965, 1.07}, LAMBERTIAN}};

Scene *initScene0() {
  Scene *scene = initScene();
  setCamera(scene, point3(3.0f, 1.0f, 0.0f), vec3(0.0f, 0.3f, 0.0f), vec3(0.0f, 1.0f, 0.0f), 40,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.1f, 0.3f, 0.5f));
  Material mat;
  mat.IOR = 1.5;
  mat.roughness = 0.1;
  mat.specularColor = color3(0.5f);

  mat.type = METAL;
  mat.diffuseColor = color3(0.8f, 0.0f, 0.0f);
  mat.fuzz = 0.3;
  addObject(scene, initSphere(point3(0, 0.25, 0), 0.25, mat));

  mat.type = DIELECTRIC;
  mat.diffuseColor = color3(0.5f, 0.f, 0.f);
  addObject(scene, initSphere(point3(0.50, 0.25, 0), .25, mat));

  mat.type = LAMBERTIAN;
  mat.diffuseColor = color3(0.f, 0.5f, 0.5f);
  addObject(scene, initSphere(point3(0, 1, 0), .25, mat));

  mat.diffuseColor = color3(0.f, 0.f, 0.5f);
  addObject(scene, initSphere(point3(0, 0.25, 0.50), .25, mat));

  mat.diffuseColor = color3(0.8f, 0.8f, 0.0f);
  addObject(scene, initSphere(point3(0.f, -10.f, 0.f), 10.f, mat));
  //addObject(scene, initPlane(vec3(0, 1, 0), 0, mat));

  addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initLight(point3(4, 10, -2), color3(1, 1, 1)));

  return scene;
}

Scene *initScene1() {
  Scene *scene = initScene();
  setCamera(scene, point3(3, 0, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.2, 0.2, 0.7));

  Material mat;
  mat.IOR = 1.12;
  mat.roughness = 0.2;
  mat.specularColor = color3(0.4f);
  mat.diffuseColor = color3(0.6f);
  mat.type = LAMBERTIAN;

  for (int i = 0; i < 10; ++i) {
    mat.diffuseColor = color3(0.301, 0.034, 0.039);
    mat.specularColor = color3(1.0, 0.992, 0.98);
    mat.IOR = 1.1382;
    mat.roughness = 0.0886;
    mat.roughness = ((float)10 - i) / (10 * 9.f);
    addObject(scene, initSphere(point3(0, 0, -1.5 + i / 9.f * 3.f), .15, mat));
  }
  for (int i = 0; i < 10; ++i) {
    mat.diffuseColor = color3(0.012, 0.036, 0.106);
    mat.specularColor = color3(1.0, 0.965, 1.07);
    mat.IOR = 1.1153;
    mat.roughness = 0.068;
    mat.roughness = ((float)i + 1) / 10.f;
    addObject(scene, initSphere(point3(0, 1, -1.5 + i / 9.f * 3.f), .15, mat));
  }
  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(1.0, 0.882, 0.786);
  mat.IOR = 2.4449;
  mat.roughness = 0.0681;
  addObject(scene, initSphere(point3(-3.f, 1.f, 0.f), 2., mat));

  mat.diffuseColor = color3(0.016, 0.073, 0.04);
  mat.specularColor = color3(1.0, 1.056, 1.146);
  mat.IOR = 1.1481;
  mat.roughness = 0.0625;
  addObject(scene, initPlane(vec3(0, 1, 0), +1, mat));

  addLight(scene, initLight(point3(10, 10, 10), color3(10, 10, 10)));
  addLight(scene, initLight(point3(4, 10, -2), color3(5, 3, 10)));
  return scene;
}

Scene *initScene2() {
  Scene *scene = initScene();
  setCamera(scene, point3(0.5, 3, 1), vec3(0, 0, 0.6), vec3(0, 0, 1), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.2, 0.2, 0.7));
  Material mat;
  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(1.0, 0.882, 0.786);
  mat.IOR = 2.4449;
  mat.roughness = 0.0681;
  mat.type = LAMBERTIAN;

  mat.diffuseColor = color3(0.05, 0.05, 0.05);
  mat.specularColor = color3(0.95);
  mat.IOR = 1.1022;
  mat.roughness = 0.0579;

  addObject(scene, initPlane(vec3(0, 0, 1), 0, mat));

  mat.diffuseColor = color3(0.005, 0.013, 0.032);
  mat.specularColor = color3(1.0, 0.748, 0.718);
  for (int i = 0; i < 4; ++i) {
    mat.IOR = 1.1051 + (-0.1 + float(i) / 3.f * 0.4);
    for (int j = 0; j < 10; ++j) {
      mat.roughness = 0.0568 + (-0.1 + float(j) / 9.f * 0.3);
      addObject(scene, initSphere(point3(-1.5 + float(j) / 9.f * 3.f, 0,
                                         0.4 + float(i) * 0.4f),
                                  .15, mat));
    }
  }

  addLight(scene, initLight(point3(-20, 5, 10), color3(30, 30, 30)));
  addLight(scene, initLight(point3(10, 10, 10), color3(30, 30, 30)));
  addLight(scene, initLight(point3(50, -100, 10), color3(1, 0.7, 2)));
  return scene;
}

Scene *initScene3() {
  Scene *scene = initScene();
  setCamera(scene, point3(4.5, .8, 4.5), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.2, 0.2, 0.7));
  Material mat;
  mat.diffuseColor = color3(0.301, 0.034, 0.039);
  mat.specularColor = color3(1.0, 0.992, 0.98);
  mat.IOR = 1.1382;
  mat.roughness = 0.0886;
  mat.type = LAMBERTIAN;

  addLight(scene, initLight(point3(0, 1.7, 1), .5f * color3(3, 3, 3)));
  addLight(scene, initLight(point3(3, 2, 3), .5f * color3(4, 4, 4)));
  addLight(scene, initLight(point3(4, 3, -1), .5f * color3(5, 5, 5)));

  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(0.7, 0.882, 0.786);
  mat.IOR = 6;
  mat.roughness = 0.0181;
  addObject(scene, initSphere(point3(0, 0.1, 0), .3, mat));

  mat.diffuseColor = color3(0.26, 0.036, 0.014);
  mat.specularColor = color3(1.0, 0.852, 1.172);
  mat.IOR = 1.3771;
  mat.roughness = 0.01589;
  addObject(scene, initSphere(point3(1, -.05, 0), .15, mat));

  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(0.7, 0.882, 0.786);
  mat.IOR = 3;
  mat.roughness = 0.00181;
  addObject(scene, initSphere(point3(3, 0.05, 2), .25, mat));

  mat.diffuseColor = color3(0.46, 0.136, 0.114);
  mat.specularColor = color3(0.8, 0.852, 0.8172);
  mat.IOR = 1.5771;
  mat.roughness = 0.01589;
  addObject(scene, initSphere(point3(1.3, 0., 2.6), 0.215, mat));

  mat.diffuseColor = color3(0.06, 0.26, 0.22);
  mat.specularColor = color3(0.70, 0.739, 0.721);
  mat.IOR = 1.3051;
  mat.roughness = 0.567;
  addObject(scene, initSphere(point3(1.9, 0.05, 2.2), .25, mat));

  mat.diffuseColor = color3(0.012, 0.036, 0.406);
  mat.specularColor = color3(1.0, 0.965, 1.07);
  mat.IOR = 1.1153;
  mat.roughness = 0.068;
  mat.roughness = 0.18;
  addObject(scene, initSphere(point3(0, 0, 1), .20, mat));

  mat.diffuseColor = color3(.2, 0.4, .3);
  mat.specularColor = color3(.2, 0.2, .2);
  mat.IOR = 1.382;
  mat.roughness = 0.05886;
  addObject(scene, initPlane(vec3(0, 1, 0), 0.2, mat));

  mat.diffuseColor = color3(.5, 0.09, .07);
  mat.specularColor = color3(.2, .2, .1);
  mat.IOR = 1.8382;
  mat.roughness = 0.886;
  addObject(scene, initPlane(vec3(1, 0.0, -1.0), 2, mat));

  mat.diffuseColor = color3(0.1, 0.3, .05);
  mat.specularColor = color3(.5, .5, .5);
  mat.IOR = 1.9382;
  mat.roughness = 0.0886;
  addObject(scene, initPlane(vec3(0.3, -0.2, 1), 4, mat));
  return scene;
}

Scene *initScene4() {
  Scene *scene = initScene();
  setCamera(scene, point3(6, 4, 4), vec3(0, 1, 0), vec3(0, 1, 0), 90,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.2, 0.2, 0.7));
  Material mat;
  mat.diffuseColor = color3(0.301, 0.034, 0.039);
  mat.specularColor = color3(1.0, 0.992, 0.98);
  mat.IOR = 1.1382;
  mat.roughness = 0.0886;
  mat.type = LAMBERTIAN;

  addLight(scene, initLight(point3(10, 10.7, 1), .5f * color3(3, 3, 3)));
  addLight(scene, initLight(point3(8, 20, 3), .5f * color3(4, 4, 4)));
  addLight(scene, initLight(point3(4, 30, -1), .5f * color3(5, 5, 5)));

  mat.diffuseColor = color3(.2, 0.4, .3);
  mat.specularColor = color3(.2, 0.2, .2);
  mat.IOR = 2.382;
  mat.roughness = 0.005886;
  //addObject(scene, initPlane(vec3(0, 1, 0), 0.2, mat));

  mat.diffuseColor = color3(.5, 0.09, .07);
  mat.specularColor = color3(.2, .2, .1);
  mat.IOR = 2.8382;
  mat.roughness = 0.00886;
  //addObject(scene, initPlane(vec3(1, 0.0, 0.0), 2, mat));

  mat.diffuseColor = color3(0.1, 0.3, .05);
  mat.specularColor = color3(.5, .5, .5);
  mat.IOR = 2.9382;
  mat.roughness = 0.00886;
  //addObject(scene, initPlane(vec3(0, 0, 1), 4, mat));

  for (int i = 0; i < 600; i++) {
    addObject(scene,
              initSphere(point3(1 + rand() % 650 / 100.0, rand() % 650 / 100.0,
                                1 + rand() % 650 / 100.0),
                         .05 + rand() % 200 / 1000.0, mat_lib[rand() % 10]));
  }
  return scene;
}

Scene *initScene5() {

  Scene *scene = initScene();
  setCamera(scene, point3(3, 1, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.1f, 0.3f, 0.5f));
  Material mat;
  mat.diffuseColor = color3(0.301, 0.034, 0.039);
  mat.specularColor = color3(1.0, 0.992, 0.98);
  mat.IOR = 1.1382;
  mat.roughness = 0.0886;
  mat.type = LAMBERTIAN;
  addLight(scene, initLight(point3(0, 1.7, 1), .5f * color3(3, 3, 3)));
  addLight(scene, initLight(point3(3, 2, 3), .5f * color3(4, 4, 4)));
  addLight(scene, initLight(point3(4, 3, -1), .5f * color3(5, 5, 5)));
  addLight(scene, initLight(point3(1, 0, 1), .5f * color3(3, 3, 3)));

  vec3 v0 = vec3(1,0,0);
  vec3 v1 = vec3(0,1,0);
  vec3 v2 = vec3(0,0,1);

  vec3 v1v0 = v1-v0;
  vec3 v2v0 = v2-v0;
  vec3 n = normalize(cross(v1v0, v2v0));

  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(0.7, 0.882, 0.786);
  mat.IOR = 6;
  mat.roughness = 0.0181;
  vec2 textures[3];
  addObject(scene, initTriangle(v0, v1, v2, n,textures,mat));

  v0 = vec3(1,0,0);
  v1 = vec3(0,1,0);
  v2 = vec3(0,0,-1);

  v1v0 = v1-v0;
  v2v0 = v2-v0;
  n = normalize(cross(v1v0, v2v0));
  addObject(scene, initTriangle(v0, v1, v2, n,textures,mat));
  
  mat.diffuseColor = color3(0.016, 0.073, 0.04);
  mat.specularColor = color3(1.0, 1.056, 1.146);
  mat.IOR = 1.1481;
  mat.roughness = 0.0625;
  addObject(scene, initPlane(vec3(0, 1, 0), +1, mat));

  addLight(scene, initLight(point3(10, 100, 100), color3(50, 50, 50)));
  return scene;
}

void addObjectsFromFile(const char *filename, Scene *scene, Material default_mat){
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;
  Material mat;

  readObjToTriangleMesh(filename, attrib, shapes, materials);

  for(size_t s = 0; s < shapes.size(); s++){
    size_t index_offset = 0;
    for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size();f++){
      int fv = shapes[s].mesh.num_face_vertices[f];

      std::vector<point3> vector;
      vec2 textures[fv];

      // Loop over vertices in the face.
      for (int v = 0; v < fv; v++) {
        // access to vertex
        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
        tinyobj::real_t vx = attrib.vertices[3*idx.vertex_index+0];
        tinyobj::real_t vy = attrib.vertices[3*idx.vertex_index+1];
        tinyobj::real_t vz = attrib.vertices[3*idx.vertex_index+2];
        //tinyobj::real_t nx = attrib.normals[3*idx.normal_index+0];
        //tinyobj::real_t ny = attrib.normals[3*idx.normal_index+1];
        //tinyobj::real_t nz = attrib.normals[3*idx.normal_index+2];
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
        vector.push_back(v0);
      }

      index_offset += fv;

      vec3 v0 = vector[0];
      vec3 v1 = vector[1];
      vec3 v2 = vector[2];

      vec3 v1v0 = v1-v0;
      vec3 v2v0 = v2-v0;
      vec3 n = normalize(cross(v1v0, v2v0));

      if(!materials.empty()){
        int matIndex = shapes[s].mesh.material_ids[f];
        mat.IOR = materials[matIndex].ior;
        mat.roughness = materials[matIndex].ior;
        mat.diffuseColor.r = materials[matIndex].diffuse[0];
        mat.diffuseColor.g = materials[matIndex].diffuse[1];
        mat.diffuseColor.b = materials[matIndex].diffuse[2];
        if(materials[matIndex].specular[0] == 0 && materials[matIndex].specular[1] == 0 && materials[matIndex].specular[1] == 0){
          mat.specularColor.r = materials[matIndex].diffuse[0];
          mat.specularColor.g = materials[matIndex].diffuse[1];
          mat.specularColor.b = materials[matIndex].diffuse[2];
        }
        else{
          mat.specularColor.r = materials[matIndex].specular[0];
          mat.specularColor.g = materials[matIndex].specular[1];
          mat.specularColor.b = materials[matIndex].specular[2];
        }
      }
      else{
        mat = default_mat;
      }
      addObject(scene, initTriangle(v0, v1, v2, n,textures,mat));
    }
  }
}

Scene *initScene6() {

  Scene *scene = initScene();
  setCamera(scene, point3(-2.5, 1.5f, -2), vec3(0, 1, 0), vec3(0, 1, 0), 40,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.6f, 0.3f, 0.5f));
  
  //addLight(scene, initLight(point3(100, 1, 1), color3(50, 50, 50)));
  //addLight(scene, initLight(point3(1, 0.1, 0.5), color3(5, 5, 5)));
  //addLight(scene, initLight(point3(0, 0, 0), color3(5, 5, 5)));
  addLight(scene, initLight(point3(-6, 3, 7), color3(3, 3, 3)));

  addLight(scene, initLight(point3(6, 3, -7), color3(3, 3, 3)));

  addLight(scene, initLight(point3(-6, 4, 0), color3(3, 3, 3)));
  addLight(scene, initLight(point3(4, 10, 10), color3(1, 1, 1)));



  Material mat;
  mat.type = LAMBERTIAN;
  mat.diffuseColor = color3(0.5f);

  addObjectsFromFile("../assets/bunny.obj", scene, mat);

  mat.type = LAMBERTIAN;
  mat.diffuseColor = color3(0.7, 0, 0);
  mat.specularColor = color3(1, 0, 0);
  mat.IOR = 1.5;
  mat.roughness = 0.0681;
  //addObject(scene, initPlane(vec3(0, 1, 0), 0.01, mat));
  addObject(scene, initSphere(point3(0, -10, 0), 10, mat));

  mat.diffuseColor = color3(0, 0.8, 0);
  mat.specularColor = color3(0.6, 0.6, 0);
  mat.IOR = 1.5;
  mat.roughness = 0.0681;
  //addObject(scene, initPlane(vec3(0, 0, 1), +10, mat));

  
  return scene;
}

Scene *initScene7() {

  Scene *scene = initScene();
  setCamera(scene, point3(3, 3, 0), vec3(0, 3, -2), vec3(0, 1, 0), 60,
            float(WIDTH) / float(HEIGHT));
  setSkyColor(scene, color3(0.2, 0.2, 0.7)); 


  Material mat;
  mat.diffuseColor = color3(0.05, 0.012, 0.6);
  mat.specularColor = color3(1.0, 0.882, 0.786);
  mat.IOR = 1.3;
  mat.roughness = 0.0681;
  mat.type = LAMBERTIAN;

  addObjectsFromFile("../assets/werewolf.obj", scene, mat);
  
  mat.diffuseColor = color3(.4, 0.8, .4);
  mat.specularColor = color3(.4, 0.6, .2);
  mat.IOR = 1.382;
  mat.roughness = 0.05886;
  addObject(scene, initPlane(vec3(100, 1, 0), +100, mat));
  

  addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initLight(point3(4, 10, -2), color3(1, 1, 1)));
  addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initLight(point3(0, 40, 0), color3(1, 1, 1)));
  return scene;
}

Scene *initScene9() {
  Scene *scene = initScene();
  setCamera(scene, point3(3, 1, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0, 0, 0));
  Material mat;
  mat.IOR = 20;
  mat.roughness = 0.001;
  mat.specularColor = color3(0.9f);
  mat.type = LAMBERTIAN;

  mat.diffuseColor = color3(0.9f, 0.9f, 0.9f);

  addObject(scene, initPlane(vec3(0, 1, 0), 0, mat));
  addObject(scene, initPlane(vec3(0.5, 0, -0.5), 0, mat));
  addObject(scene, initPlane(vec3(-0.5, 0, -0.5), 0, mat));

  addObject(scene, initPlane(vec3(0.5, 0, -0.5), -2.2, mat));
  addObject(scene, initPlane(vec3(-0.5, 0, -0.5), +2.2, mat));

  mat = mat_lib[0];
  addObject(scene, initSphere(point3(1, 0.5, 0), .25, mat));
  
  addLight(scene, initLight(point3(1, 3, 0), color3(1, 1, 1)));
  addLight(scene, initLight(point3(1.1, 3, 0), color3(1, 1, 1)));
  //addLight(scene, initLight(point3(0.9, 1, 0), color3(1, 1, 1)));

  return scene;
}

static std::mt19937 m_rnGenerator{};
static std::uniform_real_distribution<float> m_unifDistribution{0.0f, 1.0f};
static std::uniform_real_distribution<float> m_unifDistributionRange{0.5f, 1.0f};
static std::uniform_real_distribution<float> m_unifDistributionRange2{0.0f, 0.5f};

Scene *initScene8() {
  Scene *scene = initScene();
  setCamera(scene, point3(13, 2, 3), vec3(0, 0, 0), vec3(0, 1, 0), 20,
            3.0 / 2.0);
  setSkyColor(scene, color3(0, 0, 0));
  Material ground_material;
  ground_material.type = LAMBERTIAN;
  ground_material.diffuseColor = color3(0.5,0.5,0.5);

  addObject(scene, initSphere(point3(0, -1000, 0), 1000, ground_material));

  for(int a = -11; a < 11; a++){
    for(int b = -11; b < 11; b++){
      auto choose_mat = m_unifDistribution(m_rnGenerator);
      point3 center(a+ 0.9*m_unifDistribution(m_rnGenerator), 0.2, b + 0.9*m_unifDistribution(m_rnGenerator));

      if (glm::length((center - point3(4, 0.2, 0))) > 0.9f){
        Material sphere_material;
        
        if(choose_mat < 0.8){
          auto albedo = color3(m_unifDistribution(m_rnGenerator), m_unifDistribution(m_rnGenerator), m_unifDistribution(m_rnGenerator)) * color3(m_unifDistribution(m_rnGenerator), m_unifDistribution(m_rnGenerator), m_unifDistribution(m_rnGenerator));  
          sphere_material.type = LAMBERTIAN;
          sphere_material.diffuseColor = albedo;
          addObject(scene, initSphere(center, 0.2, sphere_material));
        } else if(choose_mat < 0.95){
          auto albedo = color3(m_unifDistributionRange(m_rnGenerator),m_unifDistributionRange(m_rnGenerator),m_unifDistributionRange(m_rnGenerator)); 
          auto fuzz = m_unifDistributionRange2(m_rnGenerator);
          sphere_material.type = METAL;
          sphere_material.fuzz = fuzz;
          sphere_material.diffuseColor = albedo;
          addObject(scene, initSphere(center, 0.2, sphere_material));
        } else {
          sphere_material.type = DIELECTRIC;
          sphere_material.IOR = 1.5;
          addObject(scene, initSphere(center, 0.2, sphere_material));
        }
      }
    }
  }

  Material material1;
  material1.type = DIELECTRIC;
  material1.IOR = 1.5;
  addObject(scene, initSphere(point3(0, 1, 0), 1.0, material1));

  Material material2;
  material2.type = LAMBERTIAN;
  material2.diffuseColor = color3(0.4, 0.2, 0.1);
  addObject(scene, initSphere(point3(-4, 1, 0), 1.0, material2));
  
  Material material3;
  material3.type = METAL;
  material3.fuzz = 0.0;
  material3.diffuseColor = color3(0.7, 0.6, 0.5);
  addObject(scene, initSphere(point3(4, 1, 0), 1.0, material3));
  
  return scene;
}

int main(int argc, char *argv[]) {
  //printf("Welcome to the L3 IGTAI RayTracer project\n");

  char basename[256];

  if (argc < 2 || argc > 3) {
    printf("usage : %s filename i\n", argv[0]);
    printf("        filename : where to save the result, whithout extention\n");
    printf("        i : scenen number, optional\n");
    exit(0);
  }

  strncpy(basename, argv[1], 255);

  int scene_id = 0;
  if (argc == 3) {
    scene_id = atoi(argv[2]);
  }
  const auto aspect_ratio = 3.0/2.0;
  const int image_width = 1200;
  const int image_height = static_cast<int>(image_width / aspect_ratio); 

  Image *img = initImage(image_width, image_height);
  //Image *img = initImage(WIDTH, HEIGHT);
  Scene *scene = NULL;
  switch (scene_id) {
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

  default:
    scene = initScene0();
    break;
  }

  //printf("render scene %d\n", scene_id);

  renderImage(img, scene);
  freeScene(scene);
  scene = NULL;

  printf("save image to %s\n", basename);
  saveImage(img, basename);
  freeImage(img);
  img = NULL;
  printf("done. Goodbye\n");

  return 0;
}
