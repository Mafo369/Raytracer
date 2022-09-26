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

Material mat_lib[] = {
    /* bunny glass */
    {1.05, 2.2, {1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, DIELECTRIC },

    /* specular black phenolic */
    {1.072, 0.0588, {1.0, 0.824, 0.945}, {0.002, 0.002, 0.003}},

    /* specular blue phenolic */
    {1.1051, 0.0568, {0.005, 0.013, 0.032}, {1.0, 0.748, 0.718}},

    /* specular green phenolic */
    {1.1051, 0.0567, {0.006, 0.026, 0.022}, {1.0, 0.739, 0.721}},

    /* specular white phenolic */
    {1.1022, 0.0579, {0.286, 0.235, 0.128}, {1.0, 0.766, 0.762}},

    /* marron plastic */
    {1.0893, 0.0604, {0.202, 0.035, 0.033}, {1.0, 0.857, 0.866}},

    /* purple paint */
    {1.1382, 0.0886, {0.301, 0.034, 0.039}, {1.0, 0.992, 0.98}},

    /* red specular plastic */
    {1.0771, 0.0589, {0.26, 0.036, 0.014}, {1.0, 0.852, 1.172}},

    /* green acrylic */
    {1.1481, 0.0625, {0.016, 0.073, 0.04}, {1.0, 1.056, 1.146}},

    /* blue acrylic */
    {1.1153, 0.068, {0.012, 0.036, 0.106}, {1.0, 0.965, 1.07}},

    {1.51, 2.2, {1.0, 1.0, 1.0}, {1.0, .0, .0}, DIELECTRIC },
};



Scene *initScene0() {
  Scene *scene = initScene();
  setCamera(scene, point3(3, 1, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.1f, 0.3f, 0.5f));
  Material mat;
  mat.IOR = 1.3;
  mat.roughness = 0.1;
  mat.specularColor = color3(0.5f);

  mat.diffuseColor = color3(.5f);
  addObject(scene, initSphere(point3(0, 0, 0), 0.25, mat));

  mat.diffuseColor = color3(0.5f, 0.f, 0.f);
  addObject(scene, initSphere(point3(1, 0, 0), .25, mat));

  mat.diffuseColor = color3(0.f, 0.5f, 0.5f);
  addObject(scene, initSphere(point3(0, 1, 0), .25, mat));

  mat.diffuseColor = color3(0.f, 0.f, 0.5f);
  addObject(scene, initSphere(point3(0, 0, 1), .25, mat));

  mat.diffuseColor = color3(0.6f);
  mat.m_texture = nullptr;
  addObject(scene, initPlane(vec3(0, 1, 0), 0, mat));

  addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initLight(point3(4, 10, -2), color3(1, 1, 1)));

  return scene;
}

Scene *initScene1() {

  Scene *scene = initScene();
  setCamera(scene, point3(-9, 0, 0), vec3(0, 0.3, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.2, 0.2, 0.7));

  Material mat;
  mat.IOR = 1.12;
  mat.roughness = 2.0;
  mat.specularColor = color3(0.4f);
  mat.diffuseColor = color3(0.6f);

  for (int i = 0; i < 10; ++i) {
    mat.diffuseColor = color3(0.301, 0.034, 0.039);
    mat.specularColor = color3(1.0, 0.992, 0.98);
    mat.IOR = 1.1382;
    mat.roughness = 0.0886;
    mat.roughness = ((float)10 - i) / (10 * 9.f);
    addObject(scene, initSphere(point3(0, 0, -1.5 + i / 9.f * 3.f), .15, mat_lib[i]));
  }
  for (int i = 0; i < 10; ++i) {
    mat.diffuseColor = color3(0.012, 0.036, 0.106);
    mat.specularColor = color3(1.0, 0.965, 1.07);
    mat.IOR = 1.1153;
    mat.roughness = 0.068;
    mat.roughness = ((float)i + 1) / 10.f;
    addObject(scene, initSphere(point3(0, 1, -1.5 + i / 9.f * 3.f), .15, mat_lib[i]));
  }
  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(1.0, 0.882, 0.786);
  mat.IOR = 2.4449;
  mat.roughness = 0.0681;
  addObject(scene, initSphere(point3(-3.f, 1.f, 0.f), 2., mat_lib[0]));
  
  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(1.0, 0.882, 0.786);
  mat.IOR = 2.4449;
  mat.roughness = 0.0681;
  addObject(scene, initSphere(point3(-12.f, 1.f, 0.f), 0.8, mat_lib[1]));

  mat.diffuseColor = color3(0.016, 0.073, 0.04);
  mat.specularColor = color3(1.0, 1.056, 1.146);
  mat.IOR = 1.1481;
  mat.roughness = 0.0625;
  addObject(scene, initPlane(vec3(0, 1, 0), +1, mat));

  addLight(scene, initLight(point3(10, 10.7, 1), .5f * color3(3, 3, 3)));
  addLight(scene, initLight(point3(8, 20, 3), .5f * color3(4, 4, 4)));
  addLight(scene, initLight(point3(4, 30, -1), .5f * color3(5, 5, 5)));
  //addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  //addLight(scene, initLight(point3(4, 10, -2), color3(1, 1, 1)));
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
  setCamera(scene, point3(6, 4, 6), vec3(0, 1, 0), vec3(0, 1, 0), 90,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.2, 0.2, 0.7));
  Material mat;
  mat.diffuseColor = color3(0.301, 0.034, 0.039);
  mat.specularColor = color3(1.0, 0.992, 0.98);
  mat.IOR = 1.1382;
  mat.roughness = 0.0886;

  addLight(scene, initLight(point3(10, 10.7, 1), .5f * color3(3, 3, 3)));
  addLight(scene, initLight(point3(8, 20, 3), .5f * color3(4, 4, 4)));
  addLight(scene, initLight(point3(4, 30, -1), .5f * color3(5, 5, 5)));

  mat.diffuseColor = color3(.2, 0.4, .3);
  mat.specularColor = color3(.2, 0.2, .2);
  mat.IOR = 2.382;
  mat.roughness = 0.005886;
  addObject(scene, initPlane(vec3(0, 1, 0), 0.2, mat));

  mat.diffuseColor = color3(.5, 0.09, .07);
  mat.specularColor = color3(.2, .2, .1);
  mat.IOR = 2.8382;
  mat.roughness = 0.00886;
  addObject(scene, initPlane(vec3(1, 0.0, 0.0), 2, mat));

  mat.diffuseColor = color3(0.1, 0.3, .05);
  mat.specularColor = color3(.5, .5, .5);
  mat.IOR = 2.9382;
  mat.roughness = 0.00886;
  addObject(scene, initPlane(vec3(0, 0, 1), 4, mat));

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
  //addLight(scene, initLight(point3(0, 1.7, 1), .5f * color3(3, 3, 3)));
  //addLight(scene, initLight(point3(3, 2, 3), .5f * color3(4, 4, 4)));
  //addLight(scene, initLight(point3(4, 3, -1), .5f * color3(5, 5, 5)));
  //addLight(scene, initLight(point3(1, 0, 1), .5f * color3(3, 3, 3)));

  vec3 v0 = vec3(1,0,0);
  vec3 v1 = vec3(0,1,0);
  vec3 v2 = vec3(0,0,-1);

  vec3 v1v0 = v1-v0;
  vec3 v2v0 = v2-v0;
  vec3 n = normalize(cross(v2v0, v1v0));

  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(0.7, 0.882, 0.786);
  mat.IOR = 6;
  mat.roughness = 0.0181;
  vec2 textures[3];
  addObject(scene, initTriangle(v0, v1, v2, n,textures,mat_lib[0]));

  v1 = vec3(1,0,0);
  v0 = vec3(0,1,0);
  v2 = vec3(0,0,1);

  v1v0 = v1-v0;
  v2v0 = v2-v0;
  n = normalize(cross(v2v0, v1v0));
  mat.diffuseColor = color3(0.5, 0,0);
  mat.specularColor = color3(0.5, 0, 0);
  mat.IOR = 1.01;
  mat.roughness = 2.2;
  mat.mtype = DIFFUSE;
  mat.m_texture = new image_texture("../assets/container2.png");
  textures[0] = vec2(0,0);
  textures[1] = vec2(1,0);
  textures[2] = vec2(0,1);
  addObject(scene, initTriangle(v0, v1, v2, n,textures,mat));
  v0 = vec3(1,-1, 1);
  textures[0] = vec2(1,1);
  textures[1] = vec2(1,0);
  textures[2] = vec2(0,1);
  addObject(scene, initTriangle(v0, v1, v2, n,textures,mat));
  
  mat.m_texture = nullptr;
  mat.diffuseColor = color3(0.016, 0.073, 0.04);
  mat.specularColor = color3(1.0, 1.056, 1.146);
  mat.IOR = 1.1481;
  mat.roughness = 0.0625;
  addObject(scene, initPlane(vec3(0, 1, 0), +1, mat));

  addObject(scene, initSphere(point3(-3, 1, 0), 2, mat));

  addLight(scene, initLight(point3(4, 2, 0), color3(1,1,1)));
  addLight(scene, initLight(point3(10, 100, 100), color3(1, 1, 1)));
  addLight(scene, initLight(point3(10, 10, -16), color3(1, 1, 1)));
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
        mat.IOR = materials[matIndex].ior;
        mat.roughness = materials[matIndex].ior;
        mat.diffuseColor.r = materials[matIndex].diffuse[0];
        mat.diffuseColor.g = materials[matIndex].diffuse[1];
        mat.diffuseColor.b = materials[matIndex].diffuse[2];
        //mat.m_texture = new image_texture(materials[matIndex].diffuse_texname.c_str());
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
      addObject(scene, initSmoothTriangle(v0, v1, v2, n, textures, normals[0], normals[1], normals[2],mat));
    }
  }
}

Scene *initScene6() {

  Scene *scene = initScene();
  setCamera(scene, point3(-4.5, 2.5, 3.5), vec3(0, 0, 0), vec3(0, 0.5, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.459f, 0.f, 0.878f));
  
  //addLight(scene, initLight(point3(100, 1, 1), color3(50, 50, 50)));
  //addLight(scene, initLight(point3(1, 0.1, 0.5), color3(5, 5, 5)));
  //addLight(scene, initLight(point3(0, 0, 0), color3(5, 5, 5)));
  //addLight(scene, initLight(point3(-6, 3, 7), color3(10, 10, 10)));

  //addLight(scene, initLight(point3(6, 3, -7), color3(1, 1, 1)));

  addLight(scene, initLight(point3(-6, 4, 0), color3(1, 1, 1)));
  addLight(scene, initLight(point3(4, 10, 10), color3(1, 1, 1)));

  addLight(scene, initLight(point3(-5, 4, 5), color3(1, 1, 1)));
  addLight(scene, initLight(point3(-5, 4, 4), color3(1, 1, 1)));

  //addLight(scene, initLight(point3(20, 1, -20), color3(1, 1, 1)));

  Material mat;
  mat.IOR = 1.3f;
  mat.diffuseColor = color3(.001f);
  mat.specularColor = color3(1.0f);
  mat.roughness = 1.f;

  addObjectsFromFile("../assets/new_bunny.obj", scene, mat_lib[0]);

  mat.diffuseColor = color3(0.7, 0, 0);
  mat.specularColor = color3(0.95, 0, 0);
  mat.IOR = 1.5;
  mat.roughness = 0.5;
  mat.mtype = DIFFUSE;
  addObject(scene, initSphere(point3(4, 0., -2), 0.4, mat));

  mat.diffuseColor = color3(0., 0, 0.7);
  mat.specularColor = color3(0., 0, 0.95);
  mat.IOR = 1.5;
  mat.roughness = 0.5;
  mat.mtype = DIFFUSE;
  addObject(scene, initSphere(point3(-0.5, 0., 0), 0.2, mat));

  mat.diffuseColor = color3(0.7, 0.5, 0);
  mat.specularColor = color3(0.95, 0.6, 0);
  mat.IOR = 1.5;
  mat.roughness = 0.5;
  mat.mtype = DIFFUSE;
  addObject(scene, initSphere(point3(-0.9, 0., 3), 0.2, mat));

  addObject(scene, initSphere(point3(4, 0., 3), 0.2, mat_lib[2]));

  mat.diffuseColor = color3(0.5);
  mat.specularColor = color3(0.5);
  mat.IOR = 1.5;
  mat.roughness = 0.0681;
  mat.mtype = DIFFUSE;
  mat.m_texture = new checker_texture(color3(1.0, 1.0, 1.0), color3(0, 0, 0));

  mat.diffuseColor = color3(0.5);
  mat.specularColor = color3(0.5);
  mat.IOR = 1.1;
  //mat.roughness = 0.0681;
  mat.roughness = 1.2;
  mat.mtype = DIFFUSE;
  //mat.m_texture = new checker_texture(color3(0.2, 0.3, 0.1), color3(0.9, 0.9, 0.9));
  mat.m_texture = new image_texture("../assets/chessboardtexture.png");
  addObject(scene, initPlane(vec3(0, 1, 0), 0.4, mat));
  mat.m_texture = nullptr;

  mat.diffuseColor = color3(0, 0.8, 0);
  mat.specularColor = color3(0.6, 0.6, 0);
  mat.IOR = 1.5;
  mat.roughness = 0.0681;
  //addObject(scene, initPlane(vec3(0, 0, 1), +10, mat));

  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(0.7, 0.882, 0.786);
  mat.IOR = 6;
  mat.roughness = 0.0181;
  addObject(scene, initSphere(point3(0, 0., -4), .3, mat));
  mat.m_texture = nullptr;

  mat.diffuseColor = color3(0.26, 0.036, 0.014);
  mat.specularColor = color3(1.0, 0.852, 1.172);
  mat.IOR = 1.3771;
  mat.roughness = 0.01589;
  addObject(scene, initSphere(point3(-1, -.05, 0), .15, mat));

  mat.diffuseColor = color3(0.014, 0.012, 0.012);
  mat.specularColor = color3(0.7, 0.882, 0.786);
  mat.IOR = 3;
  mat.roughness = 0.00181;
  addObject(scene, initSphere(point3(1, 0.05, 2), .25, mat));
  mat.m_texture = nullptr;

  mat.diffuseColor = color3(0.46, 0.136, 0.114);
  mat.specularColor = color3(0.8, 0.852, 0.8172);
  mat.IOR = 1.5771;
  mat.roughness = 0.01589;
  addObject(scene, initSphere(point3(1.6, 0., 2.6), 0.215, mat));

  mat.diffuseColor = color3(0.06, 0.26, 0.22);
  mat.specularColor = color3(0.70, 0.739, 0.721);
  mat.IOR = 1.3051;
  mat.roughness = 0.567;
  mat.m_texture = new image_texture("../assets/earthmap.png");
  addObject(scene, initSphere(point3(1.9, 0.05, 2.2), .25, mat));

  mat.diffuseColor = color3(0.012, 0.036, 0.406);
  mat.specularColor = color3(1.0, 0.965, 1.07);
  mat.IOR = 1.1153;
  mat.roughness = 0.068;
  mat.roughness = 0.18;
  addObject(scene, initSphere(point3(-2, 0, 1.5), .20, mat));
  
  return scene;
}

Scene *initScene7() {

  Scene *scene = initScene();
  setCamera(scene, point3(8, 5, 7), vec3(0, 0, 0), vec3(0, 1, 0), 60,
            float(WIDTH) / float(HEIGHT));
  setSkyColor(scene, color3(0.2, 0.8, 0.7)); 

  Material mat;
  mat.diffuseColor = color3(0.5);
  mat.specularColor = color3(0.5);
  mat.IOR = 1.5;
  mat.roughness = 0.0681;
  mat.mtype = DIFFUSE;
  addObjectsFromFile("../assets/werewolf.obj", scene, mat);
  
  //mat.diffuseColor = color3(.4, 0.8, .4);
  //mat.specularColor = color3(.4, 0.6, .2);
  //mat.IOR = 1.382;
  //mat.roughness = 0.05886;
  //addObject(scene, initPlane(vec3(100, 1, 0), +100, mat));
  
  //Material mat;
  //mat.diffuseColor = color3(0.5);
  //mat.specularColor = color3(0.5);
  //mat.IOR = 1.5;
  //mat.roughness = 0.0681;
  //mat.mtype = DIFFUSE;
  //mat.m_texture = new checker_texture(color3(0.2, 0.3, 0.1), color3(0.9, 0.9, 0.9));
  //addObject(scene, initPlane(vec3(0, 1, 0), 1, mat));

  //mat.m_texture = nullptr;
  //mat.diffuseColor = color3(1.0f, 0, 0);
  //mat.specularColor = color3(0.5f, 0.f, 0.f);
  //mat.roughness = 0.5f;
  //addObject(scene, initSphere(point3(-2, 0, -1), .25, mat));

  //mat.diffuseColor = color3(0, 1.0f, 0);
  //mat.specularColor = color3(0.f, 1.f, 0.f);
  //mat.roughness = 0.5f;
  //addObject(scene, initSphere(point3(0, 0, -1), .25, mat));

  //mat.diffuseColor = color3(0.5, 0 ,1.0f);
  //mat.specularColor = color3(0.5f, 0.f, 1.f);
  //mat.roughness = 0.5f;
  //addObject(scene, initSphere(point3(-1, 0, 0), .25, mat));

  //mat.diffuseColor = color3(0.5, 1.0f, 0.5);
  //mat.specularColor = color3(0.5f, 1.f, 0.5f);
  //mat.roughness = 0.5f;
  //addObject(scene, initSphere(point3(-2, 0, 0), .25, mat));

  addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initLight(point3(4, 10, -2), color3(1, 1, 1)));
  addLight(scene, initLight(point3(10, 10, 10), color3(1, 1, 1)));
  addLight(scene, initLight(point3(0, 40, 0), color3(1, 1, 1)));
  return scene;
}

Scene *initScene8() {
  Scene *scene = initScene();
  setCamera(scene, point3(28, 1, 28), vec3(0, 1, 0), vec3(0, 1, 0), 60,
            (float)WIDTH / (float)HEIGHT);
  setSkyColor(scene, color3(0.4, 0.9, 0.9));

  addLight(scene, initLight(point3(52, 10, 52), color3(1, 1, 1)));
  addLight(scene, initLight(point3(52, 10, 16), color3(1, 1, 1)));
  addLight(scene, initLight(point3(16, 10, 52), color3(1, 1, 1)));

  addLight(scene, initLight(point3(0, 50, 0), color3(1, 1, 1)));
  
  addObject(scene, initSphere(point3(24, 2, 24), 2.f, mat_lib[10]));
  addLight(scene, initLight(point3(5, 30, 5),color3(1, 1, 1)));

  Material mats;
  mats.diffuseColor = color3(0.f, 0.f, 0.5f);
  mats.specularColor = color3(0.f, 0.f, 0.7f);
  mats.roughness = 0.5f;

  addObject(scene, initSphere(point3(-5, 1.8, -5), 1.8f, mat_lib[8]));
  addObject(scene, initSphere(point3(18, 1.8, 26), 1.8f, mats));
  addObject(scene, initSphere(point3(26, 1.8, 18), 1.8f, mat_lib[8]));


  Material mat;
  mat.diffuseColor = color3(0.5);
  mat.specularColor = color3(0.5);
  mat.IOR = 1.1;
  //mat.roughness = 0.0681;
  mat.roughness = 1.2;
  mat.mtype = DIFFUSE;
  //mat.m_texture = new checker_texture(color3(0.2, 0.3, 0.1), color3(0.9, 0.9, 0.9));
  mat.m_texture = new image_texture("../assets/chessboardtexture.png");
  addObject(scene, initPlane(vec3(0, 1, 0), 0, mat));

  return scene;
}
