#include "Renderer.h"

#include "Walnut/Random.h"

#include "raytracer.h"

#define WIDTH 200
#define HEIGHT 100

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
  setSkyColor(scene, color3(0.1f, 0.3f, 0.5f));
  Material mat;
  mat.IOR = 1.3;
  mat.roughness = 0.1;
  mat.specularColor = color3(0.5f);
  mat.m_texture = new image_texture("../assets/earthmap.png");

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

namespace Utils {

	static uint32_t ConvertToRGBA(const glm::vec4& color)
	{
		uint8_t r = (uint8_t)(color.r * 255.0f);
		uint8_t g = (uint8_t)(color.g * 255.0f);
		uint8_t b = (uint8_t)(color.b * 255.0f);
		uint8_t a = (uint8_t)(color.a * 255.0f);

		uint32_t result = (a << 24) | (b << 16) | (g << 8) | r;
		return result;
	}

}
Renderer::Renderer() {

  scene = initScene8();
  tree = initKdTree(scene);

  printf("End building tree\n");

  //! \todo initialize KdTree

}

void Renderer::OnResize(uint32_t width, uint32_t height)
{
	if (m_FinalImage)
	{
		// No resize necessary
		if (m_FinalImage->GetWidth() == width && m_FinalImage->GetHeight() == height)
			return;

    m_FinalImage->Resize(width, height);
    img = initImage(m_FinalImage->GetWidth(), m_FinalImage->GetHeight());
    setCamera(scene, point3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0), 60,
        (float)m_FinalImage->GetWidth() / (float)m_FinalImage->GetHeight());
  }
  else
  {
    m_FinalImage = std::make_shared<Walnut::Image>(width, height, Walnut::ImageFormat::RGBA);
    setCamera(scene, point3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0), 60,
        (float)m_FinalImage->GetWidth() / (float)m_FinalImage->GetHeight());
  }

	delete[] m_ImageData;
	m_ImageData = new uint32_t[width * height];
}

void Renderer::Render(const CameraI& camera)
{

	for (uint32_t y = 0; y < m_FinalImage->GetHeight(); y++)
	{
#pragma omp parallel for
		for (uint32_t x = 0; x < m_FinalImage->GetWidth(); x++)
		{
	    Ray ray;
			vec3 dir = camera.GetRayDirections()[x + y * m_FinalImage->GetWidth()];
      rayInit(&ray, camera.GetPosition(), normalize(dir));
      auto colorr = trace_ray(scene, &ray, tree);
      auto color = vec4(colorr.r, colorr.g, colorr.b, 1);
			color = glm::clamp(color, glm::vec4(0.0f), glm::vec4(1.0f));
			m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(color);
		}
	}

	m_FinalImage->SetData(m_ImageData);
}

glm::vec4 Renderer::TraceRay(const Ray& ray)
{
	float radius = 0.5f;
	// rayDirection = glm::normalize(rayDirection);

	// (bx^2 + by^2)t^2 + (2(axbx + ayby))t + (ax^2 + ay^2 - r^2) = 0
	// where
	// a = ray origin
	// b = ray direction
	// r = radius
	// t = hit distance

	float a = glm::dot(ray.dir, ray.dir);
	float b = 2.0f * glm::dot(ray.orig, ray.dir);
	float c = glm::dot(ray.orig, ray.orig) - radius * radius;

	// Quadratic forumula discriminant:
	// b^2 - 4ac

	float discriminant = b * b - 4.0f * a * c;
	if (discriminant < 0.0f)
		return glm::vec4(0, 0, 0, 1);

	// Quadratic formula:
	// (-b +- sqrt(discriminant)) / 2a

	float closestT = (-b - glm::sqrt(discriminant)) / (2.0f * a);
	float t0 = (-b + glm::sqrt(discriminant)) / (2.0f * a); // Second hit distance (currently unused)

	glm::vec3 hitPoint = ray.orig + ray.dir * closestT;
	glm::vec3 normal = glm::normalize(hitPoint);

	glm::vec3 lightDir = glm::normalize(glm::vec3(-1, -1, -1));
	float lightIntensity = glm::max(glm::dot(normal, -lightDir), 0.0f); // == cos(angle)

	glm::vec3 sphereColor(1, 0, 1);
	sphereColor *= lightIntensity;
	return glm::vec4(sphereColor, 1.0f);
}
