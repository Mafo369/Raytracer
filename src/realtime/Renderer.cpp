#include "Renderer.h"

#include "Walnut/Random.h"

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

Renderer::Renderer() {}

void Renderer::OnResize(uint32_t width, uint32_t height)
{
	if (m_FinalImage)
	{
		// No resize necessary
		if (m_FinalImage->GetWidth() == width && m_FinalImage->GetHeight() == height)
			return;

    m_FinalImage->Resize(width, height);
    img = initImage(m_FinalImage->GetWidth(), m_FinalImage->GetHeight());
    setCameraFOV(scene, point3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0), 60,
        (float)m_FinalImage->GetWidth(), (float)m_FinalImage->GetHeight(), 0.01);
  }
  else
  {
    m_FinalImage = std::make_shared<Walnut::Image>(width, height, Walnut::ImageFormat::RGBA);
    setCameraFOV(scene, point3(0, 0, 0), vec3(0, 0, 0), vec3(0, 0, 0), 60,
        (float)m_FinalImage->GetWidth(), (float)m_FinalImage->GetHeight(), 0.01);
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
      rayInit(&ray, camera.GetPosition(), normalize(dir), vec2(x,y));
      Intersection intersection;
      auto colorr = trace_ray(scene, &ray, tree, &intersection);
      auto color = vec4(colorr.r, colorr.g, colorr.b, 1);
			color = glm::clamp(color, glm::vec4(0.0f), glm::vec4(1.0f));
			m_ImageData[x + y * m_FinalImage->GetWidth()] = Utils::ConvertToRGBA(color);
		}
	}

	m_FinalImage->SetData(m_ImageData);
}
