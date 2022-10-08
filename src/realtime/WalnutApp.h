#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"
#include "Walnut/Timer.h"

#include "Renderer.h"
#include "Camera.h"
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include "../example_scenes.h"
#include "../Object.h"

using namespace Walnut;

class ExampleLayer : public Walnut::Layer
{
public:
	ExampleLayer(Scene* scene)
		: m_Camera(45.0f, 0.1f, 100.0f) {
      m_Renderer.setScene(scene);
    }

	virtual void OnUpdate(float ts) override
	{
		m_Camera.OnUpdate(ts);
	}

	virtual void OnUIRender() override
	{
		ImGui::Begin("Settings");
		ImGui::Text("Last render: %.3fms", m_LastRenderTime);
		if (ImGui::Button("Render"))
		{
			Render();
		}
		ImGui::End();

    ImGui::Begin("Scene");
    const std::string cameraPos = glm::to_string(m_Camera.GetPosition());
    ImGui::Text(cameraPos.c_str());
    for(size_t i = 0; i < m_Renderer.scene->objects.size(); i++){
      ImGui::PushID(i);

      auto objMat = std::dynamic_pointer_cast<CookTorrance>(m_Renderer.scene->objects[i]->mat);
      ImGui::DragFloat("ior", &objMat->m_IOR, 0.1, 1.0, 6.0);
      ImGui::DragFloat("roughness", &objMat->m_roughness, 0.01, 0.0, 3.0);
      ImGui::ColorEdit3("dif", glm::value_ptr(objMat->m_diffuseColor));
      ImGui::ColorEdit3("spec", glm::value_ptr(objMat->m_specularColor));

      ImGui::Separator();
      ImGui::PopID();
    }
    ImGui::End();

		ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
		ImGui::Begin("Viewport");

		m_ViewportWidth = ImGui::GetContentRegionAvail().x;
		m_ViewportHeight = ImGui::GetContentRegionAvail().y;

		auto image = m_Renderer.GetFinalImage();
		if (image)
			ImGui::Image(image->GetDescriptorSet(), { (float)image->GetWidth(), (float)image->GetHeight() },
				ImVec2(0, 1), ImVec2(1, 0));

		ImGui::End();
		ImGui::PopStyleVar();

		Render();
	}

	void Render()
	{
		Timer timer;

		m_Renderer.OnResize(m_ViewportWidth, m_ViewportHeight);
		m_Camera.OnResize(m_ViewportWidth, m_ViewportHeight);
		m_Renderer.Render(m_Camera);

		m_LastRenderTime = timer.ElapsedMillis();
	}
private:
	Renderer m_Renderer;
	CameraI m_Camera;
	uint32_t m_ViewportWidth = 0, m_ViewportHeight = 0;

	float m_LastRenderTime = 0.0f;
};


Walnut::Application* Walnut::CreateApplication(int argc, char** argv)
{
	Walnut::ApplicationSpecification spec;
	spec.Name = "Ray Tracing";

	Walnut::Application* app = new Walnut::Application(spec);

  int scene_id = 0;
  if (argc == 3) {
    scene_id = atoi(argv[2]);
  }else if(argc == 4){
    scene_id = atoi(argv[3]);
  }
  Scene* scene = parseScene(scene_id);
  auto layer = std::make_shared<ExampleLayer>(scene);
  app->PushLayer(layer);
	app->SetMenubarCallback([app]()
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Exit"))
			{
				app->Close();
			}
			ImGui::EndMenu();
		}
	});
	return app;
}
