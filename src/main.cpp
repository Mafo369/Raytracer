
#include "Walnut/Application.h"
#include "Walnut/EntryPoint.h"

#include "Walnut/Image.h"
#include "Walnut/Timer.h"

#include "Renderer.h"
#include "Camera.h"

#include "example_scenes.h"

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

  default:
    scene = initScene0();
    break;
  }
  return scene;
}

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

bool g_ApplicationRunning = true;
#define WIDTH 800
#define HEIGHT 600

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 4) {
    printf("usage : %s [-rt] [filename] i\n", argv[0]);
    printf("        filename : where to save the result, whithout extention\n");
    printf("        i : scenen number, optional\n");
    printf("        -rt : real time rendering, optional, no filename needed\n");
    exit(0);
  }


  if(strcmp(argv[1], "-rt") == 0){
    while (g_ApplicationRunning)
    {
      Walnut::Application* app = Walnut::CreateApplication(argc, argv);
      app->Run();
      delete app;
    }
    return 0;
  }
  char basename[256];
  strncpy(basename, argv[1], 255);

  RenderImage *img = initImage(WIDTH, HEIGHT);
  int scene_id = 0;
  if (argc == 3) {
    scene_id = atoi(argv[2]);
  }
  Scene* scene = parseScene(scene_id);

  printf("render scene %d\n", scene_id);

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


