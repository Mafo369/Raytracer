#include "realtime/WalnutApp.h"

bool g_ApplicationRunning = true;

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


