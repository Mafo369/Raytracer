#include "realtime/WalnutApp.h"
#include "rendering/CPURenderer.h"

bool g_ApplicationRunning = true;

int main( int argc, char* argv[] ) {
    // if ( argc < 2 || argc > 4 ) {
    //     printf( "usage : %s [-rt] [filename] i\n", argv[0] );
    //     printf( "        filename : where to save the result, whithout extention\n" );
    //     printf( "        i : scenen number, optional\n" );
    //     printf( "        -rt : real time rendering, optional, no filename needed\n" );
    //     exit( 0 );
    // }

#ifdef NDEBUG
    std::cout << "######### REALEASE MODE #########" << std::endl;
#else
    std::cout << "######### DEBUG MODE #########" << std::endl;
#endif

    // if ( strcmp( argv[1], "-rt" ) == 0 ) {
    //     while ( g_ApplicationRunning ) {
    //         Walnut::Application* app = Walnut::CreateApplication( argc, argv );
    //         app->Run();
    //         delete app;
    //     }
    //     return 0;
    // }

    std::string imageName = "test.png";

    // int scene_id = 0;
    // if ( argc == 3 ) { scene_id = atoi( argv[2] ); }
    Scene* scene = parseScene( 18 );

    auto renderer = new CPURenderer();
    renderer->init( imageName, { WIDTH, HEIGHT }, scene );
    renderer->run();
    renderer->destroy();

    freeScene( scene );
    scene = NULL;

    printf( "Done. Goodbye\n" );

    return 0;
}
