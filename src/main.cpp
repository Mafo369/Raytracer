#include "realtime/WalnutApp.h"

bool g_ApplicationRunning = true;

int main( int argc, char* argv[] ) {
    if ( argc < 2 || argc > 4 ) {
        printf( "usage : %s [-rt] [filename] i\n", argv[0] );
        printf( "        filename : where to save the result, whithout extention\n" );
        printf( "        i : scenen number, optional\n" );
        printf( "        -rt : real time rendering, optional, no filename needed\n" );
        exit( 0 );
    }

#ifdef NDEBUG
    std::cout << "######### REALEASE MODE #########" << std::endl;
#else
    std::cout << "######### DEBUG MODE #########" << std::endl;
#endif

    if ( strcmp( argv[1], "-rt" ) == 0 ) {
        while ( g_ApplicationRunning ) {
            Walnut::Application* app = Walnut::CreateApplication( argc, argv );
            app->Run();
            delete app;
        }
        return 0;
    }
    char basename[256];
    strncpy( basename, argv[1], 255 );

    RenderImage* img = initImage( WIDTH, HEIGHT );
    int scene_id     = 0;
    if ( argc == 3 ) { scene_id = atoi( argv[2] ); }
    Scene* scene = parseScene( scene_id );

    printf( "render scene %d\n", scene_id );
    float posLight = -5;
    auto lightSize = 4.f;

    int i = 97;
    // for ( ; posLight < 5; posLight += 0.3 ) {
    // Transform modemMatrixE;
    // modemMatrixE.scale( lightSize, lightSize, lightSize );
    // modemMatrixE.translate( vec3( 0, 5, 0 ) );
    // auto light                = scene->objects[scene->objects.size() - 1];
    // light->geom.sphere.center = modemMatrixE.transformFrom( vec3( 0 ) );
    // light->geom.sphere.radius = lightSize;
    // light->transform          = modemMatrixE;
    // basename[4]               = char( i );
    printf( "save image to %s\n", basename );
    renderImage( img, scene );

    saveImage( img, basename );
    freeImage( img );
    img = initImage( WIDTH, HEIGHT );

    i++;
    //}
    freeImage( img );
    img = NULL;
    freeScene( scene );
    scene = NULL;

    printf( "done. Goodbye\n" );

    return 0;
}
