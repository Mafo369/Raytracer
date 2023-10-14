#include "example_scenes.h"
#include "raytracer.h"
#include "application.h"
#include <string.h>

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
        Application app;

        try {
            app.run(800, 600);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    }

    char basename[256];
    strncpy( basename, argv[1], 255 );

    RenderImage* img = initImage( WIDTH, HEIGHT );
    int scene_id     = 0;
    if ( argc == 3 ) { scene_id = atoi( argv[2] ); }
    Scene* scene = parseScene( scene_id );

    printf( "render scene %d\n", scene_id );
    printf( "save image to %s\n", basename );
    renderImage( img, scene );

    saveImage( img, basename );
    freeImage( img );
    img = initImage( WIDTH, HEIGHT );

    freeImage( img );
    img = NULL;
    freeScene( scene );
    scene = NULL;

    printf( "done. Goodbye\n" );

    return 0;
}
