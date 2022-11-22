#include "image.h"
#include <iostream>
#include <lodepng.h>
#include <stdio.h>
#include <string.h>

#define SAVE_PNG

RenderImage* loadPng( const char* filename ) { // Not used -> tried to implement textures
    unsigned int w;
    unsigned int h;
    std::vector<unsigned char> image;
    unsigned int error = lodepng::decode( image, w, h, filename );

    if ( error ) {
        printf( "Unable to load texture: %s\n", filename );
        return nullptr;
    }

    RenderImage* img = initImage( w, h );

    size_t idx = 0;
    for ( size_t j = 0; j < img->height; j++ ) {
        for ( size_t i = 0; i < img->width; i++ ) {
            color3* ptr    = getPixelPtr( img, i, j );
            unsigned int r = image[idx++];
            unsigned int g = image[idx++];
            unsigned int b = image[idx++];
            color3 color   = color3( r, g, b ) / 255.f;
            *ptr           = color;
            idx++;
        }
    }

    return img;
}

color3* getPixelPtr( RenderImage* img, size_t x, size_t y ) {
    return &( img->data[y * img->width + x] );
}

RenderImage* initImage( size_t width, size_t height ) {
    RenderImage* img = (RenderImage*)malloc( sizeof( RenderImage ) );
    img->width       = width;
    img->height      = height;
    img->data        = (color3*)malloc( sizeof( color3 ) * width * height );
    return img;
}

void freeImage( RenderImage* img ) {
    free( img->data );
    free( img );
}

void saveImage( RenderImage* img, char* basename ) {
#ifdef SAVE_PNG
    char filename[256 + 4];
    strcpy( filename, basename );
    strcat( filename, ".png" );

    size_t cpt           = 0;
    unsigned char* image = new unsigned char[img->width * img->height * 3];
    // write image to file
    for ( unsigned y = 0; y < img->height; y++ ) {
        color3* ptr = getPixelPtr( img, 0, y );
        for ( unsigned x = 0; x < img->width; x++ ) {
            ivec3 c      = clamp( ivec3( 255.f * *ptr ), 0, 255 );
            image[cpt++] = c.x;
            image[cpt++] = c.y;
            image[cpt++] = c.z;
            ++ptr;
        }
    }

    /*Encode the image*/
    unsigned error = lodepng_encode24_file( filename, image, img->width, img->height );
    delete[] image;
    /*if there's an error, display it*/
    if ( error ) printf( "error %u: %s\n", error, lodepng_error_text( error ) );
#else
    // save the image to basename.ppm
    {
        FILE* fp = NULL;
        char filename[256 + 4];
        strcpy( filename, basename );
        strcat( filename, ".ppm" );

        fp = fopen( filename, "w" );

        if ( fp ) {
            // file is created

            fprintf( fp, "P3\n" );
            fprintf( fp, "%zu %zu\n255\n", img->width, img->height );

            // write image to file
            for ( unsigned y = 0; y < img->height; y++ ) {
                color3* ptr = getPixelPtr( img, 0, img->height - y - 1 );
                for ( unsigned x = 0; x < img->width; x++ ) {
                    ivec3 c         = clamp( ivec3( 255.f * *ptr ), 0, 255 );
                    unsigned char r = c.x;
                    unsigned char g = c.y;
                    unsigned char b = c.z;
                    fprintf( fp, "%d %d %d  ", r, g, b );
                    ++ptr;
                }
                fprintf( fp, "\n" );
            }
            fclose( fp );
        }
    }
#endif
}
