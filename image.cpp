#include "image.h"
#include <lodepng.h>
#include <stdio.h>
#include <string.h>

#define SAVE_PNG

color3 *getPixelPtr(Image *img, size_t x, size_t y) {
    return &(img->data[y * img->width + x]);
}

Image *initImage(size_t width, size_t height) {
    Image *img = (Image*) malloc(sizeof(Image));
    img->width = width;
    img->height = height;
    img->data = (color3 *)malloc(sizeof(color3)*width*height);
    return img;
}

void freeImage(Image *img) {
    free(img->data);
    free(img);
}

void saveImage(Image *img, char *basename) {
#ifdef SAVE_PNG
  char filename[256+4];
  strcpy(filename, basename);
  strcat(filename, ".png");

  size_t cpt = 0;
  unsigned char *image = new unsigned char [img->width*img->height*3];
  // write image to file
  for(unsigned y = 0; y < img->height; y++) {
    color3 *ptr = getPixelPtr(img, 0, img->height-y-1);
    for(unsigned x = 0; x < img->width; x++) {
      ivec3 c = clamp(ivec3(255.f**ptr), 0, 255);
      image[cpt++] = c.x;
      image[cpt++] = c.y;
      image[cpt++] = c.z; 
      ++ptr;
    }  
  }
  
  /*Encode the image*/
  unsigned error = lodepng_encode24_file(filename, image, img->width, img->height);
  delete [] image;
  /*if there's an error, display it*/
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
#else
    // save the image to basename.ppm
    {
        FILE *fp = NULL;
        char filename[256+4];
        strcpy(filename, basename);
        strcat(filename, ".ppm");

        fp = fopen(filename, "w");

        if (fp) {
            // file is created

            fprintf(fp, "P3\n");
            fprintf(fp, "%zu %zu\n255\n", img->width, img->height);

            // write image to file
            for(unsigned y = 0; y < img->height; y++) {
                color3 *ptr = getPixelPtr(img, 0, img->height-y-1);
                for(unsigned x = 0; x < img->width; x++) {
                    ivec3 c = clamp(ivec3(255.f**ptr), 0, 255);
                    unsigned char r = c.x;
                    unsigned char g = c.y;
                    unsigned char b = c.z;
                    fprintf(fp, "%d %d %d  ", r, g, b);
                    ++ptr;
                }
                fprintf(fp, "\n");
            }
            fclose(fp);
        }
    }
#endif

}
