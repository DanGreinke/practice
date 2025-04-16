#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* ------------------------------------------------------------------------------------
Calculate the Mandelbrot set, generate and save an image of it

Declare variables
- HEIGHT and WIDTH of image in pixels
- HEIGHT and WIDTH of region in complex plane (must be same aspect ratio as image)
- maximum number of iterations
- Create region structure to assist with zooming in

Generate region of complex plane
- Take zoom, x_center, y_center as input
- Calculate x_min, x_max, y_min, y_max based on zoom and center
  - Lock the aspect ratio of the image to 16:9

Mandelbrot function
- Take x,y as input, iterate the Mandelbrot function
- ouput number of iterations before escape, or max iterations if it never escapes

Generate a 2D array in the complex plane
- for each pixel in the image, map it to a point in the complex plane
- for each point in the complex plane, iterate the Mandelbrot function
- Save the number of iterations before excape, or max iterations if it never escapes

Color Map
- Use color map to determine color of pixel based on number of iterations
- Use a simple color map (e.g. grayscale) for now

Main function
- Call 2D array generation function
  - Call Mandelbrot function to generate number of iterations
- Call image generation function
- Save image to file

Usage Notes: 

PPM format is chosen for ease of use and to minimize dependencies. I recommend exporting 
to a PNG and then deleting the PPM file after running this script.

Common 16:9 resolutions to try
- 1280x720
- 1920x1080
- 2560x1440
- 4480x2520
------------------------------------------------------------------------------------ */

#define HEIGHT 720
#define WIDTH 1280
int zoom = 120;
double x_center = -0.695;
double y_center = -0.380;
int max_iterations = 1000;

typedef struct {
    double x_min;
    double x_max;
    double y_min;
    double y_max;
} region;

region generate_region(int zoom, double x_c, double y_c) {
    region r;
    r.x_min = x_c - 2.0/zoom;   
    r.x_max = x_c + 2.0/zoom;
    r.y_min = y_c - 1.125/zoom;
    r.y_max = y_c + 1.125/zoom;
    return r;
}

int mandelbrot(double x_in, double y_in) {
    double x_0 = x_in;
    double y_0 = y_in;
    double x = 0;
    double y = 0;
    int iterations = 0;
    double x_tmp;
    while (x*x + y*y < 4.0f && iterations < max_iterations) {
        x_tmp = x*x - y*y + x_0;
        y = 2.0f*x*y + y_0;
        x = x_tmp;
        iterations++;
    }
    return iterations;
}

int* generate_iteration_array(region r) {
    int height = HEIGHT;
    int width = WIDTH;
    double x_step = (r.x_max - r.x_min) / width;
    double y_step = (r.y_max - r.y_min) / height;
    double x_in, y_in;
    // Allocate memory on the heap
    int *out_array_data = malloc(height * width * sizeof(int));

    // Check if memory allocation was successful
    if (out_array_data == NULL) {
        perror("Memory allocation failed");
        return NULL;
    }

    // Loop over row of pixels
    for (int i = 0; i < height; i++) {
        // Loop over pixel in row (column)
        for (int j = 0; j < width; j++) {      
            x_in = r.x_min + j * x_step;
            y_in = r.y_min + i * y_step;
            out_array_data[i*width + j] = mandelbrot(x_in, y_in);
        }
    }
    return out_array_data;
}

int color_map(int iterations) {
    if (iterations == max_iterations) {
        return 0; // Black
    }
    else{
        return iterations % 256; // Grayscale
    }
}

int main(void) {
    // Set the region of the complex plane to be visualized
    region r = generate_region(zoom, x_center, y_center);
    
    // Receive pointer from generate_iteration_array function
    int *iteration_data = generate_iteration_array(r);

    FILE *img;
    int i, j;
    int height = HEIGHT;
    int width = WIDTH;

    // Initialize Image and Write Header
    img = fopen("mandelbrot.ppm", "w");
    if (img == NULL) {
        perror("Error opening file");
        free(iteration_data);
        return 1;
    }
    // Writing Magic Number, dimensions, & max color value
    fprintf(img, "P2\n");  
    fprintf(img, "%d %d\n", width, height);  
    fprintf(img, "255\n");

    // Loop over row in image
    for (i = 0; i < height; i++) {
        // Loop over pixel in row (column)
        for (j = 0; j < width; j++) {
            fprintf(img, "%d ", color_map(iteration_data[i*width + j]));
        }
        fprintf(img, "\n");
    }
    fclose(img);
    // Free the allocated memory
    free(iteration_data); iteration_data = NULL;
    return 0;
}