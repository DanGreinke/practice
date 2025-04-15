#include <stdio.h>
#include <math.h>

/* ------------------------------------------------------------------------------------
Calculate the Mandelbrot set, generate and save an image of it

Declare variables
- HEIGHT and WIDTH of image in pixels
- HEIGHT and WIDTH of region in complex plane (must be same aspect ratio as image)
- maximum number of iterations

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

Usage Note: PPM format is chosen for ease of use. I recommend exporting to a PNG and
then deleting the PPM file.
------------------------------------------------------------------------------------ */

#define HEIGHT 1440
#define WIDTH 2560
double y_min = -0.390625; //-1.125;
double y_max = -0.3375; //1.125;
double x_min = -0.75; //-2.0;
double x_max = -0.65625; //2.0;
int max_iterations = 1000;

// new x -> 1.316 (center), min = 1.253, max = 1.347
// new y -> 0.783
// 800 by 470
// 860 by 504
typedef struct {
    int arr[HEIGHT][WIDTH]; // 2D array to store iterations
} iteration_array;

int mandelbrot(double x_in, double y_in){
    double x_0 = x_in;
    double y_0 = y_in;
    double x = 0;
    double y = 0;
    int iterations = 0;
    double x_tmp;
    while (x*x + y*y < 4.0f && iterations < max_iterations){
        x_tmp = x*x - y*y + x_0;
        y = 2.0f*x*y + y_0;
        x = x_tmp;
        iterations++;
    }
    return iterations;
}

iteration_array generate_iteration_array(double x_min, double x_max, double y_min, double y_max){
    iteration_array out_array;
    int height = HEIGHT;
    int width = WIDTH;
    double x_step = (x_max - x_min) / width;
    double y_step = (y_max - y_min) / height;
    double x_in, y_in;
    // Loop over row of pixels
    for (int i = 0; i < height; i++){
        // Loop over pixel in row (column)
        for (int j = 0; j < width; j++){      
            x_in = x_min + j * x_step;
            y_in = y_min + i * y_step;
            out_array.arr[i][j] = mandelbrot(x_in, y_in);
        }
    }
    return out_array;
}

int color_map(int iterations){
    if (iterations == max_iterations){
        return 0; // Black
    }
    else{
        return iterations % 256; // Grayscale
    }
}

int main(void){
    iteration_array array = generate_iteration_array(x_min, x_max, y_min, y_max);
    FILE *img;
    int i, j;
    int height = HEIGHT;
    int width = WIDTH;

    // Initialize Image and Write Header
    img = fopen("mandelbrot.ppm", "w");
    // Writing Magic Number to the File 
    fprintf(img, "P2\n");  
    // Writing WIDTH and HEIGHT 
    fprintf(img, "%d %d\n", width, height);  
    // Writing the maximum gray value 
    fprintf(img, "255\n");
    // Loop over row in image
    for (i = 0; i < height; i++){
        // Loop over pixel in row (column)
        for (j = 0; j < width; j++){
            fprintf(img, "%d ", color_map(array.arr[i][j]));
        }
        fprintf(img, "\n");
    }
    fclose(img);
    return 0;
}