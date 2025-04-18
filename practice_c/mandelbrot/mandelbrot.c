#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> // For strcmp
#include <getopt.h> // For getopt_long
#include <unistd.h> // For optind, opterr, optopt (often included by getopt.h)
#include <stdbool.h> // For bool type


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
- Parse command line arguments
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

// Structure to hold the arguments
struct arguments {
    int width;
    int height;
    int zoom;
    double x_center;
    double y_center;
    int max_iterations;
    char *output_file;
};

// Structure for area selection
typedef struct {
    double x_min;
    double x_max;
    double y_min;
    double y_max;
} region;

// --- Function Prototypes ---
static void print_usage(const char *prog_name);
region generate_region(int zoom, double x_c, double y_c, int width, int height);
int mandelbrot(double x_in, double y_in, int max_iterations);
int* generate_iteration_array(region r, int width, int height, int max_iterations);
int color_map(int iterations, int max_iterations);
// --- End Function Prototypes ---

// Function to print usage instructions
static void print_usage(const char *prog_name) {
    fprintf(stderr, "Usage: %s [options]\n\n", prog_name);
    fprintf(stderr, "Generate a PPM image of the Mandelbrot set.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -W, --width=PIXELS    Image width in pixels (default: 1280)\n");
    fprintf(stderr, "  -H, --height=PIXELS   Image height in pixels (default: 720)\n");
    fprintf(stderr, "  -z, --zoom=LEVEL      Zoom level (integer, default: 1)\n");
    fprintf(stderr, "  -x, --xcenter=COORD   X-coordinate of the center (double, default: 0.0)\n");
    fprintf(stderr, "  -y, --ycenter=COORD   Y-coordinate of the center (double, default: 0.0)\n");
    fprintf(stderr, "  -m, --maxiter=COUNT   Maximum iterations (integer, default: 1000)\n");
    fprintf(stderr, "  -o, --output=FILE     Output PPM file name (default: mandelbrot.ppm)\n");
    fprintf(stderr, "  -h, --help            Display this help message and exit\n");
    fprintf(stderr, "  --version             Display version information and exit\n");
    fprintf(stderr, "\nExample:\n");
    fprintf(stderr, "  %s -W 1920 -H 1080 -z 5000 -x -0.745 -y 0.11 -m 1500 -o my_mandelbrot.ppm\n", prog_name);
}

region generate_region(int zoom, double x_c, double y_c, int width, int height) {
    region r;

    // Ensure zoom is positive
    if (zoom <=0) {
        fprintf(stderr, "Warning: Zoom level must be positive. Using default zoom 1.\n");
        zoom = 1;
    }
    // Ensure dimensions are positive
    if (width <= 0 || height <= 0) {
        fprintf(stderr, "Error: Image dimensions must be positive. Cannot generate image.\n");
        exit(1);
    }

    // Calculate bounds of region, based on aspect ratio
    double aspect_ratio = (double)width/(double)height;
    double base_region_height=2.5; // y from -1.25 to 1.25 at zoom 1
    double region_height = base_region_height / zoom;
    double region_width = region_height*aspect_ratio;

    // Calculate corners of zoom region
    r.x_min = x_c - (region_width / 2.0);   
    r.x_max = x_c + (region_width / 2.0);
    r.y_min = y_c - (region_height / 2.0);
    r.y_max = y_c + (region_height / 2.0);

    return r;
}

int mandelbrot(double x_in, double y_in, int max_iterations) {
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

int* generate_iteration_array(region r, int width, int height, int max_iterations) {
    // check for valid dimensions passed
    if (width <= 0 || height <= 0) {
        fprintf(stderr, "Error: Image dimensions must be positive. Cannot generate image.\n");
        return NULL;
    }

    // check if region is valid
    if (r.x_min >= r.x_max || r.y_min >= r.y_max) {
        fprintf(stderr, "Error: Invalid region coordinates. Cannot generate image.\n");
        return NULL;
    }

    // check if max_iterations is valid
    if (max_iterations <= 0) {
        fprintf(stderr, "Error: Maximum iterations must be positive. Cannot generate image.\n");
        return NULL;
    }
    
    // Calculate step size to map pixels to complex plane
    double x_step = (r.x_max - r.x_min) / width;
    double y_step = (r.y_max - r.y_min) / height;
    double x_in, y_in;

    // Allocate memory on the heap
    size_t num_elements = (size_t)height * width;

    // Check for multiplication overflow before allocation
    if (height > 0 && num_elements / height != (size_t)width) {
         perror("Integer overflow calculating allocation size");
         return NULL;
    }
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
            out_array_data[i*width + j] = mandelbrot(x_in, y_in, max_iterations);
        }
    }
    return out_array_data;
}

int color_map(int iterations, int max_iterations) {
    if (iterations == max_iterations) {
        return 0; // Black for points in the Mandelbrot set
    }
    else{
        return iterations % 256; // Grayscale
    }
}

int main(int argc, char *argv[]) {
    // Default values for arguments
    struct arguments arguments;
    arguments.width = 1280;
    arguments.height = 720;
    arguments.zoom = 1;
    arguments.x_center = 0.0;
    arguments.y_center = 0.0;
    arguments.max_iterations = 1000;
    arguments.output_file = "mandelbrot.ppm"; // Default output file

    // --- getopt_long setup ---
    int opt;
    char *endptr; // For strtol/strtod error checking

    // Define long options
    // The last element must be all zeros.
    static struct option long_options[] = {
        {"width",   required_argument, 0, 'W'},
        {"height",  required_argument, 0, 'H'},
        {"zoom",    required_argument, 0, 'z'},
        {"xcenter", required_argument, 0, 'x'},
        {"ycenter", required_argument, 0, 'y'},
        {"maxiter", required_argument, 0, 'm'},
        {"output",  required_argument, 0, 'o'},
        {"help",    no_argument,       0, 'h'},
        {"version", no_argument,       0, 'v'}, // Using 'v' for version
        {0, 0, 0, 0} // Terminator
    };

    // A string listing valid short options letters.
    // Options requiring an argument are followed by a colon (:)
    const char *short_options = "W:H:z:x:y:m:o:hv";

    // Suppress getopt's default error messages
    // opterr = 0;

    // Parsing loop
    while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
        switch (opt) {
            case 'W':
                arguments.width = strtol(optarg, &endptr, 10);
                if (*endptr != '\0' || arguments.width <= 0) {
                    fprintf(stderr, "Error: Invalid width '%s'. Must be a positive integer.\n", optarg);
                    print_usage(argv[0]);
                    return 1;
                }
                break;
            case 'H':
                arguments.height = strtol(optarg, &endptr, 10);
                if (*endptr != '\0' || arguments.height <= 0) {
                    fprintf(stderr, "Error: Invalid height '%s'. Must be a positive integer.\n", optarg);
                    print_usage(argv[0]);
                    return 1;
                }
                break;
            case 'z':
                arguments.zoom = strtol(optarg, &endptr, 10);
                if (*endptr != '\0' || arguments.zoom <= 0) {
                    fprintf(stderr, "Error: Invalid zoom level '%s'. Must be a positive integer.\n", optarg);
                    print_usage(argv[0]);
                    return 1;
                }
                break;
            case 'x':
                arguments.x_center = strtod(optarg, &endptr);
                if (*endptr != '\0') {
                    fprintf(stderr, "Error: Invalid x-center '%s'. Must be a floating-point number.\n", optarg);
                    print_usage(argv[0]);
                    return 1;
                }
                break;
            case 'y':
                arguments.y_center = strtod(optarg, &endptr);
                if (*endptr != '\0') {
                    fprintf(stderr, "Error: Invalid y-center '%s'. Must be a floating-point number.\n", optarg);
                    print_usage(argv[0]);
                    return 1;
                }
                break;
            case 'm':
                arguments.max_iterations = strtol(optarg, &endptr, 10);
                if (*endptr != '\0' || arguments.max_iterations <= 0) {
                    fprintf(stderr, "Error: Invalid max iterations '%s'. Must be a positive integer.\n", optarg);
                    print_usage(argv[0]);
                    return 1;
                }
                break;
            case 'o':
                if (optarg == NULL || strcmp(optarg, "") == 0) {
                    fprintf(stderr, "Error: Output filename cannot be empty.\n");
                    print_usage(argv[0]);
                    return 1;
                }
                arguments.output_file = optarg;
                break;
            case 'h': // Help
                print_usage(argv[0]);
                return 0; // Successful exit after printing help
            case 'v': // Version
                // Use the version string defined earlier if available, or hardcode
                printf("mandelbrot 1.1 (getopt_long version)\n");
                // Add copyright/license info if desired
                return 0; // Successful exit after printing version
            case '?': // Unknown option or missing argument
                // getopt_long already prints an error message by default if opterr != 0
                fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
                return 1;
            default:
                // Should not happen with the options defined
                fprintf(stderr, "Internal error: Unhandled option character code 0%o\n", opt);
                print_usage(argv[0]);
                return 1;
        }
    }

    // Check for any non-option arguments left (we don't expect any)
    if (optind < argc) {
        fprintf(stderr, "Error: Unexpected non-option arguments: ");
        while (optind < argc) {
            fprintf(stderr, "%s ", argv[optind++]);
        }
        fprintf(stderr, "\n");
        print_usage(argv[0]);
        return 1;
    }

    // --- End getopt_long setup ---

    // Print the settings being used
    printf("Generating Mandelbrot set with:\n");
    printf("  Resolution: %d x %d\n", arguments.width, arguments.height);
    printf("  Center:     (%.10f, %.10f)\n", arguments.x_center, arguments.y_center);
    printf("  Zoom:       %d\n", arguments.zoom);
    printf("  Max Iter:   %d\n", arguments.max_iterations);
    printf("  Output:     %s\n", arguments.output_file);

    // Set the region of the complex plane to be visualized
    region r = generate_region(arguments.zoom, arguments.x_center, arguments.y_center,
                               arguments.width, arguments.height);

    printf("  Region:     x[%.10f, %.10f], y[%.10f, %.10f]\n", r.x_min, r.x_max, r.y_min, r.y_max);

    // Receive pointer from generate_iteration_array function
    int *iteration_data = generate_iteration_array(r, arguments.width, arguments.height,
                                                  arguments.max_iterations);

    if (iteration_data == NULL) {
        return 1; // Error already be printed by generate_iteration_array, or malloc.
    }

    FILE *img;
    int i, j;
    int height = arguments.height;
    int width = arguments.width;

    // Initialize Image and Write Header
    img = fopen(arguments.output_file, "w");
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
            // Pass iterations to color map
            fprintf(img, "%d ", color_map(iteration_data[i*width + j], arguments.max_iterations));
        }
        // Newline after each row
        fprintf(img, "\n");
    }
    // Close the file
    if (fclose(img) != 0) {
        perror("Error closing output file");
        // Continue to free memory, but report the error
    } else {
       printf("Successfully wrote image to %s\n", arguments.output_file);
    }

    // Free the allocated memory
    free(iteration_data); iteration_data = NULL;
    return 0; // Success
}