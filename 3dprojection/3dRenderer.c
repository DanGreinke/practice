#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
// Credit to Gemini for helping me vibecode this based on 3dProjPerspective.py
// This code is a simple 3D renderer using SDL2 in C.

// Simple 3D vector structure
// x, y, z coordinates, w for homogeneous coordinates
typedef struct {
    float x;
    float y;
    float z;
    float w;
} vec3d;

// Triangle structure
typedef struct {
    vec3d p[3]; // Three vertices of the triangle
    vec3d c;    // Color of the triangle
} triangle;

// Mesh structure
typedef struct {
    triangle* tris; // Array of triangles
    int num_tris;   // Number of triangles in the mesh
} mesh;

// Matrix 4x4 structure
typedef struct {
    float m[4][4]; // 4x4 matrix
} mat4x4;

typedef struct {
    void* data;     // Pointer to the dynamically allocated array (now void*)
    size_t element_size; // Size of each element
    int size;      // Current number of elements
    int capacity;  // Maximum number of elements
} DynamicArray;

// Key State Array
const int NUM_KEYS = 512; // Assuming you need to track up to 512 keys
bool key_states[NUM_KEYS] = {false};

// Function to get key state
bool GetKey(int key) {
    if (key >= 0 && key < NUM_KEYS) {
        return key_states[key];
    }
    return false; // Key out of range
}

// Function Prototypes
vec3d vec3d_mult_by_mat4x4(vec3d vec, mat4x4 m);
mat4x4 mat4x4_mult_by_mat4x4(mat4x4 m_1, mat4x4 m_2);
vec3d vec3d_add(vec3d a, vec3d b);
vec3d vec3d_sub(vec3d a, vec3d b);
vec3d vec3d_mult(vec3d a, float k);
vec3d vec3d_div(vec3d a, float k);
float vec3d_dot(vec3d a, vec3d b);
vec3d vec3d_cross(vec3d a, vec3d b);
float vec3d_length(vec3d v);
mat4x4 matrix_point_at(vec3d pos, vec3d target, vec3d up);
vec3d vec3d_normalize(vec3d v);
triangle tri_div_by_w(triangle tri);
triangle project_triangle(triangle tri);
triangle transform_triangle(triangle tri, mat4x4 world_mat);
void make_translation_matrix(float x_offset, float y_offset, float z_offset, mat4x4* translation_matrix);
void make_rotation_matrix(float theta_x, float theta_y, float theta_z, mat4x4* rotation_matrix);
float adjust(float p, float scale, float offset);
float tri_depth(triangle tri);
vec3d get_surface_normal(triangle tri);
bool check_visibility(triangle tri, vec3d cam, vec3d* tri_normal);
float illuminate_triangle(vec3d normalized_light_vec, vec3d tri_normal);
mesh load_obj(const char* obj_file);
void draw_polygon(SDL_Renderer* renderer, int x1, int y1, int x2, int y2, int x3, int y3, int r, int g, int b);
int compare_triangles(const void* a, const void* b);
bool GetKey(int key); // For checking keyboard input state
mat4x4 Matrix_QuickInverse(mat4x4 m); // For inverting the camera matrix
void draw_triangle_edges(SDL_Renderer* renderer, float x1, float y1, float x2, float y2, float x3, float y3, int r, int g, int b); // For drawing wireframes
vec3d vector_intersect_plane(vec3d plane_p, vec3d plane_n, vec3d line_start, vec3d line_end); // Used in clipping
float dist(vec3d p, vec3d plane_p, vec3d plane_n); // Used in clipping
int triangle_clip_against_plane(vec3d plane_p, vec3d plane_n, triangle *in_tri, triangle *out_tri1, triangle *out_tri2); // Clipping function
void dynamicArrayInit(DynamicArray* arr, size_t element_size, int initialCapacity); // Dynamic array helper
void dynamicArrayPushBack(DynamicArray* arr, const void* value); // Dynamic array helper
void* dynamicArrayPopFront(DynamicArray* arr); // Dynamic array helper
void dynamicArrayFree(DynamicArray* arr); // Dynamic array helper


// Global variables (similar to Python script)
const float f_near = 0.1f;          // Near clipping plane
const float f_far = 1000.0f;        // Far clipping plane
const float f_fov = 90.0f;          // Field of view in degrees
const float f_aspect_ratio = 1.0f;  // Aspect ratio (width/height)
float f_fov_rad;                    // Field of view in radians
vec3d viewer_Camera = {0.0f, 0.0f, 0.0f, 1.0f}; // Camera position
vec3d v_look_dir;
vec3d light_direction = {0.0f, 0.0f, -1.0f, 1.0f}; // Direction of the light source
vec3d light_color = {255.0f, 255.0f, 255.0f, 1.0f}; // Color of the light source
const float scale = 200.0f;         // Scaling factor for rendering
float angle_x = 0.0f, angle_y = 0.0f, angle_z = 0.0f; // Rotation angles
const int WINDOW_SIZE = 800;        // Size of the window
const float TICK = 60.0f;           // Target frames per second
const float ROTATE_SPEED = 1.0f / TICK; // Rotation speed per frame
const float obj_distance = 8.0f;    // Distance of the object from the camera

// Function to multiply a vector by a 4x4 matrix
vec3d vec3d_mult_by_mat4x4(vec3d vec, mat4x4 m) {
    // Convert the vector to a 4-element array for matrix multiplication
    float tmp_vec[4] = {vec.x, vec.y, vec.z, vec.w};
    float result[4] = {0};
    // Perform matrix multiplication
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i] += tmp_vec[j] * m.m[j][i];
        }
    }
    return (vec3d){result[0], result[1], result[2], result[3]};
}

mat4x4 mat4x4_mult_by_mat4x4(mat4x4 m_1, mat4x4 m_2) {
    mat4x4 result;
    // Perform matrix multiplication
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.m[i][j] = 0;
            for (int k = 0; k < 4; k++) {
                result.m[i][j] += m_1.m[i][k] * m_2.m[k][j];
            }
        }
    }
    return result;
}

// Vector Operations
vec3d vec3d_add(vec3d a, vec3d b) {
    return (vec3d){a.x + b.x, a.y + b.y, a.z + b.z, 1.0f};
}

vec3d vec3d_sub(vec3d a, vec3d b) {
    return (vec3d){a.x - b.x, a.y - b.y, a.z - b.z, 1.0f};
}

vec3d vec3d_mult(vec3d a, float k) {
    return (vec3d){a.x * k, a.y * k, a.z * k, 1.0f};
}

vec3d vec3d_div(vec3d a, float k) {
    return (vec3d){a.x / k, a.y / k, a.z / k, 1.0f};
}

float vec3d_dot(vec3d a, vec3d b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3d vec3d_cross(vec3d a, vec3d b) {
    return (vec3d){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
        1.0f
    };
}

float vec3d_length(vec3d v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

mat4x4 matrix_point_at(vec3d pos, vec3d target, vec3d up) {
    // Calculate the forward vector
    vec3d new_forward = vec3d_sub(target, pos);
    new_forward = vec3d_normalize(new_forward);
    // Calculate the up direction
    vec3d a = vec3d_mult(new_forward, vec3d_dot(up, new_forward));
    vec3d new_up = vec3d_sub(up, a);
    new_up = vec3d_normalize(new_up);
    // Calculate the right direction
    vec3d new_right = vec3d_cross(new_up, new_forward);
    // Construct Dimensioning & Translation Matrix
    mat4x4 matrix;
    matrix.m[0][0] = new_right.x;   matrix.m[0][1] = new_right.y;   matrix.m[0][2] = new_right.z;   matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = new_up.x;      matrix.m[1][1] = new_up.y;      matrix.m[1][2] = new_up.z;      matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = new_forward.x; matrix.m[2][1] = new_forward.y; matrix.m[2][2] = new_forward.z; matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = pos.x;         matrix.m[3][1] = pos.y;         matrix.m[3][2] = pos.z;         matrix.m[3][3] = 1.0f;
    return matrix;
}

mat4x4 Matrix_QuickInverse(mat4x4 m) 
// Only for Rotation/Translation Matrices
{
    mat4x4 matrix;
    matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
    matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
    matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
    matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
    matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
    matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
    matrix.m[3][3] = 1.0f;
    return matrix;
}

// Function to normalize a vector
vec3d vec3d_normalize(vec3d v) {
    // Calculate the magnitude (length) of the vector
    float l = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    // Divide each component by the magnitude to get a unit vector
    return (vec3d){v.x / l, v.y / l, v.z / l, v.w};
}

// Function to divide a triangle by w (homogeneous divide)
triangle tri_div_by_w(triangle tri) {
    triangle out_tri;
    // Divide each vertex by its w component

    for (int i = 0; i < 3; i++) {
        // Return input triangle if any component is zero
        if (tri.p[i].w == 0) {
            return tri;
        }
        // Perform the division
        out_tri.p[i] = (vec3d){tri.p[i].x / tri.p[i].w, tri.p[i].y / tri.p[i].w, tri.p[i].z / tri.p[i].w, tri.p[i].w};
    }
    return out_tri;
}

// Function to project a triangle from 3D to 2D
triangle project_triangle(triangle tri) {
    // Define the projection matrix
    mat4x4 projection_matrix;
    projection_matrix.m[0][0] = -f_aspect_ratio * f_fov_rad; projection_matrix.m[0][1] = 0; projection_matrix.m[0][2] = 0; projection_matrix.m[0][3] = 0; // inverting x axis
    projection_matrix.m[1][0] = 0; projection_matrix.m[1][1] = -f_fov_rad; projection_matrix.m[1][2] = 0; projection_matrix.m[1][3] = 0; // inverting y axis
    projection_matrix.m[2][0] = 0; projection_matrix.m[2][1] = 0; projection_matrix.m[2][2] = f_far / (f_far - f_near); projection_matrix.m[2][3] = 1;
    projection_matrix.m[3][0] = 0; projection_matrix.m[3][1] = 0; projection_matrix.m[3][2] = (-f_far * f_near) / (f_far - f_near); projection_matrix.m[3][3] = 0;
    triangle projected_tri;
    // Project each vertex of the triangle
    for (int i = 0; i < 3; i++) {
        projected_tri.p[i] = vec3d_mult_by_mat4x4(tri.p[i], projection_matrix);
    }
    return projected_tri;
}

// Function to transform a triangle by a world matrix
triangle transform_triangle(triangle tri, mat4x4 world_mat) {
    triangle transformed_tri;
    // Transform each vertex of the triangle
    for (int i = 0; i < 3; i++) {
        transformed_tri.p[i] = vec3d_mult_by_mat4x4(tri.p[i], world_mat);
    }
    return transformed_tri;
}

// Function to create a translation matrix
void make_translation_matrix(float x_offset, float y_offset, float z_offset, mat4x4* translation_matrix) {
    // Initialize the translation matrix
    translation_matrix->m[0][0] = 1; translation_matrix->m[0][1] = 0; translation_matrix->m[0][2] = 0; translation_matrix->m[0][3] = 0;
    translation_matrix->m[1][0] = 0; translation_matrix->m[1][1] = 1; translation_matrix->m[1][2] = 0; translation_matrix->m[1][3] = 0;
    translation_matrix->m[2][0] = 0; translation_matrix->m[2][1] = 0; translation_matrix->m[2][2] = 1; translation_matrix->m[2][3] = 0;
    translation_matrix->m[3][0] = x_offset; translation_matrix->m[3][1] = y_offset; translation_matrix->m[3][2] = z_offset; translation_matrix->m[3][3] = 1;
}

// Function to create a rotation matrix
void make_rotation_matrix(float theta_x, float theta_y, float theta_z, mat4x4* rotation_matrix) {
    mat4x4 rot_x = {
        {
            {1, 0, 0, 0},
            {0, cosf(theta_x), -sinf(theta_x), 0},
            {0, sinf(theta_x), cosf(theta_x), 0},
            {0, 0, 0, 1}
        }
    };
    mat4x4 rot_y = {
        {
            {cosf(theta_y), 0, sinf(theta_y), 0},
            {0, 1, 0, 0},
            {-sinf(theta_y), 0, cosf(theta_y), 0},
            {0, 0, 0, 1}
        }
    };
    mat4x4 rot_z = {
        {
            {cosf(theta_z), -sinf(theta_z), 0, 0},
            {sinf(theta_z), cosf(theta_z), 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}
        }
    };

    mat4x4 temp_matrix = mat4x4_mult_by_mat4x4(rot_x, rot_y);
    *rotation_matrix = mat4x4_mult_by_mat4x4(temp_matrix, rot_z);
}

// Function to adjust coordinates for rendering
float adjust(float p, float scale, float offset) {
    // Scale and offset the coordinate
    return p * scale + offset;
}

// Function to calculate the average depth of a triangle
float tri_depth(triangle tri) {
    // Calculate the average z-coordinate of the triangle's vertices
    return (tri.p[0].z + tri.p[1].z + tri.p[2].z) / 3.0f;
}

// Function to get the surface normal of a triangle
vec3d get_surface_normal(triangle tri) {
    // Calculate two vectors along the edges of the triangle
    vec3d line_1 = {tri.p[1].x - tri.p[0].x, tri.p[1].y - tri.p[0].y, tri.p[1].z - tri.p[0].z, 1.0f};
    vec3d line_2 = {tri.p[2].x - tri.p[0].x, tri.p[2].y - tri.p[0].y, tri.p[2].z - tri.p[0].z, 1.0f};
    // Calculate the cross product to get the surface normal
    vec3d surface_normal = {
        line_1.y * line_2.z - line_1.z * line_2.y,
        line_1.z * line_2.x - line_1.x * line_2.z,
        line_1.x * line_2.y - line_1.y * line_2.x,
        1.0f
    };
    // Normalize the surface normal
    float l = sqrtf(surface_normal.x * surface_normal.x + surface_normal.y * surface_normal.y + surface_normal.z * surface_normal.z);
    surface_normal.x /= l;
    surface_normal.y /= l;
    surface_normal.z /= l;
    return surface_normal;
}

// Function to check if a triangle is visible
bool check_visibility(triangle tri, vec3d cam, vec3d* tri_normal) {
    // Get the surface normal of the triangle
    *tri_normal = get_surface_normal(tri);
    // Calculate the vector from the camera to a vertex of the triangle
    vec3d cam_to_vertex = {tri.p[0].x - cam.x, tri.p[0].y - cam.y, tri.p[0].z - cam.z, 1.0f};
    // Calculate the dot product of the surface normal and the camera-to-vertex vector
    float dot_prod = tri_normal->x * cam_to_vertex.x + tri_normal->y * cam_to_vertex.y + tri_normal->z * cam_to_vertex.z;
    // If the dot product is negative, the triangle is visible
    return dot_prod < 0;
}

// Function to illuminate a triangle
float illuminate_triangle(vec3d normalized_light_vec, vec3d tri_normal) {
    // Calculate the dot product of the light direction and the triangle normal
    return normalized_light_vec.x * tri_normal.x + normalized_light_vec.y * tri_normal.y + normalized_light_vec.z * tri_normal.z;
}

// Function to load an OBJ file
mesh load_obj(const char* obj_file) {
    // Open the OBJ file
    FILE* file = fopen(obj_file, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening OBJ file: %s\n", obj_file);
        exit(1);
    }

    // Dynamic arrays for vertices and triangles
    vec3d* vectors = NULL;
    int num_vectors = 0;
    int vectors_capacity = 0;

    int* triangles = NULL;
    int num_triangles = 0;
    int triangles_capacity = 0;

    // Read the file line by line
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        // Parse vertex data
        if (strncmp(line, "v ", 2) == 0) {
            // Resize the vectors array if necessary
            if (num_vectors >= vectors_capacity) {
                vectors_capacity = (vectors_capacity == 0) ? 16 : vectors_capacity * 2;
                vectors = (vec3d*)realloc(vectors, vectors_capacity * sizeof(vec3d));
            }
            // Read vertex coordinates
            sscanf(line, "v %f %f %f", &vectors[num_vectors].x, &vectors[num_vectors].y, &vectors[num_vectors].z);
            vectors[num_vectors].w = 1.0f;
            num_vectors++;
        // Parse face data
        } else if (strncmp(line, "f ", 2) == 0) {
            // Resize the triangles array if necessary
            if (num_triangles >= triangles_capacity) {
                triangles_capacity = (triangles_capacity == 0) ? 16 : triangles_capacity * 2;
                triangles = (int*)realloc(triangles, triangles_capacity * 3 * sizeof(int));
            }
            // Read triangle vertex indices
            sscanf(line, "f %d %d %d", &triangles[num_triangles * 3], &triangles[num_triangles * 3 + 1], &triangles[num_triangles * 3 + 2]);
            // Adjust indices to be 0-based
            triangles[num_triangles * 3]--;
            triangles[num_triangles * 3 + 1]--;
            triangles[num_triangles * 3 + 2]--;
            num_triangles++;
        }
    }
    fclose(file);

    // Create the mesh
    mesh out_mesh;
    out_mesh.num_tris = num_triangles;
    out_mesh.tris = (triangle*)malloc(num_triangles * sizeof(triangle));
    // Populate the mesh with triangles
    for (int t = 0; t < num_triangles; t++) {
        out_mesh.tris[t].p[0] = vectors[triangles[t * 3]];
        out_mesh.tris[t].p[1] = vectors[triangles[t * 3 + 1]];
        out_mesh.tris[t].p[2] = vectors[triangles[t * 3 + 2]];
        out_mesh.tris[t].c = (vec3d){0,0,0,1}; // Default color
    }

    // Free the temporary arrays
    free(vectors);
    free(triangles);
    return out_mesh;
}

// Function to draw a polygon
void draw_polygon(SDL_Renderer* renderer, int x1, int y1, int x2, int y2, int x3, int y3, int r, int g, int b) {
    // Set the drawing color
    SDL_SetRenderDrawColor(renderer, r, g, b, 255);
    // Define the vertices of the triangle
    SDL_Vertex vertices[3] = {
        {{x1, y1}, {r, g, b, 255}, {0, 0}},
        {{x2, y2}, {r, g, b, 255}, {0, 0}},
        {{x3, y3}, {r, g, b, 255}, {0, 0}}
    };
    // Draw the triangle
    SDL_RenderGeometry(renderer, NULL, vertices, 3, NULL, 0);
}

void draw_triangle_edges(SDL_Renderer* renderer, float x1, float y1, float x2, float y2, float x3, float y3, int r, int g, int b) {
    // Adjust the coordinates to the center of the screen
    x1 = adjust(x1, scale, WINDOW_SIZE / 2);
    y1 = adjust(y1, scale, WINDOW_SIZE / 2);
    x2 = adjust(x2, scale, WINDOW_SIZE / 2);
    y2 = adjust(y2, scale, WINDOW_SIZE / 2);
    x3 = adjust(x3, scale, WINDOW_SIZE / 2);
    y3 = adjust(y3, scale, WINDOW_SIZE / 2);

    // Set the draw color
    SDL_SetRenderDrawColor(renderer, r, g, b, 255);

    // Draw the lines
    SDL_RenderDrawLine(renderer, x1, y1, x2, y2);
    SDL_RenderDrawLine(renderer, x2, y2, x3, y3);
    SDL_RenderDrawLine(renderer, x3, y3, x1, y1);
}

// Comparison function for qsort
int compare_triangles(const void* a, const void* b) {
    // Cast the void pointers to triangle pointers
    const triangle* triA = (const triangle*)a;
    const triangle* triB = (const triangle*)b;
    // Calculate the depth of each triangle
    float depthA = tri_depth(*triA);
    float depthB = tri_depth(*triB);

    // Sort in descending order (farther first)
    if (depthA < depthB) return 1;
    if (depthA > depthB) return -1;
    return 0;
}

vec3d vector_intersect_plane(vec3d plane_p, vec3d plane_n, vec3d line_start, vec3d line_end) {
    plane_n = vec3d_normalize(plane_n);
    float plane_d = -vec3d_dot(plane_n, plane_p);
    float ad = vec3d_dot(line_start, plane_n);
    float bd = vec3d_dot(line_end, plane_n);
    float t = (-plane_d - ad) / (bd - ad);
    vec3d line_start_to_end = vec3d_sub(line_end, line_start);
    vec3d line_to_intersect = vec3d_mult(line_start_to_end, t);
    vec3d intersection = vec3d_add(line_start, line_to_intersect);
    return intersection;
}

// Function to calculate the signed distance from a point to a plane
float dist(vec3d p, vec3d plane_p, vec3d plane_n) {
    // vec3d n = vec3d_normalize(p); // This line is incorrect
    // return (plane_n.x * n.x + plane_n.y * n.y + plane_n.z * n.z) - vec3d_dot(plane_n, plane_p);
    return vec3d_dot(plane_n, vec3d_sub(p, plane_p)); // Corrected distance calculation
}

int triangle_clip_against_plane(vec3d plane_p, vec3d plane_n, triangle *in_tri, triangle *out_tri1, triangle *out_tri2) {
    // Make sure plane normal is normalized
    plane_n = vec3d_normalize(plane_n);

    // Create two temporary storage arrays to classify points on either side of the plane
    // If the distance sign is positive, the point is on the "inside" of the plane
    vec3d* inside_points[3]; int n_inside_points_count = 0;
    vec3d* outside_points[3]; int n_outside_points_count = 0;

    // Get signed distance of each point in triangle to plane
    float d0 = dist(in_tri->p[0], plane_p, plane_n);
    float d1 = dist(in_tri->p[1], plane_p, plane_n);
    float d2 = dist(in_tri->p[2], plane_p, plane_n);

    if (d0 >= 0) { inside_points[n_inside_points_count++] = &in_tri->p[0]; }
    else { outside_points[n_outside_points_count++] = &in_tri->p[0]; }
    if (d1 >= 0) { inside_points[n_inside_points_count++] = &in_tri->p[1]; }
    else { outside_points[n_outside_points_count++] = &in_tri->p[1]; }
    if (d2 >= 0) { inside_points[n_inside_points_count++] = &in_tri->p[2]; }
    else { outside_points[n_outside_points_count++] = &in_tri->p[2]; }

    // Classify the triangle points and split the input triangle into smaller output triangles
    // if needed. There are four possible outcomes:
    //   1. All points are inside the plane (no clipping)
    //   2. All points are outside the plane, and the whole triangle is clipped (deleted)
    //   3. Two points are inside and one point is outside (triangle becomes a quad, split into two new triangles)
    //   4. One point is inside the plane and two are outside (triangle becomes a smaller triangle)
    if (n_inside_points_count == 0) {
        // All points are outside the plane, and we delete the whole triangle
        return 0;
    }
    if (n_inside_points_count == 3) {
        // All points are inside the plane, we keep the whole triangle
        *out_tri1 = *in_tri;
        return 1;
    }
    if (n_inside_points_count == 1 && n_outside_points_count == 2) {
        // One point is inside and two points are outside
    
        // copy color info
        out_tri1->c = in_tri->c;
        // Keep inside point
        out_tri1->p[0] = *inside_points[0];
        // Find two new points where the triangle intersects the plane
        out_tri1->p[1] = vector_intersect_plane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
        out_tri1->p[2] = vector_intersect_plane(plane_p, plane_n, *inside_points[0], *outside_points[1]);
    
        return 1;
    
    }
    if (n_inside_points_count == 2 && n_outside_points_count == 1) {
        // Two points are inside and one point is outside, therefore we have a quad
    
        // copy color info for two new triangles
        out_tri1->c = in_tri->c;
        out_tri2->c = in_tri->c;
    
        // Create the first triangle with the two inside points, and the intercept
        // between the outside point and the plane
        out_tri1->p[0] = *inside_points[0];
        out_tri1->p[1] = *inside_points[1];
        out_tri1->p[2] = vector_intersect_plane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
        // Create the second triangle using the second inside point, the newly generated point on triangle 1,
        // and the intercept between the other outside point and the plane
        out_tri2->p[0] = *inside_points[1];
        out_tri2->p[1] = out_tri1->p[2];
        out_tri2->p[2] = vector_intersect_plane(plane_p, plane_n, *inside_points[1], *outside_points[0]);
        return 2;
    }
}

// Function to initialize the dynamic array
void dynamicArrayInit(DynamicArray* arr, size_t element_size, int initialCapacity) {
    arr->data = malloc(initialCapacity * element_size);
    if (arr->data == NULL) {
        fprintf(stderr, "Memory allocation failed in dynamicArrayInit\n");
        exit(1);
    }
    arr->element_size = element_size;
    arr->size = 0;
    arr->capacity = initialCapacity;
}

// Function to add an element to the end of the array (like push_back)
void dynamicArrayPushBack(DynamicArray* arr, const void* value) {
    if (arr->size == arr->capacity) {
        arr->capacity *= 2;
        void* newData = realloc(arr->data, arr->capacity * arr->element_size);
        if (newData == NULL) {
            fprintf(stderr, "Memory reallocation failed in dynamicArrayPushBack\n");
            exit(1);
        }
        arr->data = newData;
    }

    // Calculate the memory location of the new element
    char* target = (char*)arr->data + arr->size * arr->element_size;

    // Copy the new element into the array
    memcpy(target, value, arr->element_size);

    arr->size++;
}

// Function to remove and return the last element from the array
void* dynamicArrayPopFront(DynamicArray* arr) {
    // Check if the array is empty
    if (arr == NULL || arr->size == 0) {
        fprintf(stderr, "Error: Cannot pop from an empty or NULL array\n");
        return NULL;
    }

    // 1. Allocate memory for the element to be returned
    void* return_val = malloc(arr->element_size);
    if (return_val == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for popped element\n");
        return NULL; // Memory allocation failed
    }

    // 2. Copy the first element's data into the return buffer
    memcpy(return_val, arr->data, arr->element_size);

    // 3. Shift all subsequent elements one position to the left
    // Calculate the starting address of the second element
    char* source = (char*)arr->data + arr->element_size;
    // Calculate the number of bytes to move (all elements except the first)
    size_t bytes_to_move = (arr->size - 1) * arr->element_size;

    // Use memmove for safe shifting (handles overlapping memory)
    // The destination is the start of the array data
    if (arr->size > 1) { // Only shift if there's more than one element
         memmove(arr->data, source, bytes_to_move);
    }

    // 4. Decrease the array size
    arr->size--;

    // Optional: Consider shrinking the array if size becomes much smaller than capacity
    // (This adds complexity and might not always be desired due to realloc cost)
    // if (arr->size > 0 && arr->size <= arr->capacity / 4) { ... shrink logic ... }

    // 5. Return the pointer to the copied first element
    return return_val;
}

// Function to free the allocated memory
void dynamicArrayFree(DynamicArray* arr) {
    free(arr->data);
    arr->data = NULL;
    arr->size = 0;
    arr->capacity = 0;
    arr->element_size = 0;
}

int main(int argc, char* argv[]) {
    // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "Could not initialize SDL: %s\n", SDL_GetError());
        return 1;
    }

    // Create a window
    SDL_Window* window = SDL_CreateWindow("3D Renderer", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WINDOW_SIZE, WINDOW_SIZE, SDL_WINDOW_SHOWN);
    if (window == NULL) {
        fprintf(stderr, "Could not create window: %s\n", SDL_GetError());
        return 1;
    }

    // Create a renderer
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
    if (renderer == NULL) {
        fprintf(stderr, "Could not create renderer: %s\n", SDL_GetError());
        return 1;
    }

    // Load the OBJ file
    mesh cube = load_obj("mountains.obj");

    // Normalize the light direction
    vec3d normalized_light_direction = vec3d_normalize(light_direction);

    // Calculate f_fov_rad
    f_fov_rad = 1.0f / tanf(f_fov * 0.5f / 180.0f * M_PI);

    bool quit = false;
    SDL_Event event;

    // Allocate memory for triangles to raster
    triangle* triangles_to_raster = (triangle*)malloc(cube.num_tris * sizeof(triangle));
    if (triangles_to_raster == NULL) {
        fprintf(stderr, "Failed to allocate memory for triangles_to_raster\n");
        return 1;
    }

    // Main loop
    while (!quit) {
        int num_triangles_to_raster = 0;
        // Handle events
        while (SDL_PollEvent(&event) != 0) {
            if (event.type == SDL_QUIT) {
                quit = true;
            } else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym < NUM_KEYS) {
                    key_states[event.key.keysym.sym] = true;
                }
            } else if (event.type == SDL_KEYUP) {
                if (event.key.keysym.sym < NUM_KEYS) {
                    key_states[event.key.keysym.sym] = false;
                }
            }
        }
        // Create rotation and translation matrices
        mat4x4 rotation_matrix;
        mat4x4 translation_matrix;
        mat4x4 world_matrix;
        mat4x4 camera_rotation_matrix;

        make_rotation_matrix(angle_x, angle_y, angle_z, &rotation_matrix);
        make_translation_matrix(0, 0, obj_distance, &translation_matrix);

        world_matrix = mat4x4_mult_by_mat4x4(rotation_matrix, translation_matrix);

        // Create Point At matrix for the camera
        // vec3d v_look_dir = {0, 0, 1};
        vec3d v_up = {0, 1, 0, 1};
        vec3d v_target = {0, 0, 1, 1};
        make_rotation_matrix(angle_x, angle_y, angle_z, &camera_rotation_matrix);
        vec3d v_look_dir = vec3d_mult_by_mat4x4(v_target, camera_rotation_matrix);
        // vec3d v_target = vec3d_add(viewer_Camera, v_look_dir);
        mat4x4 camera_matrix = matrix_point_at(viewer_Camera, v_target, v_up);
        // Invert the camera matrix
        mat4x4 view_matrix = Matrix_QuickInverse(camera_matrix);  

        // Camera Movement
        float camera_speed = 8.0f; // Adjust this for faster/slower movement
        float fElapsedTime = 1.0f / TICK; // Time since last frame

        vec3d vector_forward = vec3d_mult(v_look_dir, 8.0f * fElapsedTime);

		if (GetKey(SDLK_w)) {
			viewer_Camera = vec3d_add(viewer_Camera, vector_forward);
		}

		if (GetKey(SDLK_s)) {
			viewer_Camera = vec3d_sub(viewer_Camera, vector_forward);
		}

		if (GetKey(SDLK_a)) {
            angle_y -= 2.0f * fElapsedTime;
        }

		if (GetKey(SDLK_d)) {
			angle_y += 2.0f * fElapsedTime;
        }

        if (GetKey(SDLK_r)) {
            viewer_Camera.y += camera_speed * fElapsedTime; // Move up
        }

        if (GetKey(SDLK_f)) {
            viewer_Camera.y -= camera_speed * fElapsedTime; // Move down
        }
        // Process each triangle in the mesh
        for (int i = 0; i < cube.num_tris; i++) {
            // Transform the triangle by the world matrix
            triangle transformed_triangle = transform_triangle(cube.tris[i], world_matrix);
            vec3d tri_normal;
            // Check if the triangle is visible
            bool visible = check_visibility(transformed_triangle, viewer_Camera, &tri_normal);
            // Convert World space to View space
            triangle viewed_triangle = transform_triangle(transformed_triangle, view_matrix);
            // Clip the triangle against the near plane
            int num_clipped_triangles = 0;
            triangle clipped[2];

            num_clipped_triangles = triangle_clip_against_plane((vec3d){0.0f, 0.0f, f_near}, // Use f_near
                                                                (vec3d){0.0f, 0.0f, 1.0f}, // Normal points +Z (into scene)
                                                                &viewed_triangle,
                                                                &clipped[0],
                                                                &clipped[1]);
              
            for (int n = 0; n < num_clipped_triangles; n++) {

                // Project the triangle to 2D
                triangle projected_triangle = project_triangle(clipped[n]);
                // Divide the triangle by w
                projected_triangle = tri_div_by_w(projected_triangle);

                // If the triangle is visible
                if (visible) {
                    // Calculate the illumination of the triangle
                    float illumination = illuminate_triangle(normalized_light_direction, tri_normal);
                    // Set the color of the triangle based on illumination
                    projected_triangle.c = (vec3d){light_color.x * illumination, light_color.y * illumination, light_color.z * illumination, 1.0f};
                    // Add the triangle to the raster list
                    triangles_to_raster[num_triangles_to_raster++] = projected_triangle;
                }
            }
        }
        // printf("Sorting Triangles\n");
        // printf("Drawing %d triangles\n", num_triangles_to_raster);
        // Sort triangles by depth
        qsort(triangles_to_raster, num_triangles_to_raster, sizeof(triangle), compare_triangles);

        // Clear the screen
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // Draw each triangle
        // printf("Drawing %d triangles\n", num_triangles_to_raster);


        for (int i = 0; i < num_triangles_to_raster; i++) {
            triangle clipped[2]; // Temporary storage for clipped triangles
            DynamicArray list_triangles;

            // Initialize the list for the current triangle from the raster list
            dynamicArrayInit(&list_triangles, sizeof(triangle), 2); // Start small, it will grow
            dynamicArrayPushBack(&list_triangles, &triangles_to_raster[i]);

            // Clip against screen edges using the queue-like approach
            for (int p = 0; p < 4; p++) {
                int nTrisToAdd = 0;
                // Get the number of triangles to process *for this specific plane*
                int triangles_to_process_this_plane = list_triangles.size;
                // printf("Processing %d triangles for plane %d\n", triangles_to_process_this_plane, p);
                // Process only the triangles currently in the list for this plane
                while (triangles_to_process_this_plane > 0) {
                    // printf("Processing triangle %d\n", triangles_to_process_this_plane);
                    // Take triangle from front of queue (popFront returns a malloc'd copy)
                    triangle* test_tri = (triangle*)dynamicArrayPopFront(&list_triangles);
                    if (test_tri == NULL) {
                        // This should ideally not happen if triangles_to_process_this_plane > 0
                        fprintf(stderr, "Error: dynamicArrayPopFront returned NULL unexpectedly during clipping.\n");
                        triangles_to_process_this_plane--; // Decrement anyway to avoid infinite loop
                        continue;
                    }
                    triangles_to_process_this_plane--; // Decrement the count for this plane pass

                    // Define the clipping plane based on the loop variable 'p'
                    vec3d plane_p, plane_n;
                    switch (p) {
                        case 0: // Bottom edge (Plane Y=0, Normal points "inwards" -Y, axis is reversed)
                            plane_p = (vec3d){0.0f, 2.0f, 0.0f, 1.0f};
                            plane_n = (vec3d){0.0f, -1.0f, 0.0f, 1.0f};
                            break;
                        case 1: // Top edge (Plane Y=WINDOW_SIZE-1, Normal points "inwards" +Y, axis is reversed)
                            plane_p = (vec3d){0.0f, -2.0f, 0.0f, 1.0f};
                            plane_n = (vec3d){0.0f, +1.0f, 0.0f, 1.0f};
                            break;
                        case 2: // Left edge (Plane X=0, Normal points "inwards" +X)
                            plane_p = (vec3d){-2.0f, 0.0f, 0.0f, 1.0f};
                            plane_n = (vec3d){1.0f, 0.0f, 0.0f, 1.0f};
                            break;
                        case 3: // Right edge (Plane X=WINDOW_SIZE-1, Normal points "inwards" -X)
                            plane_p = (vec3d){2.0f, 0.0f, 0.0f, 1.0f};
                            plane_n = (vec3d){-1.0f, 0.0f, 0.0f, 1.0f};
                            break;
                        default: // Should not happen
                            free(test_tri); // Free memory before continuing
                            continue;
                    }
                    // printf("Num Tris to add: %d\n", nTrisToAdd);
                    // Clip the popped triangle against the current plane
                    nTrisToAdd = triangle_clip_against_plane(plane_p, plane_n,
                                                             test_tri, // Pass the pointer to the popped triangle data
                                                             &clipped[0],
                                                             &clipped[1]);

                    // Add the resulting clipped triangle(s) back to the *end* of the list
                    // These will be processed by the *next* plane (p+1)
                    for (int w = 0; w < nTrisToAdd; w++) {
                        dynamicArrayPushBack(&list_triangles, &clipped[w]);
                    }

                    // Free the memory allocated by dynamicArrayPopFront for the triangle we just processed
                    free(test_tri);
                }
                // After the while loop, all triangles that existed *before* this plane's
                // clipping pass have been processed. The list now contains the results,
                // ready for the next plane.
                // If the list becomes empty after clipping against a plane, we can stop early.
                 if (list_triangles.size == 0) {
                    printf("All triangles clipped out at plane %d\n", p);
                    break; // Stop clipping this original triangle if it's fully outside a plane
                 }
            } // End of plane clipping loop (p)

            // Draw the final resulting triangles from the list after all clipping planes
            // printf("Drawing %d triangles\n", list_triangles.size);
            for (int k = 0; k < list_triangles.size; k++) {
                // Access the triangle data correctly using pointer arithmetic
                triangle* tri_to_draw = (triangle*)((char*)list_triangles.data + k * list_triangles.element_size);
                // Draw the filled polygon
                draw_polygon(renderer,
                    adjust(tri_to_draw->p[0].x, scale, WINDOW_SIZE / 2),
                    adjust(tri_to_draw->p[0].y, scale, WINDOW_SIZE / 2),
                    adjust(tri_to_draw->p[1].x, scale, WINDOW_SIZE / 2),
                    adjust(tri_to_draw->p[1].y, scale, WINDOW_SIZE / 2),
                    adjust(tri_to_draw->p[2].x, scale, WINDOW_SIZE / 2),
                    adjust(tri_to_draw->p[2].y, scale, WINDOW_SIZE / 2),
                    (int)fmax(0.0f, fmin(255.0f, tri_to_draw->c.x)), // Clamp color values
                    (int)fmax(0.0f, fmin(255.0f, tri_to_draw->c.y)),
                    (int)fmax(0.0f, fmin(255.0f, tri_to_draw->c.z)));

                // Draw the edges (optional)
                // Note: Edges are drawn *before* the adjust function is applied internally
                draw_triangle_edges(renderer,
                                    tri_to_draw->p[0].x, tri_to_draw->p[0].y,
                                    tri_to_draw->p[1].x, tri_to_draw->p[1].y,
                                    tri_to_draw->p[2].x, tri_to_draw->p[2].y,
                                    255, 255, 255); // White edges
            }
            dynamicArrayFree(&list_triangles); // Free the list's internal data for this original triangle
        } // End of raster loop (i)


        // Update the screen
        SDL_RenderPresent(renderer);

        // Simulate frame update
        SDL_Delay(1000 / TICK);
    }

    // Clean up
    free(triangles_to_raster);
    free(cube.tris);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
