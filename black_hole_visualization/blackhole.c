#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Window Settings
#define WIDTH 100
#define HEIGHT 75

// Physics Constants
#define RS 2.0              // Schwarzschild Radius
#define MAX_STEPS 2000      // Max integration steps
#define STEP_SIZE 0.05      // RK4 Delta
#define ESCAPE_RADIUS 25.0  // Distance considered "safe"

// Accretion Disk Constants
#define DISK_INNER 3.0      // Innermost stable orbit approx
#define DISK_OUTER 8.0
#define DISK_Color_R 1.0
#define DISK_Color_G 0.6
#define DISK_Color_B 0.1

// Camera State (Spherical Coordinates)
double cam_dist = 14.0;
double cam_theta = 1.57; // Start at equator (~PI/2)
double cam_phi   = 0.0;

// Pixel Buffer
float *pixel_buffer;

// Structs
typedef struct { float r, g, b; } Color;
typedef struct { double x, y, z; } Vec3;

// State: Spherical Position and Velocity
typedef struct {
    double r, theta, phi;       // Pos
    double pr, ptheta, pphi;    // Velocity (derivatives w.r.t affine param)
} State;

// ---------------------------------------------------------------------------
// Math Helpers
// ---------------------------------------------------------------------------
Vec3 normalize(Vec3 v) {
    double m = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    Vec3 r = {v.x/m, v.y/m, v.z/m};
    return r;
}

Vec3 cross(Vec3 a, Vec3 b) {
    Vec3 r;
    r.x = a.y*b.z - a.z*b.y;
    r.y = a.z*b.x - a.x*b.z;
    r.z = a.x*b.y - a.y*b.x;
    return r;
}

// ---------------------------------------------------------------------------
// Physics Engine (Schwarzschild Metric)
// ---------------------------------------------------------------------------

State get_derivatives(State s) {
    State d;
    double r = s.r;
    double r_minus_rs = r - RS;
    if (r_minus_rs < 1e-5) r_minus_rs = 1e-5;

    double sin_t = sin(s.theta);
    double cos_t = cos(s.theta);
    
    // Velocities
    d.r = s.pr;
    d.theta = s.ptheta;
    d.phi = s.pphi;

    // Time component constraint (null geodesic ds^2=0)
    double time_metric = 1.0 - (RS / r);
    double term_r = (s.pr * s.pr) / time_metric;
    double term_ang = (r*r) * (s.ptheta*s.ptheta + sin_t*sin_t * s.pphi*s.pphi);
    double dt_sq = (term_r + term_ang) / time_metric;

    // Accelerations (Christoffel Symbols)
    // r-acceleration
    d.pr = (RS / (2.0 * r * r_minus_rs)) * (s.pr * s.pr)
         + r_minus_rs * (s.ptheta*s.ptheta + sin_t*sin_t * s.pphi*s.pphi)
         - (RS * r_minus_rs / (2.0 * r*r*r)) * dt_sq;

    // theta-acceleration
    d.ptheta = -(2.0/r)*s.pr*s.ptheta + sin_t*cos_t*(s.pphi*s.pphi);

    // phi-acceleration
    d.pphi = -(2.0/r)*s.pr*s.pphi - 2.0*(cos_t/sin_t)*s.ptheta*s.pphi;

    return d;
}

State rk4_step(State s, double h) {
    State k1, k2, k3, k4, tmp;
    
    k1 = get_derivatives(s);

    tmp = s; 
    tmp.r += 0.5*h*k1.r; tmp.theta += 0.5*h*k1.theta; tmp.phi += 0.5*h*k1.phi;
    tmp.pr += 0.5*h*k1.pr; tmp.ptheta += 0.5*h*k1.ptheta; tmp.pphi += 0.5*h*k1.pphi;
    k2 = get_derivatives(tmp);

    tmp = s;
    tmp.r += 0.5*h*k2.r; tmp.theta += 0.5*h*k2.theta; tmp.phi += 0.5*h*k2.phi;
    tmp.pr += 0.5*h*k2.pr; tmp.ptheta += 0.5*h*k2.ptheta; tmp.pphi += 0.5*h*k2.pphi;
    k3 = get_derivatives(tmp);

    tmp = s;
    tmp.r += h*k3.r; tmp.theta += h*k3.theta; tmp.phi += h*k3.phi;
    tmp.pr += h*k3.pr; tmp.ptheta += h*k3.ptheta; tmp.pphi += h*k3.pphi;
    k4 = get_derivatives(tmp);

    State res = s;
    res.r += (h/6.0)*(k1.r + 2*k2.r + 2*k3.r + k4.r);
    res.theta += (h/6.0)*(k1.theta + 2*k2.theta + 2*k3.theta + k4.theta);
    res.phi += (h/6.0)*(k1.phi + 2*k2.phi + 2*k3.phi + k4.phi);
    
    res.pr += (h/6.0)*(k1.pr + 2*k2.pr + 2*k3.pr + k4.pr);
    res.ptheta += (h/6.0)*(k1.ptheta + 2*k2.ptheta + 2*k3.ptheta + k4.ptheta);
    res.pphi += (h/6.0)*(k1.pphi + 2*k2.pphi + 2*k3.pphi + k4.pphi);
    return res;
}

// ---------------------------------------------------------------------------
// Ray Tracing
// ---------------------------------------------------------------------------

void trace_rays() {
    printf("Rendering... (Theta: %.2f, Phi: %.2f)\n", cam_theta, cam_phi);

    // 1. Calculate Camera Basis Vectors in Cartesian
    double cx = cam_dist * sin(cam_theta) * cos(cam_phi);
    double cy = cam_dist * sin(cam_theta) * sin(cam_phi);
    double cz = cam_dist * cos(cam_theta);
    
    // Forward Vector (Looking at origin)
    Vec3 F = normalize((Vec3){-cx, -cy, -cz});
    
    // Up Vector
    Vec3 globalUp = {0,0,1};
    if (fabs(fabs(cam_theta) - 0.0) < 0.1 || fabs(fabs(cam_theta) - M_PI) < 0.1) {
        globalUp.x = 1; globalUp.z = 0; // Fix gimbal lock at poles
    }

    Vec3 R = normalize(cross(F, globalUp)); // Right
    Vec3 U = normalize(cross(R, F));        // Up

    for(int y = 0; y < HEIGHT; y++) {
        for(int x = 0; x < WIDTH; x++) {
            
            // Normalized Device Coords
            double u = (x - WIDTH/2.0) / (WIDTH/2.0);
            double v = (y - HEIGHT/2.0) / (HEIGHT/2.0);

            // Ray Direction in Cartesian
            double fov = 1.0;
            double ratio = (double)WIDTH/HEIGHT;
            
            Vec3 dir;
            dir.x = F.x + R.x * u * fov * ratio + U.x * v * fov;
            dir.y = F.y + R.y * u * fov * ratio + U.y * v * fov;
            dir.z = F.z + R.z * u * fov * ratio + U.z * v * fov;
            dir = normalize(dir);

            // Initial State (Spherical)
            State s;
            s.r = cam_dist;
            s.theta = cam_theta;
            s.phi = cam_phi;

            // Convert Cartesian Velocity to Spherical Velocity
            double rho = cx*cx + cy*cy; 
            double r_sq = rho + cz*cz;
            double r_val = sqrt(r_sq);
            double sqrt_rho = sqrt(rho);

            if (sqrt_rho < 1e-5) sqrt_rho = 1e-5;

            // v_r
            s.pr = (cx*dir.x + cy*dir.y + cz*dir.z) / r_val;
            
            // v_theta
            double num_theta = cz*(cx*dir.x + cy*dir.y) - rho*dir.z;
            s.ptheta = num_theta / (r_sq * sqrt_rho);

            // v_phi
            s.pphi = (cx*dir.y - cy*dir.x) / rho;

            // Trace Loop
            int hitType = 0; // 0=Space, 1=BlackHole, 2=Disk

            for(int k=0; k<MAX_STEPS; k++) {
                State prev = s;
                s = rk4_step(s, STEP_SIZE);

                // Check Accretion Disk (Crossing the Equator plane)
                if ((prev.theta - M_PI_2) * (s.theta - M_PI_2) < 0) {
                    double f = fabs(prev.theta - M_PI_2) / fabs(prev.theta - s.theta);
                    double r_cross = prev.r + f * (s.r - prev.r);

                    if (r_cross > DISK_INNER && r_cross < DISK_OUTER) {
                        hitType = 2;
                        break;
                    }
                }

                // Check Event Horizon
                if (s.r < RS * 1.01) {
                    hitType = 1;
                    break;
                }

                // Check Escape
                if (s.r > ESCAPE_RADIUS) {
                    hitType = 0;
                    break;
                }
            }

            // Coloring
            int p = (y * WIDTH + x) * 3;
            if (hitType == 1) {
                // Black Hole (Black)
                pixel_buffer[p] = 0.0f; pixel_buffer[p+1] = 0.0f; pixel_buffer[p+2] = 0.0f;
            } 
            else if (hitType == 2) {
                // Disk (Orange) with gradient
                float hotness = (DISK_OUTER - s.r) / (DISK_OUTER - DISK_INNER);
                pixel_buffer[p]   = DISK_Color_R * hotness + 0.2; 
                pixel_buffer[p+1] = DISK_Color_G * hotness; 
                pixel_buffer[p+2] = DISK_Color_B * 0.1;
            } 
            else {
                // Background (Dark Grey, No Stars)
                pixel_buffer[p]   = 0.1f;
                pixel_buffer[p+1] = 0.1f;
                pixel_buffer[p+2] = 0.15f;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// OpenGL / Input
// ---------------------------------------------------------------------------

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_FLOAT, pixel_buffer);
    glFlush();
}

void processKeys(int key, int x, int y) {
    const double angle_step = 0.15;
    switch(key) {
        case GLUT_KEY_LEFT:  cam_phi -= angle_step; break;
        case GLUT_KEY_RIGHT: cam_phi += angle_step; break;
        case GLUT_KEY_UP:    cam_theta -= angle_step; if(cam_theta < 0.1) cam_theta = 0.1; break;
        case GLUT_KEY_DOWN:  cam_theta += angle_step; if(cam_theta > 3.0) cam_theta = 3.0; break;
    }
    trace_rays();
    glutPostRedisplay();
}

void normalKeys(unsigned char key, int x, int y) {
    if (key == 'r') {
        cam_theta = 1.57; cam_phi = 0.0;
        trace_rays();
        glutPostRedisplay();
    }
}

void init() {
    pixel_buffer = (float*)malloc(WIDTH * HEIGHT * 3 * sizeof(float));
    trace_rays();
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Black Hole Visualization");
    
    init();
    
    glutDisplayFunc(display);
    glutSpecialFunc(processKeys);
    glutKeyboardFunc(normalKeys);
    
    printf("Controls:\n  Arrow Keys: Orbit\n  R: Reset\n");
    
    glutMainLoop();
    free(pixel_buffer);
    return 0;
}
