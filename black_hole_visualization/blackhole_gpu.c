// Define this for Linux to access Shader functions
#define GL_GLEXT_PROTOTYPES 

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glext.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Window Settings (GPU is fast, we can do high res now!)
#define WIDTH 800
#define HEIGHT 600

// Camera State (CPU Side)
double cam_dist = 14.0;
double cam_theta = 1.57; 
double cam_phi   = 0.0;

// OpenGL handles
GLuint shaderProgram;
GLint loc_res, loc_camPos, loc_camDir, loc_camRight, loc_camUp;

// ---------------------------------------------------------------------------
// THE FRAGMENT SHADER (This runs on the GPU for every pixel)
// ---------------------------------------------------------------------------
const char* fragmentShaderSource = 
"#version 120\n"
"uniform vec2 u_resolution;\n"
"uniform vec3 u_camPos;    // Spherical: r, theta, phi\n"
"uniform vec3 u_camFwd;\n"
"uniform vec3 u_camRight;\n"
"uniform vec3 u_camUp;\n"
"\n"
"const float RS = 2.0;\n"
"const int MAX_STEPS = 500;\n" // Can lower steps slightly for GPU, or keep high
"const float STEP_SIZE = 0.05;\n"
"const float ESCAPE_RADIUS = 25.0;\n"
"const float DISK_INNER = 3.0;\n"
"const float DISK_OUTER = 8.0;\n"
"const float PI = 3.14159265;\n"
"\n"
"// Struct to hold state\n"
"struct State {\n"
"    float r, theta, phi;\n"
"    float pr, ptheta, pphi;\n"
"};\n"
"\n"
"// Get derivatives (Schwarzschild Geodesic)\n"
"State get_derivatives(State s) {\n"
"    State d;\n"
"    float r = s.r;\n"
"    float r_minus_rs = r - RS;\n"
"    if (r_minus_rs < 0.001) r_minus_rs = 0.001;\n"
"\n"
"    float sin_t = sin(s.theta);\n"
"    float cos_t = cos(s.theta);\n"
"\n"
"    // Velocities\n"
"    d.r = s.pr;\n"
"    d.theta = s.ptheta;\n"
"    d.phi = s.pphi;\n"
"\n"
"    // Time component dt/dlambda\n"
"    float time_metric = 1.0 - (RS / r);\n"
"    float term_r = (s.pr * s.pr) / time_metric;\n"
"    float term_ang = (r*r) * (s.ptheta*s.ptheta + sin_t*sin_t * s.pphi*s.pphi);\n"
"    float dt_sq = (term_r + term_ang) / time_metric;\n"
"\n"
"    // Accelerations\n"
"    d.pr = (RS / (2.0 * r * r_minus_rs)) * (s.pr * s.pr)\n"
"         + r_minus_rs * (s.ptheta*s.ptheta + sin_t*sin_t * s.pphi*s.pphi)\n"
"         - (RS * r_minus_rs / (2.0 * r*r*r)) * dt_sq;\n"
"\n"
"    d.ptheta = -(2.0/r)*s.pr*s.ptheta + sin_t*cos_t*(s.pphi*s.pphi);\n"
"    d.pphi = -(2.0/r)*s.pr*s.pphi - 2.0*(cos_t/sin_t)*s.ptheta*s.pphi;\n"
"\n"
"    return d;\n"
"}\n"
"\n"
"void main() {\n"
"    // 1. Calculate Ray Direction\n"
"    vec2 uv = (gl_FragCoord.xy - 0.5 * u_resolution.xy) / (u_resolution.y * 0.5);\n"
"    vec3 dir = normalize(u_camFwd + uv.x * u_camRight + uv.y * u_camUp);\n"
"\n"
"    // 2. Initial State (Spherical)\n"
"    State s;\n"
"    s.r = u_camPos.x;     // r\n"
"    s.theta = u_camPos.y; // theta\n"
"    s.phi = u_camPos.z;   // phi\n"
"\n"
"    // 3. Convert Cartesian Velocity (dir) to Spherical Velocity\n"
"    // Reconstruct Cartesian Cam Pos for basis conversion\n"
"    float cx = s.r * sin(s.theta) * cos(s.phi);\n"
"    float cy = s.r * sin(s.theta) * sin(s.phi);\n"
"    float cz = s.r * cos(s.theta);\n"
"\n"
"    float rho = cx*cx + cy*cy;\n"
"    float r_sq = rho + cz*cz;\n"
"    float r_val = sqrt(r_sq);\n"
"    float sqrt_rho = sqrt(rho);\n"
"    if(sqrt_rho < 0.001) sqrt_rho = 0.001;\n"
"\n"
"    s.pr = (cx*dir.x + cy*dir.y + cz*dir.z) / r_val;\n"
"    float num_theta = cz*(cx*dir.x + cy*dir.y) - rho*dir.z;\n"
"    s.ptheta = num_theta / (r_sq * sqrt_rho);\n"
"    s.pphi = (cx*dir.y - cy*dir.x) / rho;\n"
"\n"
"    // 4. Trace Loop\n"
"    int hitType = 0; // 0=Space, 1=BH, 2=Disk\n"
"    float final_r = 0.0;\n"
"\n"
"    for(int i=0; i<MAX_STEPS; i++) {\n"
"        State k1 = get_derivatives(s);\n"
"        \n"
"        State tmp = s;\n"
"        tmp.r += 0.5*STEP_SIZE*k1.r; tmp.theta += 0.5*STEP_SIZE*k1.theta; tmp.phi += 0.5*STEP_SIZE*k1.phi;\n"
"        tmp.pr += 0.5*STEP_SIZE*k1.pr; tmp.ptheta += 0.5*STEP_SIZE*k1.ptheta; tmp.pphi += 0.5*STEP_SIZE*k1.pphi;\n"
"        State k2 = get_derivatives(tmp);\n"
"\n"
"        tmp = s;\n"
"        tmp.r += 0.5*STEP_SIZE*k2.r; tmp.theta += 0.5*STEP_SIZE*k2.theta; tmp.phi += 0.5*STEP_SIZE*k2.phi;\n"
"        tmp.pr += 0.5*STEP_SIZE*k2.pr; tmp.ptheta += 0.5*STEP_SIZE*k2.ptheta; tmp.pphi += 0.5*STEP_SIZE*k2.pphi;\n"
"        State k3 = get_derivatives(tmp);\n"
"\n"
"        tmp = s;\n"
"        tmp.r += STEP_SIZE*k3.r; tmp.theta += STEP_SIZE*k3.theta; tmp.phi += STEP_SIZE*k3.phi;\n"
"        tmp.pr += STEP_SIZE*k3.pr; tmp.ptheta += STEP_SIZE*k3.ptheta; tmp.pphi += STEP_SIZE*k3.pphi;\n"
"        State k4 = get_derivatives(tmp);\n"
"\n"
"        // Save prev for disk check\n"
"        float prev_theta = s.theta;\n"
"        float prev_r = s.r;\n"
"\n"
"        // Advance\n"
"        s.r += (STEP_SIZE/6.0)*(k1.r + 2.0*k2.r + 2.0*k3.r + k4.r);\n"
"        s.theta += (STEP_SIZE/6.0)*(k1.theta + 2.0*k2.theta + 2.0*k3.theta + k4.theta);\n"
"        s.phi += (STEP_SIZE/6.0)*(k1.phi + 2.0*k2.phi + 2.0*k3.phi + k4.phi);\n"
"        s.pr += (STEP_SIZE/6.0)*(k1.pr + 2.0*k2.pr + 2.0*k3.pr + k4.pr);\n"
"        s.ptheta += (STEP_SIZE/6.0)*(k1.ptheta + 2.0*k2.ptheta + 2.0*k3.ptheta + k4.ptheta);\n"
"        s.pphi += (STEP_SIZE/6.0)*(k1.pphi + 2.0*k2.pphi + 2.0*k3.pphi + k4.pphi);\n"
"\n"
"        // Check Disk (Crossing PI/2)\n"
"        if ((prev_theta - 1.5707) * (s.theta - 1.5707) < 0.0) {\n"
"            float f = abs(prev_theta - 1.5707) / abs(prev_theta - s.theta);\n"
"            float r_cross = prev_r + f * (s.r - prev_r);\n"
"            if (r_cross > DISK_INNER && r_cross < DISK_OUTER) {\n"
"                hitType = 2;\n"
"                final_r = r_cross;\n"
"                break;\n"
"            }\n"
"        }\n"
"\n"
"        // Check Horizon\n"
"        if (s.r < RS * 1.01) {\n"
"            hitType = 1;\n"
"            break;\n"
"        }\n"
"        // Check Escape\n"
"        if (s.r > ESCAPE_RADIUS) {\n"
"            hitType = 0;\n"
"            break;\n"
"        }\n"
"    }\n"
"\n"
"    // 5. Coloring\n"
"    if (hitType == 1) {\n"
"        gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);\n"
"    } else if (hitType == 2) {\n"
"        float hotness = (DISK_OUTER - final_r) / (DISK_OUTER - DISK_INNER);\n"
"        vec3 col = vec3(1.0 * hotness + 0.2, 0.6 * hotness, 0.1 * hotness);\n"
"        gl_FragColor = vec4(col, 1.0);\n"
"    } else {\n"
"        gl_FragColor = vec4(0.05, 0.05, 0.1, 1.0);\n"
"    }\n"
"}\n";

// ---------------------------------------------------------------------------
// Host Side Helpers
// ---------------------------------------------------------------------------

typedef struct { double x, y, z; } Vec3;

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

void initShader() {
    // 1. Create Shader Objects
    GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    
    // 2. Load Source
    glShaderSource(fragShader, 1, &fragmentShaderSource, NULL);
    
    // 3. Compile
    glCompileShader(fragShader);
    
    // Debug Compilation
    GLint success;
    glGetShaderiv(fragShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(fragShader, 512, NULL, infoLog);
        printf("Shader Compile Error:\n%s\n", infoLog);
    }
    
    // 4. Create Program
    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, fragShader);
    glLinkProgram(shaderProgram);
    
    // 5. Get Uniform Locations
    loc_res = glGetUniformLocation(shaderProgram, "u_resolution");
    loc_camPos = glGetUniformLocation(shaderProgram, "u_camPos");
    loc_camDir = glGetUniformLocation(shaderProgram, "u_camFwd");
    loc_camRight = glGetUniformLocation(shaderProgram, "u_camRight");
    loc_camUp = glGetUniformLocation(shaderProgram, "u_camUp");
    
    printf("Shader initialized.\n");
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    
    glUseProgram(shaderProgram);
    
    // Update Resolution
    glUniform2f(loc_res, (float)WIDTH, (float)HEIGHT);
    
    // Calculate Camera Vectors on CPU
    double cx = cam_dist * sin(cam_theta) * cos(cam_phi);
    double cy = cam_dist * sin(cam_theta) * sin(cam_phi);
    double cz = cam_dist * cos(cam_theta);
    
    Vec3 camPosCart = {cx, cy, cz};
    Vec3 F = normalize((Vec3){-cx, -cy, -cz});
    
    Vec3 globalUp = {0,0,1};
    if (fabs(fabs(cam_theta) - 0.0) < 0.1 || fabs(fabs(cam_theta) - M_PI) < 0.1) {
        globalUp.x = 1; globalUp.z = 0;
    }
    Vec3 R = normalize(cross(F, globalUp));
    Vec3 U = normalize(cross(R, F));
    
    // Send Uniforms
    glUniform3f(loc_camPos, (float)cam_dist, (float)cam_theta, (float)cam_phi);
    glUniform3f(loc_camDir, (float)F.x, (float)F.y, (float)F.z);
    glUniform3f(loc_camRight, (float)R.x, (float)R.y, (float)R.z);
    glUniform3f(loc_camUp, (float)U.x, (float)U.y, (float)U.z);

    // Draw Full Screen Quad
    glRectf(-1.0f, -1.0f, 1.0f, 1.0f);
    
    glUseProgram(0);
    glFlush();
}

void processKeys(int key, int x, int y) {
    const double angle_step = 0.05; // Smoother controls
    switch(key) {
        case GLUT_KEY_LEFT:  cam_phi -= angle_step; break;
        case GLUT_KEY_RIGHT: cam_phi += angle_step; break;
        case GLUT_KEY_UP:    cam_theta -= angle_step; if(cam_theta < 0.1) cam_theta = 0.1; break;
        case GLUT_KEY_DOWN:  cam_theta += angle_step; if(cam_theta > 3.0) cam_theta = 3.0; break;
    }
    glutPostRedisplay();
}

void normalKeys(unsigned char key, int x, int y) {
    if (key == 'r') {
        cam_theta = 1.57; cam_phi = 0.0;
        glutPostRedisplay();
    }
}

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("GPU Black Hole Raytracer");
    
    initShader();
    
    glutDisplayFunc(display);
    glutSpecialFunc(processKeys);
    glutKeyboardFunc(normalKeys);
    
    printf("GPU Accelerated.\nControls:\n  Arrow Keys: Orbit\n  R: Reset\n");
    
    glutMainLoop();
    return 0;
}
