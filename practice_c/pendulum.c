#include <stdio.h>
#include <math.h>

#define PI 3.14159265
#define G 9.8
#define L 2.0
#define mu 0.1

#define THETA_DOT_0 0.0
#define THETA_0 PI / 3.0

double theta_dot_dot(double theta, double theta_dot){
    return -mu * theta_dot - (G / L) * sin(theta);
}

int main(void){
    double t;
    double theta_double_dot;

    double theta = THETA_0;
    double theta_dot = THETA_DOT_0;
    double delta_t = 0.01;
    
    for(t=0.0; t<=2.0; t+=delta_t){
        theta_double_dot = theta_dot_dot(theta, theta_dot);
        theta += theta_dot * delta_t;
        theta_dot += theta_double_dot * delta_t;
    }

    printf("theta = %f\n", theta);
    return 0;
}






