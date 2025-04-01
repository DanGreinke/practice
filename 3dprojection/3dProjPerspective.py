import pygame
import numpy as np
import math
from math import tan, sin, cos
from dataclasses import dataclass
from typing import List # enable us to specify datatype in mesh.tris

"""
What this script does is render a 3D cube in a pygame window, and the user may rotate it by pressing the WASD keys.

Implementing a slightly different take on cube rendering, based on this tutorial: https://www.youtube.com/watch?v=ih20l3pJoeU
    - Modified the file, 3dproj.py

I subbed in dataclasses for the structs in the original c++ code, used numpy for matrix multiplication, and pygame for rendering.
In addition, I decided to hard-code the projection and rotation matrices, to make them easier to see.
"""
scale = 100
angle_x = angle_y = angle_z = 0
WINDOW_SIZE = 800
TICK = 60
ROTATE_SPEED = 1 / TICK
window = pygame.display.set_mode((WINDOW_SIZE, WINDOW_SIZE))
clock = pygame.time.Clock()

@dataclass
class vec3d:
    x: float
    y: float
    z: float
    
@dataclass
class triangle:
    p: tuple[vec3d, vec3d, vec3d]

@dataclass
class mesh:
    tris: List[triangle]

cube = mesh([
    # South
    triangle([vec3d(-1.0, -1.0, -1.0), vec3d(-1.0, 1.0, -1.0), vec3d(1.0, 1.0, -1.0)]),
    triangle([vec3d(-1.0, -1.0, -1.0), vec3d(1.0, 1.0, -1.0), vec3d(1.0, -1.0, -1.0)]),
    
    # East
    triangle([vec3d(1.0, -1.0, -1.0), vec3d(1.0, 1.0, -1.0), vec3d(1.0, 1.0, 1.0)]),
    triangle([vec3d(1.0, -1.0, -1.0), vec3d(1.0, 1.0, 1.0), vec3d(1.0, -1.0, 1.0)]),

    # North
    triangle([vec3d(1.0, -1.0, 1.0), vec3d(1.0, 1.0, 1.0), vec3d(-1.0, -1.0, 1.0)]),
    triangle([vec3d(1.0, 1.0, 1.0), vec3d(-1.0, 1.0, 1.0), vec3d(-1.0, -1.0, 1.0)]),

    # West
    triangle([vec3d(-1.0, -1.0, 1.0), vec3d(-1.0, 1.0, 1.0), vec3d(-1.0, 1.0, -1.0)]),
    triangle([vec3d(-1.0, -1.0, 1.0), vec3d(-1.0, 1.0, -1.0), vec3d(-1.0, -1.0, -1.0)]),

    # Top
    triangle([vec3d(-1.0, 1.0, -1.0), vec3d(-1.0, 1.0, 1.0), vec3d(1.0, 1.0, 1.0)]),
    triangle([vec3d(-1.0, 1.0, -1.0), vec3d(1.0, 1.0, 1.0), vec3d(1.0, 1.0, -1.0)]),

    # Bottom
    triangle([vec3d(1.0, -1.0, 1.0), vec3d(-1.0, -1.0, 1.0), vec3d(-1.0, -1.0, -1.0)]),
    triangle([vec3d(1.0, -1.0, 1.0), vec3d(-1.0, -1.0, -1.0), vec3d(1.0, -1.0, -1.0)])
    ])

# Projection Matrix
f_near = 0.1
f_far = 1000
f_fov = 90
f_aspect_ratio = WINDOW_SIZE / WINDOW_SIZE
f_fov_rad = 1 / tan(f_fov * 0.5 / 180.0 * math.pi)

def vec3d_mult_by_4x4_matrix(vec, m):
    """
    Using a 4x4 matrix to capture the perspective adjustments that are necessary
    to make the far side of the cube appear smaller than the near side. This 
    requires us to add a forth dimension to our vectors, and multiply them by the
    projection matrix.

    The rotation matrices have been similarly adjusted to make this function work
    for the rotation calculations as well.
    """
    tmp_vec = np.array([vec.x, vec.y, vec.z, 1])
    tmp_vec = np.matmul(tmp_vec, m)

    if tmp_vec[3] != 0:
        tmp_vec /= tmp_vec[3]
    return vec3d(tmp_vec[0], tmp_vec[1], tmp_vec[2])

def project_triangle(tri):
    """
    Project shape onto screen"
    """
    projection_matrix = np.array([[f_aspect_ratio * f_fov_rad, 0,         0,                                    0],
                                  [0,                          f_fov_rad, 0,                                    0],
                                  [0,                          0,         f_far / (f_far - f_near),             1],
                                  [0,                          0,         (-f_far * f_near) / (f_far - f_near), 0]
                                ])
    out_tri = []
    for p in tri.p:
        p = vec3d_mult_by_4x4_matrix(p, projection_matrix)
        out_tri.append(p)
    out_tri = triangle(tuple(out_tri))
    return out_tri

def translate_triangle(tri, x_offset, y_offset, z_offset):
    translated_tri = triangle([vec3d(tri.p[0].x + x_offset, tri.p[0].y + y_offset, tri.p[0].z + z_offset),
                               vec3d(tri.p[1].x + x_offset, tri.p[1].y + y_offset, tri.p[1].z + z_offset),
                               vec3d(tri.p[2].x + x_offset, tri.p[2].y + y_offset, tri.p[2].z + z_offset)])
    return translated_tri

def draw_triangle(x1, y1, x2, y2, x3, y3):
    x1, y1 = adjust(x1, scale, WINDOW_SIZE/2), adjust(y1, scale, WINDOW_SIZE/2)
    x2, y2 = adjust(x2, scale, WINDOW_SIZE/2), adjust(y2, scale, WINDOW_SIZE/2)
    x3, y3 = adjust(x3, scale, WINDOW_SIZE/2), adjust(y3, scale, WINDOW_SIZE/2)
    pygame.draw.line(window, (255, 0, 0), (x1, y1), (x2, y2))
    pygame.draw.line(window, (255, 0, 0), (x2, y2), (x3, y3))
    pygame.draw.line(window, (255, 0, 0), (x3, y3), (x1, y1))


def rotate_triangle(tri, theta_x, theta_y, theta_z):
    #print(tri)
    rot_x = np.array([[1, 0,                0,               0],
                      [0, np.cos(theta_x), -np.sin(theta_x), 0],
                      [0, np.sin(theta_x),  np.cos(theta_x), 0],
                      [0, 0,                0,               1]])
    
    rot_y = np.array([[ np.cos(theta_y), 0, np.sin(theta_y), 0],
                      [ 0,               1, 0,               0],
                      [-np.sin(theta_y), 0, np.cos(theta_y), 0],
                      [ 0,               0, 0,               1]])
    
    rot_z = np.array([[np.cos(theta_z), -np.sin(theta_z), 0, 0],
                      [np.sin(theta_z),  np.cos(theta_z), 0, 0],
                      [0,                0,               1, 0],
                      [0,                0,               0, 1]])
    rotated_points = []
    for p in tri.p:
        p = vec3d_mult_by_4x4_matrix(p, rot_x)
        p = vec3d_mult_by_4x4_matrix(p, rot_y)
        p = vec3d_mult_by_4x4_matrix(p, rot_z)
        rotated_points.append(p)
    return triangle(tuple(rotated_points))

def adjust(p, scale, offset):
    return p*scale + offset

while True:
    clock.tick(TICK)
    window.fill((0, 0, 0))

    # Multiply cube points by rotation and projection matrices, then draw the cube
    for tri in cube.tris:
        rotated_triangle = rotate_triangle(tri, angle_x, angle_y, angle_z)
        translated_triangle = translate_triangle(rotated_triangle, 0, 0, 3)
        projected_triangle = project_triangle(translated_triangle)
        draw_triangle(projected_triangle.p[0].x, 
                      projected_triangle.p[0].y, 
                      projected_triangle.p[1].x, 
                      projected_triangle.p[1].y,
                      projected_triangle.p[2].x,
                      projected_triangle.p[2].y)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            #print("Length of points array is: " + str(len(points)))
            pygame.quit()
            quit
    keys = pygame.key.get_pressed()
    if keys[pygame.K_r]:
        angle_x = angle_y = angle_z = 0
    if keys[pygame.K_a]:
        angle_y -= ROTATE_SPEED
    if keys[pygame.K_d]:
        angle_y += ROTATE_SPEED
    if keys[pygame.K_w]:
        angle_x += ROTATE_SPEED
    if keys[pygame.K_s]:
        angle_x -= ROTATE_SPEED
    if keys[pygame.K_q]:
        angle_z -= ROTATE_SPEED
    if keys[pygame.K_e]:
        angle_z += ROTATE_SPEED
    pygame.display.update()
