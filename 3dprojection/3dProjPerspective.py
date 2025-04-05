import pygame
import numpy as np
import math
from math import tan, sin, cos
from dataclasses import dataclass
from typing import List # enable us to specify datatype in mesh.tris
# import os
from pathlib import Path
from datetime import datetime
"""
What this script does is render a 3D cube in a pygame window, and the user may rotate it by pressing the WASD keys.

Implementing a slightly different take on cube rendering, based on this tutorial: https://www.youtube.com/watch?v=ih20l3pJoeU
    - Modified the file, 3dproj.py
    - See https://github.com/OneLoneCoder/Javidx9/blob/master/ConsoleGameEngine/BiggerProjects/Engine3D/OneLoneCoder_olcEngine3D_Part1.cpp

I subbed in dataclasses for the structs in the original c++ code, used numpy for matrix multiplication, and pygame for rendering.
In addition, I decided to hard-code the projection and rotation matrices, to make them easier to see.

Note: I decided to get Gemini to help vibecode a version of this in C, 3dRenderer.c, because this one hits 5 fps with the teapot.obj file.
"""


# Define data structures to support my "3D Engine"
@dataclass
class vec3d:
    x: float
    y: float
    z: float
    w: float = 1.0
    
@dataclass
class triangle:
    p: tuple[vec3d, vec3d, vec3d]
    c: vec3d = vec3d(0, 0, 0)

@dataclass
class mesh:
    tris: List[triangle]

# Declare global variables
obj_distance = 8
f_near = 0.1
f_far = 1000
f_fov = 90
f_fov_rad = 1 / tan(f_fov * 0.5 / 180.0 * math.pi)
viewer_Camera = vec3d(0, 0, 0)
light_direction = vec3d(0, 0, -1)
light_color = vec3d(255, 0, 0)
scale = 200
angle_x = angle_y = angle_z = 0
WINDOW_SIZE = 800
f_aspect_ratio = WINDOW_SIZE / WINDOW_SIZE
TICK = 60
ROTATE_SPEED = 1 / TICK
window = pygame.display.set_mode((WINDOW_SIZE, WINDOW_SIZE))
clock = pygame.time.Clock()

# Initialize projection matrix
projection_matrix = np.array([[f_aspect_ratio * f_fov_rad, 0,         0,                                    0],
                                [0,                          f_fov_rad, 0,                                    0],
                                [0,                          0,         f_far / (f_far - f_near),             1],
                                [0,                          0,         (-f_far * f_near) / (f_far - f_near), 0]
                            ])

# Normalize the light direction
l = math.sqrt(light_direction.x**2 + light_direction.y**2 + light_direction.z**2)
normalized_light_direction = vec3d(light_direction.x/l, light_direction.y/l, light_direction.z/l)

def load_obj(obj_file):
    p = Path(__file__).parent / obj_file
    vectors = []
    triangles = []
    out_tris = []
    with open(obj_file, 'r') as f:
        lines = f.readlines()
        for l in range(len(lines)):
            if lines[l].startswith('v'):
                vec = vec3d(float(lines[l].split()[1]), 
                            float(lines[l].split()[2]), 
                            float(lines[l].split()[3]))
                vectors.append(vec)
            elif lines[l].startswith('f'):
                tri_indices = (int(lines[l].split()[1]) - 1, 
                               int(lines[l].split()[2]) - 1, 
                               int(lines[l].split()[3]) - 1)
                triangles.append(tri_indices)
        for t in range(len(triangles)):
            i1,i2,i3 = triangles[t]
            tri = triangle((vectors[i1], vectors[i2], vectors[i3]))
            out_tris.append(tri)
    return mesh(out_tris)

cube = load_obj("VideoShip.obj") # "Cube"


def vec3d_mult_by_4x4_matrix(vec, m):
    """
    Using a 4x4 matrix to capture the perspective adjustments that are necessary
    to make the far side of the cube appear smaller than the near side. This 
    requires us to add a forth dimension to our vectors, and multiply them by the
    projection matrix.

    The rotation matrices have been similarly adjusted to make this function work
    for the rotation calculations as well.
    """
    tmp_vec = np.array([vec.x, vec.y, vec.z, vec.w])
    tmp_vec = np.matmul(tmp_vec, m)
    return vec3d(tmp_vec[0], tmp_vec[1], tmp_vec[2], tmp_vec[3])

def tri_div_by_w(tri):
    """
    Divide the x and y, and z coordinates by w to get the correct perspective
    """
    out_tri = []
    for p in tri.p:
        p = vec3d(p.x / p.w, p.y / p.w, p.z / p.w, p.w)
        out_tri.append(p)
    out_tri = triangle(tuple(out_tri))
    return out_tri

def project_triangle(tri):
    """
    Project shape onto screen"
    """
    # projection_matrix = np.array([[f_aspect_ratio * f_fov_rad, 0,         0,                                    0],
    #                               [0,                          f_fov_rad, 0,                                    0],
    #                               [0,                          0,         f_far / (f_far - f_near),             1],
    #                               [0,                          0,         (-f_far * f_near) / (f_far - f_near), 0]
    #                             ])
    points = np.array([[p.x, p.y, p.z, p.w] for p in tri.p])
    projected_points = np.matmul(points, projection_matrix)
    projected_tri = triangle(tuple(vec3d(*p) for p in projected_points))
    return projected_tri

def transform_triangle(tri, world_mat):
    # out_tri = []
    # for p in tri.p:
    #     p = vec3d_mult_by_4x4_matrix(p, world_mat)
    #     out_tri.append(p)
    # out_tri = triangle(tuple(out_tri))
    # return out_tri
    points = np.array([[p.x, p.y, p.z, p.w] for p in tri.p])
    transformed_points = np.matmul(points, world_mat)
    transformed_tri = triangle(tuple(vec3d(*p) for p in transformed_points))
    return transformed_tri


def translate_triangle(tri, translation_matrix):
    # translated_tri = []
    # for p in tri.p:
    #     p = vec3d_mult_by_4x4_matrix(p, translation_matrix)
    #     translated_tri.append(p)
    # translated_tri = triangle(tuple(translated_tri))
    # return translated_tri
    points = np.array([[p.x, p.y, p.z, p.w] for p in tri.p])
    translated_points = np.matmul(points, translation_matrix)
    translated_tri = triangle(tuple(vec3d(*p) for p in translated_points))
    return translated_tri

def draw_triangle(x1, y1, x2, y2, x3, y3):
    x1, y1 = adjust(x1, scale, WINDOW_SIZE/2), adjust(y1, scale, WINDOW_SIZE/2)
    x2, y2 = adjust(x2, scale, WINDOW_SIZE/2), adjust(y2, scale, WINDOW_SIZE/2)
    x3, y3 = adjust(x3, scale, WINDOW_SIZE/2), adjust(y3, scale, WINDOW_SIZE/2)
    pygame.draw.line(window, (255, 255, 0), (x1, y1), (x2, y2))
    pygame.draw.line(window, (255, 255, 0), (x2, y2), (x3, y3))
    pygame.draw.line(window, (255, 255, 0), (x3, y3), (x1, y1))

def draw_polygon(x1, y1, x2, y2, x3, y3, color):
    # print(color)
    r, g, b = abs(round(color.x)), abs(round(color.y)), abs(round(color.z))
    # print(r,g,b)
    x1, y1 = adjust(x1, scale, WINDOW_SIZE/2), adjust(y1, scale, WINDOW_SIZE/2)
    x2, y2 = adjust(x2, scale, WINDOW_SIZE/2), adjust(y2, scale, WINDOW_SIZE/2)
    x3, y3 = adjust(x3, scale, WINDOW_SIZE/2), adjust(y3, scale, WINDOW_SIZE/2)
    pygame.draw.polygon(window, (r, g, b), [(x1, y1), (x2, y2), (x3, y3)])

def get_surface_normal(tri):
    """
    Get the surface normal of input triangle
    """
    # line_1 = np.array([tri.p[1].x - tri.p[0].x, tri.p[1].y - tri.p[0].y, tri.p[1].z - tri.p[0].z])
    # line_2 = np.array([tri.p[2].x - tri.p[0].x, tri.p[2].y - tri.p[0].y, tri.p[2].z - tri.p[0].z])

    # cross_product = np.cross(line_1, line_2)
    # surface_normal = cross_product / np.linalg.norm(cross_product)
    # normal_vec = vec3d(surface_normal[0], surface_normal[1], surface_normal[2])
    # return normal_vec
    line_1 = np.array([tri.p[1].x - tri.p[0].x, tri.p[1].y - tri.p[0].y, tri.p[1].z - tri.p[0].z])
    line_2 = np.array([tri.p[2].x - tri.p[0].x, tri.p[2].y - tri.p[0].y, tri.p[2].z - tri.p[0].z])
    surface_normal = np.cross(line_1, line_2)
    surface_normal /= np.linalg.norm(surface_normal)
    return vec3d(*surface_normal)

def make_translation_matrix(x_offset, y_offset, z_offset):
    """
    Make a translation matrix for the given offsets
    """
    translation_matrix = np.array([[1, 0, 0, 0],
                                   [0, 1, 0, 0],
                                   [0, 0, 1, 0],
                                   [x_offset, y_offset, z_offset, 1]])
    return translation_matrix

def make_rotation_matrix(theta_x, theta_y, theta_z):
    """
    Make a rotation matrix for the given angles
    """
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
    
    # Combine the rotation matrices
    rotation_matrix = np.matmul(rot_x, rot_y)
    rotation_matrix = np.matmul(rotation_matrix, rot_z)
    return rotation_matrix

# def rotate_triangle(tri, rot_mat):
#     rotated_points = []
#     for p in tri.p:
#         p = vec3d_mult_by_4x4_matrix(p, rot_mat)
#         rotated_points.append(p)
#     return triangle(tuple(rotated_points))

def adjust(p, scale, offset):
    return p*scale + offset

def tri_depth(tri):
    """
    Get the average depth of a triangle
    """
    return (tri.p[0].z + tri.p[1].z + tri.p[2].z) / 3.0

def check_visibility(tri, cam=viewer_Camera):
    """
    When we look at the edge of an object, there's a roughly 90 degree angle from the surface normal
    our line of sight.

    Here we check the visibility of a surface by:
        1. Identify the surface normal for the triangle
        2. Find the vector from the 0th triangle vertex to the viewer (default at [0,0,0])
        3. Calculate the dot product b/w that surface normal vector and the line of sight vector
        4. Values less than 0 mean that the surface is visible, and values greater than 0 mean
           that the surface is facing away. 
    """
    # tri_normal = get_surface_normal(tri)
    # dot_prod = np.dot([tri_normal.x, tri_normal.y, tri_normal.z], 
    #                   [tri.p[0].x - cam.x, tri.p[0].y - cam.y, tri.p[0].z - cam.z]
    #                 )
    # return dot_prod < 0, tri_normal
    tri_normal = get_surface_normal(tri)
    tri_normal_np = np.array([tri_normal.x, tri_normal.y, tri_normal.z])
    cam_to_vertex = np.array([tri.p[0].x - cam.x, tri.p[0].y - cam.y, tri.p[0].z - cam.z])
    dot_prod = np.dot(tri_normal_np, cam_to_vertex)
    return dot_prod < 0, tri_normal

# def illuminate_triangle(light_vec, tri_normal):
#     l = math.sqrt(light_vec.x**2 + light_vec.y**2 + light_vec.z**2)
#     normalized_light_vec = vec3d(light_vec.x/l, light_vec.y/l, light_vec.z/l)

#     dp = np.dot([normalized_light_vec.x, normalized_light_vec.y, normalized_light_vec.z],
#                 [tri_normal.x, tri_normal.y, tri_normal.z])
#     return dp

def illuminate_triangle(normalized_light_vec, tri_normal):
    """
    Calculate the illumination of a triangle.
    """
    dp = np.dot([normalized_light_vec.x, normalized_light_vec.y, normalized_light_vec.z],
                [tri_normal.x, tri_normal.y, tri_normal.z])
    return dp

while True:
    tic = datetime.now()
    clock.tick(TICK)
    window.fill((0, 0, 0))
    angle_x += ROTATE_SPEED
    angle_y += ROTATE_SPEED * 0.3
    angle_z += ROTATE_SPEED * 0.5

    rotation_matrix = make_rotation_matrix(angle_x, angle_y, angle_z)
    translation_matrix = make_translation_matrix(0, 0, obj_distance)
    world_matrix = np.matmul(rotation_matrix, translation_matrix)
    # Multiply cube points by rotation and projection matrices, then draw the cube
    triangles_to_raster = []

    for tri in cube.tris:
        transformed_triangle = transform_triangle(tri, world_matrix)
        visible, tri_normal = check_visibility(transformed_triangle) 
        projected_triangle = project_triangle(transformed_triangle)
        projected_triangle = tri_div_by_w(projected_triangle)
        if visible:
            # Compute Illumination for triangles
            illumination = illuminate_triangle(light_direction, tri_normal) 
            projected_triangle.c = vec3d(light_color.x * illumination, light_color.y * illumination, light_color.z * illumination)
            triangles_to_raster.append(projected_triangle)

    triangles_to_raster.sort(key=tri_depth, reverse=True)

    for tri in triangles_to_raster:
        # Fill in triangles
        draw_polygon(tri.p[0].x, tri.p[0].y, 
                     tri.p[1].x, tri.p[1].y, 
                     tri.p[2].x, tri.p[2].y,
                     tri.c)
        # Draw triangle edges
        # draw_triangle(projected_triangle.p[0].x, 
        #             projected_triangle.p[0].y, 
        #             projected_triangle.p[1].x, 
        #             projected_triangle.p[1].y,
        #             projected_triangle.p[2].x,
        #             projected_triangle.p[2].y)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            #print("Length of points array is: " + str(len(points)))
            pygame.quit()
            quit
    # keys = pygame.key.get_pressed()
    # if keys[pygame.K_r]:
    #     angle_x = angle_y = angle_z = 0
    # if keys[pygame.K_a]:
    #     angle_y -= ROTATE_SPEED
    # if keys[pygame.K_d]:
    #     angle_y += ROTATE_SPEED
    # if keys[pygame.K_w]:
    #     angle_x += ROTATE_SPEED
    # if keys[pygame.K_s]:
    #     angle_x -= ROTATE_SPEED
    # if keys[pygame.K_q]:
    #     angle_z -= ROTATE_SPEED
    # if keys[pygame.K_e]:
    #     angle_z += ROTATE_SPEED
    pygame.display.update()
    toc = datetime.now()
    print("FPS: " + str(1/(toc - tic).total_seconds()), end="\r")
