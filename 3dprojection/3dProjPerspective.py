import pygame
import numpy as np
from math import sin, cos, tan
from dataclasses import dataclass
from typing import List # enable us to specify datatype in mesh.tris

"""
What this script does is render a 3D cube in a pygame window, and the user may rotate it by pressing the WASD keys.

Implementing a slightly different take on cube rendering, based on this tutorial: https://www.youtube.com/watch?v=ih20l3pJoeU
    - Modified the file, 3dproj.py

Current state
    - The cube mesh is implemented and can be rendered

Next steps
    - Change the perspective so the viewer isn't at an infinite distance.
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

# @dataclass
# class mat4x4:


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

def mesh_to_nparray(mesh):
    points = []
    for i in range(len(cube.tris)):
        for j in range(3):
            #print(cube.tris[i].p[j].x, cube.tris[i].p[j].y, cube.tris[i].p[j].z)
            points.append([cube.tris[i].p[j].x, cube.tris[i].p[j].y, cube.tris[i].p[j].z])
    return points

def vec3d_to_list(vec3d):
    return [vec3d.x, vec3d.y, vec3d.z]
    
def get_triangle_edges(mesh):
    edges = []
    for i in range(len(cube.tris)):
        for j in range(3):
            edges.append((vec3d_to_list(cube.tris[i].p[j]), vec3d_to_list(cube.tris[i].p[j-1])))
        #edges.append((i, j))
    return edges

print(mesh_to_nparray(cube))
print("-"*50)
print(get_triangle_edges(cube))



# Projection matrix for 3D to 2D projection
projection_matrix = np.array([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 0]])


def project_shape(pm, shape):
    """
    Project shape onto screen

    Steps:
        1. Transpose the projection matrix
        2. Multiply the shape matrix by the transposed projection matrix
        3. Transpose the resulting matrix
        4. Delete the 3rd (index 2) column (axis=1) of the resulting matrix, whih maps to the z-axis.
    """
    return np.delete(np.transpose(np.matmul(shape,np.transpose(pm))),2,axis=1)   

def rotate_shape(shape, theta_x, theta_y, theta_z):
    rot_x = np.array([[1, 0, 0],
            [0, np.cos(theta_x), -np.sin(theta_x)],
            [0, np.sin(theta_x), np.cos(theta_x)]])
    rot_y = np.array([[np.cos(theta_y), 0, np.sin(theta_y)],
            [0, 1, 0],
            [-np.sin(theta_y), 0, np.cos(theta_y)]])
    rot_z = np.array([[np.cos(theta_z), -np.sin(theta_z), 0],
            [np.sin(theta_z), np.cos(theta_z), 0],
            [0, 0, 1]])
    
    shape = np.matmul(shape, rot_x)
    shape = np.matmul(shape, rot_y)
    shape = np.matmul(shape, rot_z)
    return shape

def adjust(p, scale, offset):
    return p*scale + offset

def connect_points(p_0, p_1):
    pygame.draw.line(window, (255, 0, 0), p_0, p_1)

while True:
    clock.tick(TICK)
    window.fill((0, 0, 0))

    mesh_points = mesh_to_nparray(cube)
    mesh_edges = get_triangle_edges(cube)
    
    # Generate rotated and projected point pairs for each triangle edge
    rotated_edges = []
    for i in range(len(mesh_edges)):
        rotated_point_0 = rotate_shape(mesh_edges[i][0], angle_x, angle_y, angle_z)
        rotated_point_1 = rotate_shape(mesh_edges[i][1], angle_x, angle_y, angle_z)
        # print(rotated_point_0)

        projected_point_0 = project_shape([rotated_point_0], projection_matrix)
        projected_point_1 = project_shape([rotated_point_1], projection_matrix)
        
        rotated_edges.append((projected_point_0, projected_point_1))
    # print("Rotated_edges: " + str(rotated_edges))

    # Generate rotated and projected points for each corner of the cube
    # rotated_shape = rotate_shape(mesh_points, angle_x, angle_y, angle_z)
    # projection = project_shape(rotated_shape, projection_matrix)
    
    # for point in projection:
    #     x, y = adjust(point[0], scale, WINDOW_SIZE/2), adjust(point[1], scale, WINDOW_SIZE/2)
    #     pygame.draw.circle(window, (255, 0, 0), (x, y), 5)

    for edge in range(len(rotated_edges)):
        # print(rotated_edges[edge][0][0][1])
        p_0 = adjust(rotated_edges[edge][0][0][0], scale, WINDOW_SIZE/2), adjust(rotated_edges[edge][0][0][1], scale, WINDOW_SIZE/2)
        p_1 = adjust(rotated_edges[edge][1][0][0], scale, WINDOW_SIZE/2), adjust(rotated_edges[edge][1][0][1], scale, WINDOW_SIZE/2)
        connect_points(p_0, p_1)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            print("Length of points array is: " + str(len(points)))
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
