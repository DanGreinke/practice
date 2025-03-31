import pygame
import numpy as np
from math import sin, cos

"""
What this script does is render a 3D cube in a pygame window, and the user may rotate it by pressing the WASD keys.

I loosely followed along with this youtube tutorial: https://www.youtube.com/watch?v=sQDFydEtBLE

Where the tutorial involved a custom function for matrix operations, I used numpy instead.
"""
scale = 100
angle_x = angle_y = angle_z = 0
WINDOW_SIZE = 800
TICK = 60
ROTATE_SPEED = 1 / TICK
window = pygame.display.set_mode((WINDOW_SIZE, WINDOW_SIZE))
clock = pygame.time.Clock()

# Projection matrix for 3D to 2D projection
projection_matrix = np.array([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 0]])

# Define cube points and edges
cube_points = np.array([[-1, -1, -1],
                        [-1, 1, -1],
                        [-1, 1, 1],
                        [-1, -1, 1],
                        [1, -1, 1],
                        [1, 1, 1],
                        [1, 1, -1],
                        [1, -1, -1]])

cube_edges = [
    (0, 1), (1, 2), (2, 3), (3, 0),  # Back face
    (4, 5), (5, 6), (6, 7), (7, 4),  # Front face
    (0, 7), (1, 6), (2, 5), (3, 4)   # Connecting edges
]

# Define 3D object class
class Object_3D:
    def __init__(self, points, edges):
        self.points = points
        self.edges = edges
    def points(self):
        return self.points

    def edges(self):
        return self.edges

# Instantiate cube object
cube = Object_3D(cube_points, cube_edges)

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

def connect_points(i, j, points):
    pygame.draw.line(window, (255, 0, 0), (points[i][0], points[i][1]), (points[j][0], points[j][1]))

while True:
    clock.tick(TICK)
    window.fill((0, 0, 0))

    rotated_shape = rotate_shape(cube.points, angle_x, angle_y, angle_z)
    projection = project_shape(rotated_shape, projection_matrix)
    
    points = [0 for _ in range(len(cube.points))]
    i = 0
    
    for point in projection:
        adjust = lambda p, scale, offset: (p*scale + offset)
        x, y = adjust(point[0], scale, WINDOW_SIZE/2), adjust(point[1], scale, WINDOW_SIZE/2)
        #pygame.draw.circle(window, (255, 0, 0), (x, y), 5)
        
        points[i] = (x, y)
        i += 1

    for edge in cube_edges:
        connect_points(edge[0], edge[1], points)

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
