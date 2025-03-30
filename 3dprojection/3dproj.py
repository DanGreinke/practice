import pygame
import numpy as np
from math import sin, cos

WINDOW_SIZE = 800
window = pygame.display.set_mode((WINDOW_SIZE, WINDOW_SIZE))
clock = pygame.time.Clock()

# projection_matrix = [[1, 0, 0],
#                      [0, 1, 0],
#                      [0, 0, 0]]

# cube_points = [n for n in range(8)]
# cube_points[0] = [[-1], [-1], [-1]]
# cube_points[1] = [[1], [-1], [-1]]
# cube_points[2] = [[-1], [1], [-1]]
# cube_points[3] = [[-1], [-1], [1]]
# cube_points[4] = [[1], [1], [-1]]
# cube_points[5] = [[1], [-1], [1]]
# cube_points[6] = [[-1], [1], [1]]
# cube_points[7] = [[1], [1], [1]]

projection_matrix = np.array([[1, 0, 0],
                              [0, 1, 0],
                              [0, 0, 0]])

cube_points = np.array([[-1, -1, -1],
                        [-1, 1, -1],
                        [-1, 1, 1],
                        [-1, -1, 1],
                        [1, -1, 1],
                        [1, 1, 1],
                        [1, 1, -1],
                        [1, -1, -1]])

scale = 100
angle_x = angle_y = angle_z = 0

def multiply_m(mask, shape):
    """
    Multiply two matrices

    Steps:
        1. Transpose Matrix B -> the projection matrix
        2. Multiply Matrix A by the transposed Matrix B
        3. Transpose the resulting matrix
        4. Delete the 3rd (index 2) column (axis=1) of the resulting matrix, whih maps to the z-axis.
    """
    return np.delete(np.transpose(np.matmul(shape,np.transpose(mask))),2,axis=1)   

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
    pygame.draw.line(window, (255, 255, 255), (points[i][0], points[i][1]), (points[j][0], points[j][1]))

while True:
    clock.tick(60)
    window.fill((0, 0, 0))
    angle_x += 0.01
    angle_y += 0.01
    angle_z += 0.01

    rotated_shape = rotate_shape(cube_points, angle_x, angle_y, angle_z)
    projection = multiply_m(rotated_shape, projection_matrix)
    
    points = [0 for _ in range(len(cube_points))]
    i = 0
    
    for point in projection:
        adjust = lambda p, scale, offset: (p*scale + offset)
        x, y = adjust(point[0], scale, WINDOW_SIZE/2), adjust(point[1], scale, WINDOW_SIZE/2)
        pygame.draw.circle(window, (255, 0, 0), (x, y), 5)
        
        points[i] = (x, y)
        i += 1

    
    for i in range(len(points)):
        connect_points(i-1, i, points)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            print("Length of points array is: " + str(len(points)))
            pygame.quit()
            quit
        pygame.display.update()
