import pygame

pygame.init()
window = pygame.display.set_mode((400, 400))
clock = pygame.time.Clock()

# Define some points for a triangle
points = [(100, 100), (300, 100), (200, 300)]

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    window.fill((0, 0, 0))  # Clear the screen

    # 1. Filled Polygon (No Anti-Aliasing)
    pygame.draw.polygon(window, (255, 0, 0), points)  # Red, filled

    # 2. Anti-Aliased Polygon Outline
    pygame.draw.aalines(window, (0, 255, 0), True, points)  # Green, anti-aliased outline

    pygame.display.flip()
    clock.tick(60)

pygame.quit()
