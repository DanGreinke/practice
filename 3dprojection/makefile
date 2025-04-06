CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c99 `pkg-config --cflags sdl2`
LDFLAGS = `pkg-config --libs sdl2`

3dRenderer: 3dRenderer.o
	$(CC) $(LDFLAGS) 3dRenderer.o -o 3dRenderer

3dRenderer.o: 3dRenderer.c
	$(CC) $(CFLAGS) -c 3dRenderer.c -o 3dRenderer.o

clean:
	rm -f 3dRenderer 3dRenderer.o

run: 3dRenderer
	./3dRenderer
