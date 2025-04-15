# --- Variables ---
#
# CC: The C compiler to use (gcc).
CC = gcc

# CFLAGS: Compiler flags.
#   -Wall: Enable all warnings.
#   -Wextra: Enable extra warnings.
#   -O2: Enable optimization level 2.
#   -std=c99: Use the C99 standard.
#   -lm: Link the math library.
CFLAGS = -Wall -Wextra -O2 -std=c99

# TARGET: The name of the executable file.
TARGET = mandelbrot

# SRC: The C source file.
SRC = mandelbrot.c

# OBJ: The name of the object file.
OBJ = mandelbrot.o

# --- Targets ---

# Default target: Build the executable.
all: $(TARGET)

# Link the object file to create the executable.
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ)

# Compile the source file to an object file.
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c $(SRC)

# Clean up object files and the executable.
clean:
	rm -f $(OBJ) $(TARGET)

# Run the executable
run: $(TARGET)
	./$(TARGET)