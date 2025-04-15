# Makefile for compiling pendulum.c
#
# This Makefile provides a simple way to compile and run the pendulum.c program.
# It also includes a clean target to remove intermediate and executable files.
#
# --- Explanation of Targets ---
#
# all: This is the default target. When you run 'make' without any arguments,
#      this target will be executed. It depends on the $(TARGET) executable, 
#      which means it will compile and link the code if necessary.
#
# $(TARGET): This target builds the executable file (pendulum) from the 
#            object file (pendulum.o). It depends on the $(OBJ) target, so 
#            the object file will be created first if it doesn't exist or if 
#            the source file has been modified.
#            The command executed is:
#            $(CC) $(CFLAGS) -o $(TARGET) $(OBJ)
#            - $(CC): The C compiler (gcc in this case).
#            - $(CFLAGS): The compiler flags (e.g., warnings, optimization, math library).
#            - -o $(TARGET): Specifies the output executable file name.
#            - $(OBJ): The object file to link.
#
# $(OBJ): This target compiles the source file (pendulum.c) into an object file
#         (pendulum.o). It depends on the $(SRC) file, so it will only be executed
#         if the source file is newer than the object file or if the object file
#         doesn't exist.
#         The command executed is:
#         $(CC) $(CFLAGS) -c $(SRC)
#         - $(CC): The C compiler.
#         - $(CFLAGS): The compiler flags.
#         - -c: This flag tells the compiler to compile the source file but not
#               link it, resulting in an object file.
#         - $(SRC): The C source file.
#
# clean: This target removes the object file and the executable file.
#        It's useful for starting a fresh build or cleaning up after testing.
#        The command executed is:
#        rm -f $(OBJ) $(TARGET)
#        - rm -f: The command to remove files forcefully.
#        - $(OBJ): The object file to remove.
#        - $(TARGET): The executable file to remove.
#
# run: This target first builds the executable (if needed) and then runs it.
#      The command executed is:
#      ./$(TARGET)
#      - ./$(TARGET): Executes the pendulum executable in the current directory.
#
# --- Instructions ---
# 1. Save this file as 'makefile' (or 'Makefile') in the same directory as 'pendulum.c'.
# 2. Open a terminal and navigate to that directory.
# 3. To compile and link the code, run:
#      make
# 4. To compile, link, and then run the code, run:
#      make run
# 5. To remove the compiled files (object and executable), run:
#      make clean
#
# -----------------
#
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
CFLAGS = -Wall -Wextra -O2 -std=c99 -lm

# TARGET: The name of the executable file.
TARGET = pendulum

# SRC: The C source file.
SRC = pendulum.c

# OBJ: The name of the object file.
OBJ = pendulum.o

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
