# CMakeLists.txt for compiling pendulum.c
# ------------------
# Usage Instructions:
# 1. Create build directory
#    mkdir build
#    cd build
# 2. Generate build system, pointing at parent dir
#    cmake ..
# 3. build
#    make
# 4. Run
#    ./pendulum
# 5. Clean
#    make clean
# ------------------
# Minimum required CMake version
cmake_minimum_required(VERSION 3.10) # You can adjust this to your required version

# Project name
project(pendulum_project)

# C Language standard
set(CMAKE_C_STANDARD 99)

# Add compiler flags
add_compile_options(-Wall -Wextra -O2)

# Add executable
add_executable(pendulum pendulum.c)

# Add math library
target_link_libraries(pendulum m)