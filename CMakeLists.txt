cmake_minimum_required(VERSION 3.16)
project(LB2D)

set(CMAKE_CXX_STANDARD 20)

add_compile_definitions(USE_OPENGL_VISUALIZATION)

add_executable(LB2D src/main.cpp)
add_subdirectory(./third_party/EasyBMP)

find_package(GLEW REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)

if (MSVC)
    set_target_properties(
            LB2D PROPERTIES
            VS_DEBUGGER_WORKING_DIRECTORY "./"
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp /O2 /std:c++17 /arch:AVX2")
    target_link_libraries(LB2D PRIVATE GLUT::GLUT OpenGL::GLU OpenGL::GL GLEW::GLEW EasyBMP)
endif (MSVC)

if (UNIX)
    find_package(X11 REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -O3 -march=native")
    target_link_libraries(LB2D PRIVATE GLUT::GLUT OpenGL::GLU OpenGL::GL GLEW::GLEW X11 EasyBMP)
endif (UNIX)