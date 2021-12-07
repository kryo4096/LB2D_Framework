cmake_minimum_required(VERSION 3.21)
project(LB2D)

set(CMAKE_CXX_STANDARD 20)

add_compile_definitions(USE_OPENGL_VISUALIZATION)

add_executable(LB2D src/main.cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")

find_package(GLEW REQUIRED)
target_link_libraries(LB2D PRIVATE GLEW::GLEW)

find_package(GLUT REQUIRED)
target_link_libraries(LB2D PRIVATE GLUT::GLUT)

add_subdirectory(./third_party/EasyBMP)
target_link_libraries(LB2D PRIVATE EasyBMP)

