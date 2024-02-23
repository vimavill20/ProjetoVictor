# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)

project (MyProject)

add_executable(TestVictor main.cpp)


add_subdirectory(TIF3dVTK)

add_subdirectory(TestOpenCV)
add_subdirectory(TestBinaryImage)
add_subdirectory(Raw3dVTK)


