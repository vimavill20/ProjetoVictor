# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)

project (MyProject)

# linking with GMSH
find_library(GMSH_LIB gmsh HINTS "/home/pedro/projetos/library/gmsh-install/lib" "../lib/gmsh-install/lib" "../gmsh-4.10.5-MacOSX-sdk/lib" "../gmsh-4.12.2-MacOSARM-sdk/lib" "../gmsh_install/lib" NO_DEFAULT_PATH)
find_path(GMSH_INC gmsh.h HINTS "/home/pedro/projetos/library/gmsh-install/include" "../lib/gmsh-install/include" "../gmsh-4.10.5-MacOSX-sdk/include" "../gmsh-4.12.2-MacOSARM-sdk/include" "../gmsh_install/include" NO_DEFAULT_PATH)

link_libraries(${GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
# include_directories(${GMSH_INC})

find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)

find_package(OpenCV REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../opencv_install/ ${CMAKE_SOURCE_DIR}/opencv_install/)

add_subdirectory(TIF3dVTK)

add_subdirectory(TestNeoPz)

#add_subdirectory(TestOpenCV)
#add_subdirectory(TestBinaryImage)
add_subdirectory(TestBinaryImage3D)


