# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)


add_executable(VisualMatVtkPZ main.cpp)

# ----- Find PZ. By default searches one directory up -----
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)

find_package(OpenCV REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../opencv_install/ ${CMAKE_SOURCE_DIR}/opencv_install/)

#find_package(Gmsh REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../gmsh_install/ ${CMAKE_SOURCE_DIR}/gmsh_install/)

file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
#add_custom_target( GlobalFiles SOURCES ${PZ_SRC})


# linking with GMSH
find_library(GMSH_LIB gmsh HINTS "/home/victor/gmsh-install/lib" "../lib/gmsh-install/lib" "../gmsh-4.10.5-MacOSX-sdk/lib" "../gmsh_install/lib" NO_DEFAULT_PATH)
find_path(GMSH_INC gmsh.h HINTS "/home/victor//gmsh-install/include" "../lib/gmsh-install/include" "../gmsh-4.10.5-MacOSX-sdk/include" "../gmsh_install/include" NO_DEFAULT_PATH)
link_libraries(${GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
include_directories(${GMSH_INC})


set(All_Dirs  NeoPZ::pz ${OpenCV_LIBS}  ${OpenCV_INCLUDE_DIRS} ${GMSH_INC} ${CMAKE_CURRENT_SOURCE_DIR} )
#set(All_Dirs   ${OpenCV_LIBS} ${OpenCV_INCLUDE_DIRS} ${GMSH_INC} ${CMAKE_CURRENT_SOURCE_DIR} )


#add_executable(TestGmsh1 main1.cpp) # Cambia esta línea


#target_link_libraries(TestVictor PUBLIC NeoPZ::pz ${OpenCV_LIBS})

target_include_directories(VisualMatVtkPZ  PUBLIC All_Dirs)
target_link_libraries(VisualMatVtkPZ  PUBLIC NeoPZ::pz ${OpenCV_LIBS} ${GMSH_LIB})
#target_include_directories(TestBinaryImag PUBLIC ${PZ_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
