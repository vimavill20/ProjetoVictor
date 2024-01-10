# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)

project (MyProject)

add_executable(TestVictor main.cpp)


# ----- Find PZ. By default searches one directory up -----
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)

find_package(OpenCV REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../opencv_install/ ${CMAKE_SOURCE_DIR}/opencv_install/)

find_package(Gmsh REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../gmsh_install/ ${CMAKE_SOURCE_DIR}/gmsh_install/)

file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})


#target_link_libraries(TestVictor PUBLIC NeoPZ::pz ${OpenCV_LIBS})
target_link_libraries(TestVictor PUBLIC NeoPZ::pz ${OpenCV_LIBS} ${GMSH_LIBS})


#target_include_directories(TestVictor PUBLIC ${PZ_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
target_include_directories(TestVictor PUBLIC ${PZ_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} ${Gmsh_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
