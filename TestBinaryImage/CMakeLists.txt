# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)


add_executable(TestBinaryImage main.cpp)


set(All_Dirs  NeoPZ::pz ${OpenCV_LIBS}  ${OpenCV_INCLUDE_DIRS} ${GMSH_INC} ${CMAKE_CURRENT_SOURCE_DIR} )
#set(All_Dirs   ${OpenCV_LIBS} ${OpenCV_INCLUDE_DIRS} ${GMSH_INC} ${CMAKE_CURRENT_SOURCE_DIR} )


target_include_directories(TestBinaryImage PUBLIC All_Dirs)
target_link_libraries(TestBinaryImage PUBLIC NeoPZ::pz ${OpenCV_LIBS} ${GMSH_LIB})
#target_include_directories(TestBinaryImag PUBLIC ${PZ_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
