# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 3.11.0)


add_executable(TestVTK3D TestVTK3D.cpp TestVTK3D.h)



set(All_Dirs  NeoPZ::pz ${OpenCV_LIBS}  ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} )
#set(All_Dirs   ${OpenCV_LIBS} ${OpenCV_INCLUDE_DIRS} ${GMSH_INC} ${CMAKE_CURRENT_SOURCE_DIR} )


#add_executable(TestGmsh1 main1.cpp) # Cambia esta línea


#target_link_libraries(TestVictor PUBLIC NeoPZ::pz ${OpenCV_LIBS})

target_include_directories(TestVTK3D PUBLIC All_Dirs)
target_link_libraries(TestVTK3D PUBLIC NeoPZ::pz ${OpenCV_LIBS} ${GMSH_LIB})
#target_include_directories(TestVictor PUBLIC ${PZ_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
