# @file ProjetoVictor/CMakeList.txt  -- First file to read for CMake 



add_executable(TestBinaryImage3D main.cpp)


set(All_Dirs  NeoPZ::pz  ${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} )

target_include_directories(TestBinaryImage3D PUBLIC All_Dirs)
target_link_libraries(TestBinaryImage3D PUBLIC NeoPZ::pz ${OpenCV_LIBS})
