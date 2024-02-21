#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <bitset>
#include "pzcmesh.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "pzmanvector.h"
#include "TPZGeoMeshTools.h"
#include "TPZCompMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzvec.h"
//#include "TPZEigenSolver"
#include <gmsh.h>
#include <fstream>
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "TPZDarcyFlow.h"
#include "TPZLinearAnalysis.h"
#include "TPZExtendGridDimension.h"
#include "pzvisualmatrix.h"

typedef unsigned char BYTE;
std::vector<BYTE> readFile(const char* filename)
{
    // open the file:
    std::streampos fileSize;
    std::ifstream file(filename, std::ios::binary);

    // get its size:
    file.seekg(0, std::ios::end);
    fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    // read the data:
    std::vector<BYTE> fileData(fileSize);
    file.read((char*) &fileData[0], fileSize);
    return fileData;
}
int main (){
//
    
    TPZVec<int> nx(2, 10);
    //What does it mean the line 46?
    
    nx[0]=5;
    nx[1]=5;
    const TPZVec<REAL> x0(3, 0.);
    const TPZVec<REAL> x1(3, 0.);
    x1[0]=1.0;
    x1[1]=1.0;
    auto msh = TPZGenGrid2D(nx, x0, x1);
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    msh.Read(gmesh);
    msh.SetElementType(MMeshType::EQuadrilateral);
    std::ofstream file2("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file2);
    
    REAL w= 0.2;
    TPZExtendGridDimension extend(gmesh, w);
    
    auto gmsh3D = extend.ExtendedMesh(5);
    std::ofstream file3("TestGeoMesh3D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmsh3D, file3);
    
    return 0;
}
