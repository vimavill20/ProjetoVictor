#ifndef TestVTK3D_H
#define TestVTK3D_H
#include <iostream>
#include <math.h>
#include "pzcmesh.h"
//#include <opencv2/opencv.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
#include "pzmanvector.h"
#include "TPZGeoMeshTools.h"
#include "TPZCompMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
#include "pzvec.h"
//#include "pzvisualmatrix.h"
//#include "TPZEigenSolver"
#include <gmsh.h>
#include <fstream>
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "TPZDarcyFlow.h"
#include "TPZLinearAnalysis.h"

class TestVTK3D
{
public:


//Function to reconstruct 3d image on .vtk format using a TPZ vector of matrixs
    void VisualMatrix3DVTK(TPZVec<TPZFMatrix<double>> & matrixVector, const std::string &outfilename);

TPZFMatrix<double> create_binary_matrix(const cv::Mat& img);
//
//    struct AllSimulationData;
TPZFMatrix<double> Case10TIF(struct AllSimulationData &alldata, int &idata);
};

#endif // ARCHIVO1_H
