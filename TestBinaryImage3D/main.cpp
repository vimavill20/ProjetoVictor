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
//#include "gmsh.h"
#include <fstream>
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZLinearAnalysis.h"
#include "TPZExtendGridDimension.h"
#include "pzvisualmatrix.h"


//CODIGO PARA OBTENER LOS BITS EN UN TXT DE UNA IMAGEN .TIF

int main (){
  
    std::ofstream outFile("Sample0.txt");
    std::cout << std::endl;
    cv::Mat img = cv::imread("/Users/victorvillegassalabarria/Downloads/Sample000.tif", cv::IMREAD_GRAYSCALE);
    
            if (img.empty()) {
                std::cout << "No se pudo abrir la imagen\n";
                return -1;
            }
            std::cout << "rows " << img.rows << " cols " << img.cols << std::endl;
            for (int i = 0; i < img.rows; ++i) {
                for (int j = 0; j < img.cols; ++j) {
                    // Convertir cada píxel a binario y imprimirlo
                    outFile << std::bitset<1>(img.at<uchar>(i, j))<< ' ';
    //                    if ((i * img.cols + j + 1) % 676 == 0) {
                                outFile << '\n';
    //                            }
    //                    std::cout << (int)img.at<uchar>(i, j) << ' ';
                }
                //std::cout << '\n';
            }

}


