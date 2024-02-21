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
#include "pzvisualmatrix.h"

//#include "TPZGeoMeshTools.h"
//#include "TPZCompMeshTools.h"
//#include "TPZRefPattern.h"
//#include "TPZGenGrid2D.h"
//#include "TPZVTKGeoMesh.h"
//#include "pzvec.h"
////#include "TPZEigenSolver"
//#include <gmsh.h>
//#include <fstream>
//#include "pzstepsolver.h"
//#include "TPZMaterial.h"
//#include "TPZDarcyFlow.h"
//#include "TPZLinearAnalysis.h"
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
// Función para crear una matriz binaria a partir de una imagen en escala de grises
TPZFMatrix<double> create_binary_matrix(const cv::Mat& img) {
    // Obtener las dimensiones de la imagen
    int rows = 676;//img.rows 676;
    int cols = 616;//img.cols 616;

    // Crear una matriz de 0 y 1 con el mismo tamaño que la imagen
    TPZFMatrix<double> matrix(rows,cols);

    // Recorrer la imagen y asignar el valor binario a la matriz
   
    for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
        REAL valtom = std::bitset<1>((img.at<uchar>(row, col) )).to_ulong();
//            matrix(rows - row -1, col) =valtom;
        matrix(rows-row-1 ,col) =valtom;
        }
    }
    
    std::ofstream out("MatrixTest.txt");
    matrix.Print(out);
    std::string outfilename2 = "outputMatrixPZ2.vtk"; // Nombre del archivo de salida
    VisualMatrixVTK(matrix, outfilename2); //
    // Devolver la matriz binaria
    
 //   std::cout<<matrix<<std::endl;
    return matrix;
}

void VisualMatrixVTK(std::vector<TPZFMatrix<double>> & matrixVector, const std::string &outfilename)
{
    for(int k = 0; k < matrixVector.size(); k++)
    {
        TPZFMatrix<double> & matrix = matrixVector[k];
        const int nelx = matrix.Cols();
        const int nely = matrix.Rows();
        const int neltotal = nelx * nely;
        int i,j;
        ofstream out(outfilename + std::to_string(k) + ".vtk");
        out << "# vtk DataFile Version 3.0\n";
        out << "Generated by PZ\n";
        out << "ASCII\n";
        out << "DATASET RECTILINEAR_GRID\n";
        out << "DIMENSIONS " << (nelx+1) << " " <<  (nely+1) << " 1\n";
        out << "X_COORDINATES " << nelx+1 << " float\n";
        for (i=0; i<=nelx; i++) {
            out << i << " ";
        }
        out << std::endl;
        out << "Y_COORDINATES " << nely+1 << " float\n";
        for (j=0; j<=nely; j++) {
            out << j << " ";
        }
        out << std::endl;
        out << "Z_COORDINATES " << //nelz+1 << " float\n0.\n";
        out << "CELL_DATA " << nelx*nely << std::endl;
        out << "SCALARS mat_value float 1\n";
        out << "LOOKUP_TABLE default\n";

        for (j=0; j<nely; j++) {
            for (i=0;i<nelx;i++) {
                out << matrix[k](j,i) << std::endl;
            }
        }
        out.close();
    }
}


int main (){
//
    std::ifstream file("/Users/victorvillegassalabarria/Downloads/Sample_Labels_3D_RAW1.raw", std::ios::binary | std::ios::ate);
    if (file.is_open()) {
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);

        std::vector<char> buffer(size);
        file.read(buffer.data(), size);
//        if(!result) DebugStop();
        int count = 0;
        for(char c : buffer) {
//                for(int i = 7; i >= 0; --i) {
//                    std::cout << ((c >> i) & 1);
//                }
            std::cout << (int) c;
            std::cout << ' ';
            count ++;
            if(count > 50) break;
        }
            /* The entire file is now loaded into the buffer vector */
    } else {
        std::cout << "Unable to open file";
    }
    cv::Mat img = cv::imread("/Users/victorvillegassalabarria/Downloads/Sample000.tif", cv::IMREAD_GRAYSCALE);

    // Comprobar si se pudo abrir la imagen
    if (img.empty()) {
        std::cout << "No se pudo abrir la imagen\n";
        return -1;
    }

    // Crear una matriz binaria a partir de la imagen
    TPZFMatrix<double> matrix = create_binary_matrix(img);
    
  
    return 0;
}
//CODIGO PARA OBTENER EN UN TXT COORDENADAS DEL .TIF
//std::ofstream outFile("Sample0.txt");
//std::cout << std::endl;
//cv::Mat img = cv::imread("/Users/victorvillegassalabarria/Downloads/Sample000.tif", cv::IMREAD_GRAYSCALE);
//
//        if (img.empty()) {
//            std::cout << "No se pudo abrir la imagen\n";
//            return -1;
//        }
//        std::cout << "rows " << img.rows << " cols " << img.cols << std::endl;
//        for (int i = 0; i < img.rows; ++i) {
//            for (int j = 0; j < img.cols; ++j) {
//                // Convertir cada píxel a binario y imprimirlo
//                outFile << std::bitset<1>(img.at<uchar>(i, j))<< ' ';
////                    if ((i * img.cols + j + 1) % 676 == 0) {
//                            outFile << '\n';
////                            }
////                    std::cout << (int)img.at<uchar>(i, j) << ' ';
//            }
//            //std::cout << '\n';
//        }

