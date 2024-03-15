
#include <iostream>
#include <math.h>
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
//#include "pzvisualmatrix.h"
//#include "TPZEigenSolver"
#include <fstream>
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZLinearAnalysis.h"
#include "TestVTK3D.h"

//using std::cout;
//using std::endl;
//using std::cin;
using namespace cv;
using namespace std;
//Function to reconstruct 3d image on .vtk format using a TPZ vector of matrixs
void VisualMatrix3DVTK(TPZVec<TPZFMatrix<double>> & matrixVector, const std::string &outfilename)
{
    const int nelz=matrixVector.size();
    TPZFMatrix<double> & matrix1 = matrixVector[1];
    const int nelx = matrix1.Cols();
    const int nely = matrix1.Rows();
    const int neltotal = nelx * nely;
    int i,j;
    
    ofstream out(outfilename.c_str());
    out << "# vtk DataFile Version 3.0\n";
    out << "Generated by PZ\n";
    out << "ASCII\n";
    out << "DATASET RECTILINEAR_GRID\n";
    out << "DIMENSIONS " << (nelx+1) << " " <<  (nely+1) << " " << (nelz+1) <<"\n";
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
    out << "Z_COORDINATES " << nelz+1 << " float\n";
    for (int k=0; k<=nelz; k++) {
        out << k << " ";
    }
    out << std::endl;
    out << "CELL_DATA " << nelx*nely*nelz << std::endl;
    out << "SCALARS mat_value float 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int k = 0; k < matrixVector.size(); k++)
    {
        
        TPZFMatrix<double> & matrix = matrixVector[k];
        for (j=0; j<nely; j++) {
            for (i=0;i<nelx;i++) {
                out << matrix(j,i) << std::endl;
                }
        }float avance=static_cast<float>(k+1)/nelz;
        std::cout<<"Avance: "<<(avance)*100<<"%"<<std::endl;
        }
}
TPZFMatrix<double> create_binary_matrix(const cv::Mat& img) {
    // Obtener las dimensiones de la imagen
    int rows = 676;//img.rows 676;
    int cols = 616;//img.cols 616;
    // Crea una matriz de 0s y 1s con el mismo tamaño que la imagen
    TPZFMatrix<double> matrix(rows,cols);
    // Recorrer la imagen y asignar el valor binario a la matriz
    for (int row = 0; row < rows; row++) {
    for (int col = 0; col < cols; col++) {
        REAL valtom = std::bitset<1>((img.at<uchar>(row, col) )).to_ulong();
//            matrix(rows - row -1, col) =valtom;
        matrix(rows-row-1 ,col) =valtom;
        }
    }
    //std::ofstream out("MatrixTest1.txt");
    //matrix.Print(out);
    //std::string outfilename2 = "outputMatrixPZ2.vtk"; // Nombre del archivo de salida
    // Devolver la matriz binaria
    
 //   std::cout<<matrix<<std::endl;
    return matrix;
}

//

struct AllSimulationData
{
  
    std::vector<std::string> TIFnames;

};
TPZFMatrix<double> Case10TIF(AllSimulationData alldata, int idata){
    std::string filename=alldata.TIFnames[idata];
    filename=filename;
    cv::Mat img=cv::imread(filename, cv::IMREAD_GRAYSCALE);
    TPZFMatrix<double> matrix=create_binary_matrix(img);
    return matrix;
}
//
//MAIN ORIGINAL
//
int mainf (){
    AllSimulationData test;
    
    std::string common_name="/Users/philippedevloo/GitHub/CoreSampleResearch/ProjetoVictor/DADOS_TIF_10Layers/";

    std::vector<std::string> TIFnames;
    int nsim=10;//only to visualize the first 10 layers of the rocks
    for (int i = 0; i < nsim; ++i) {
        std::ostringstream oss;
        oss << common_name << "Sample" << std::setw(3) <<std::setfill('0') << i << ".tif";
        test.TIFnames.push_back(oss.str());
            }
    for (const auto & name : TIFnames) {
            std::cout << name << std::endl;
        }
    TPZVec<TPZFMatrix<double>>Vecmat10(nsim);
    for(int i=0;i<nsim;i++){
        auto matrizi= Case10TIF(test,i);
        Vecmat10[i]=matrizi;

    }
    std::ostringstream filenameStream;
    filenameStream << "Matrix3D" << nsim << ".vtk";
    std::string outfilename = filenameStream.str(); // Nombre del archivo de salida
    VisualMatrix3DVTK(Vecmat10, outfilename);
    return 0;
}


void generarArchivoGeo(const std::string& nombreArchivo, const std::vector<std::vector<int>>& mascaraBinaria) {
    std::ofstream archivoGeo(nombreArchivo);

    // Escribir la cabecera del archivo .geo
    archivoGeo << "lc = 0.1;\n";
    archivoGeo << "Point(1) = {0, 0, 0, lc};\n";  // Punto de origen

    // Definir puntos y líneas basadas en la máscara binaria
    int puntoID = 2;  // Comenzar desde 2 para evitar conflictos con el punto de origen

    for (int i = 0; i < mascaraBinaria.size(); ++i) {
        for (int j = 0; j < mascaraBinaria[i].size(); ++j) {
            if (mascaraBinaria[i][j] != 0) {
                // Crear un punto por cada píxel del objeto
                archivoGeo << "Point(" << puntoID << ") = {" << i << ", " << j << ", 0, lc};\n";
                ++puntoID;
            }
        }
    }

    // Crear líneas conectando los puntos en la misma fila
    for (int i = 0; i < mascaraBinaria.size(); ++i) {
        archivoGeo << "Spline(1) = {";
        for (int j = 0; j < mascaraBinaria[i].size(); ++j) {
            if (mascaraBinaria[i][j] != 0) {
                archivoGeo << puntoID - 1 << ", ";
            }
        }
        archivoGeo << puntoID - 1 << "};\n";
    }

    // Escribir el final del archivo .geo
    archivoGeo << "Transfinite Line {1} = " << puntoID - 1 << " Using Progression 1;\n";
    archivoGeo << "Transfinite Surface {1};\n";
    archivoGeo << "Recombine Surface {1};\n";

    archivoGeo.close();
}

int main() {
    // Supongamos que tienes una máscara binaria que representa tu objeto
    std::vector<std::vector<int>> mascaraBinaria = {
        {1, 1, 0, 0},
        {0, 1, 1, 0},
        {0, 0, 1, 1},
        {0, 0, 1, 1}
    };

    // Generar el archivo .geo
    generarArchivoGeo("modelo_3d.geo", mascaraBinaria);

    std::cout << "Archivo .geo generado exitosamente." << std::endl;

    return 0;
}
