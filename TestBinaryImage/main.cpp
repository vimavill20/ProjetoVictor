#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "pzcmesh.h"
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "pzmanvector.h"
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

    std::cout << std::endl;
    cv::Mat img = cv::imread("/Users/victorvillegassalabarria/Downloads/Sample000.tif", cv::IMREAD_COLOR);
    
            if (img.empty()) {
                std::cout << "No se pudo abrir la imagen\n";
                return -1;
            }
            std::cout << "rows " << img.rows << " cols " << img.cols << std::endl;
            for (int i = 0; i < 1; ++i) {
                for (int j = 0; j < 50; ++j) {
                    // Convertir cada píxel a binario y imprimirlo
                    std::cout << (int)img.at<uchar>(i, j) << ' ';
                }
                std::cout << '\n';
            }

    return 0;
}
