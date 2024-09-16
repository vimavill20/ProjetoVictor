//#include <fstream>
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
#include <fstream>
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "DarcyFlow/TPZDarcyFlow.h"
#include "TPZLinearAnalysis.h"
#include "TPZExtendGridDimension.h"
#include "pzvisualmatrix.h"
#include "EditImage.h"
//#include "TestVTK3D.h"
using namespace cv;
using namespace std;
//FUNCION PARA CREAR UNA MATRIZ DE ARCHIVO BINARIO .RAW
TPZVec<TPZFMatrix<double>> create_raw_Vecmatrix(std::string rutaArchivo, int filas, int columnas, int layers){
    std::ifstream file(rutaArchivo, std::ios::binary | std::ios::ate);
    TPZFMatrix<double> matrix(filas,columnas);
    TPZVec<TPZFMatrix<double>> Vectormatrix(layers);
//    int layers=3;
    if (file.is_open()) {
        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);

        std::vector<char> buffer(size);
        file.read(buffer.data(), size);

        int count = 0;
        for(int lay=0;lay < layers;lay++){
        for(int i = 0; i < filas; i++) {
            for(int j = 0; j <columnas ; ++j) {
                if(count >= size) break;
                matrix(i,j) = buffer[count];
                count++;
            }
        }Vectormatrix[lay]=matrix;
        }
    } else {
        std::cout << "Unable to open file";
        DebugStop();
    }

    return Vectormatrix;
}


//CODIGO PARA CREAR VTK A PARTIR DE ARCHIVO .RAW
void VisualMatrix3DVTK(TPZVec<TPZFMatrix<double>> & matrixVector, const std::string &outfilename)
{
    const int nelz = matrixVector.size();
    TPZFMatrix<double> & matrix1 = matrixVector[0];
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

//CODIGO PARA CREAR VTK A PARTIR DE ARCHIVO .RAW
void VisualMatrix3DVTK(Image3D & image, const std::string &outfilename)
{
    const int nelz = image.Depth();
    const int nelx = image.Height();
    const int nely = image.Width();
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
    out << "SCALARS " << image.VarName() << " float 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int k = 0; k < nelz; k++)
    {
        for (j=0; j<nely; j++) {
            for (i=0;i<nelx;i++) {
                out << image.getPixel(k,j,i) << std::endl;
            }
        }
    }
    out << std::endl;
}

//CODIGO PARA CREAR VTK A PARTIR DE ARCHIVO .RAW
void AddDataVTK(Image3D & image, const std::string &outfilename)
{
    const int nelz = image.Depth();
    const int nelx = image.Height();
    const int nely = image.Width();
    const int neltotal = nelx * nely * nelz;
    int i,j;
    //ofstream out(outfilename.c_str(),std::ios::app);
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
    out << "SCALARS " << image.VarName() << " float 1\n";
    out << "LOOKUP_TABLE default\n";
    for(int k = 0; k < nelz; k++) {
        for (int j=0; j<nely; j++) {
            for (int i=0;i<nelx;i++) {
                out << image.getPixel(k,j,i) << std::endl;
            }
        }
    }
    out << std::endl;
}

int mainfake();
int maintestfraturaaislada();
int maintestObjectsinPlane();
int mainnewImages();
int mainpython();
int main() {
    return mainpython();
    //return mainfake();
    //return maintestfraturaaislada();
    //return maintestObjectsinPlane();
//  //  return mainnewImages();
}

int mainnewImages(){
//    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/ImagesICPSC2015/Bentheimer_1000c_3p0035um.raw";
    
//    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/ImagesICPSC2015/Estaillades_1000c_3p31136um.raw";
//        std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/ImagesICPSC2015/Doddington_1000c_2p6929um.raw";
    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/ImagesICPSC2015/Ketton_1000c_3p00006um.raw";

        int layers=100;
        int rows=100;
        int cols=100;
        TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
//        std::cout << "Vector of matrixs is created " << outfilename1 << "\n";

    
//        Image3D input("pixeldata",VecOfMat);
//        Image3D output("objects",input.Depth(), input.Width(), input.Height());
//        Image3D ordered("ordered",input.Depth(), input.Width(), input.Height());
//        int count = input.identifyObjects(output);
////        std::cout << "Identified " << count << " objects." <<"\n"<< std::endl;
//        output.orderObjectsBySize(ordered, count);
////        std::cout << count << "\n";
//        std::string outfilename ="RAWTestBentheimer.vtk";
//        std::string outfilename1 ="RAWTestBentheimeroutput.vtk";
//        std::string outfilename2 ="RAWTestBentheimerordered.vtk";
//
//        VisualMatrix3DVTK(input, outfilename);
//        std::cout << "adding object data to " << outfilename1 << "\n";
//        AddDataVTK(output, outfilename1);
        Image3D input("pixeldata",VecOfMat);
        Image3D output("objects",input.Depth(), input.Width(), input.Height());
        Image3D ordered("ordered",input.Depth(), input.Width(), input.Height());
        int count = input.identifyObjects(output);
        output.orderObjectsBySize(ordered, count);
        std::cout << count << "\n";
        std::string outfilename ="RAWTest.vtk";
        std::string outfilename1 ="RAWTestoutput.vtk";
        std::string outfilename2 ="RAWTestordered.vtk";
        VisualMatrix3DVTK(input, outfilename);
        std::cout << "adding object data to " << outfilename1 << "\n";
        AddDataVTK(output, outfilename1);
        std::cout << "adding ordered data to " << outfilename2 << "\n";
        AddDataVTK(ordered, outfilename2);
    }






//CODIGO PRINCIPAL
int mainfake (){
//    std::string ImagenRaw="Sample_Labels_3D_RAW.raw";
    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/Sample_Labels_3D_RAW1.raw";
    
    int layers=30;
    int rows=676;//676;
    int cols=616;//616;
    TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
    std::cout << "Vecofmat created " << std::endl;

//    TPZFMatrix<double> matrix(rows,cols);
//    TPZVec<TPZFMatrix<double>> VecOfMat(layers,matrix);
//    for (int layer = 0; layer < layers; ++layer) {
//            for (int row = 0; row < rows; ++row) {
//                for (int col = 0; col < cols; ++col) {
//                    if(col==2&row==2){
//                        VecOfMat[1](row,col) = 1;}
//
//                    else{
//                        VecOfMat[layer](row,col) = 0;}
//                }
//            }
//        }
//    VecOfMat[0](4,4) = 1;
//    VecOfMat[0](4,3) = 1;

//    std::cout<<VecOfMat<<std::endl;
    Image3D input("pixeldata",VecOfMat);
    std::cout << "Image 3d created " << std::endl;

    Image3D output("objects",input.Depth(), input.Width(), input.Height());
    std::cout << "ImageOUT 3d created " << std::endl;

    Image3D ordered("ordered",input.Depth(), input.Width(), input.Height());
    std::cout << "ImageOrd 3d created " << std::endl;

    int count = input.identifyObjects(output);
    std::cout << "Identified " << count << " objects." <<"\n"<< std::endl;
    output.orderObjectsBySize(ordered, count);
    std::cout << count << "\n";
    std::string outfilename ="RAWTest.vtk";
    std::string outfilename1 ="RAWTestoutput.vtk";
    std::string outfilename2 ="RAWTestordered.vtk";


    
    VisualMatrix3DVTK(input, outfilename);
    std::cout << "adding object data to " << outfilename1 << "\n";
    AddDataVTK(output, outfilename1);
    std::cout << "adding ordered data to " << outfilename2 << "\n";
    AddDataVTK(ordered, outfilename2);
    
    int etiquetaObjeto = 100; // Puedes cambiar esto con la etiqueta del objeto que deseas analizar
    int numPixeles = output.getPixelsInObject(etiquetaObjeto);
    std::cout << "El objeto con etiqueta " << etiquetaObjeto << " tiene " << numPixeles << " píxeles." << std::endl;
    int etiquetaObjetoOrd = 100; // Puedes cambiar esto con la etiqueta del objeto que deseas analizar
    int numPixelesOrd = ordered.getPixelsInObject(etiquetaObjetoOrd);
    std::cout << "El objeto con etiqueta( objetos ordenados) " << etiquetaObjeto << " tiene " << numPixelesOrd << " píxeles." << std::endl;
    int numPixeles1 = output.getPixelsInObject(1);
    std::cout << "El objeto con etiqueta " << "1" << " tiene " << numPixeles1 << " píxeles." << std::endl;
    int numberBuracos=0;
//
    //

//    std::vector<int> ObjectsPixels(2,0);
//    for(int i=1;i<count+1;i++){
//        numPixeles=output.getPixelsInObject(i);
//        numberBuracos+=numPixeles;
//        numPixeles=0;
//    }
//    std::cout << "El numero  total de pixeles en objetos reconocidos es " << numberBuracos << std::endl;
//
    
        // Imprimir el contenido del vector
//        std::cout << "Contenido del vector:" << std::endl;
//    TPZVec<double> VecPixelsbyObjects(count);
//    VecPixelsbyObjects= output.obtenerObjetosYPixeles(output, count);
//    //std::cout << VecPixelsbyObjects <<std::endl;
//    for (int i = 0; i < VecPixelsbyObjects.size(); i++) {
//        std::cout << "Objeto " << i << ", Píxeles " << VecPixelsbyObjects[i] << std::endl;
//        }
   
        // Call the countFacesByObject method

   


//
    //output.countFacesByObject(facesCount, input);

    return 0;
}

int maintestfraturaaislada(){
    int layers=3;
    int rows=5;//676;
    int cols=5;
    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/Sample_Labels_3D_RAW1.raw";
    
    //TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
    
    TPZFMatrix<double> matrix(rows,cols);
    
    TPZVec<TPZFMatrix<double>> VecOfMat(layers,matrix);
    for (int layer = 0; layer < layers; ++layer) {
            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    if(col==2&row==2){
                        VecOfMat[1](row,col) = 1;}
                    
                    else{
                        VecOfMat[layer](row,col) = 0;}
                }
            }
        }
    Image3D input("pixeldata",VecOfMat);
    Image3D output("objects",input.Depth(), input.Width(), input.Height());
    Image3D ordered("ordered",input.Depth(), input.Width(), input.Height());
    Image3D inputbranco("objects",input.Depth(), input.Width(), input.Height());

    int count = input.identifyObjects(output);
    output.orderObjectsBySize(ordered, count);

    std::string outfilename ="RAWTest.vtk";
    std::string outfilename1 ="RAWTestoutput.vtk";
    std::string outfilename2 ="RAWTestordered.vtk";
    VisualMatrix3DVTK(input, outfilename);
    std::cout << "adding object data to " << outfilename1 << "\n";
    AddDataVTK(output, outfilename1);
    std::cout << "adding ordered data to " << outfilename2 << "\n";
    AddDataVTK(ordered, outfilename2);
    //Cubo encerado
    //TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
    //setPixel layer,col,row
    //inputbranco.setPixel(1, 2, 2, 1);
    std::string outfilenamebr="ArquivoBranco.vtk";
    std::cout << "adding ordered data to " << outfilenamebr << "\n";
    std::string outfilenamebr1="ArquivoBrancoComObjeto1.vtk";
    AddDataVTK(inputbranco, outfilenamebr);
    ordered.highlightObject(ordered,inputbranco,1);
    AddDataVTK(inputbranco, outfilenamebr1);
    return 0;
    
}
int maintestObjectsinPlane(){
    int layers=20;
    int rows=20;//676;
    int cols=20;
    //std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/Sample_Labels_3D_RAW1.raw";
    
    //TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
    
    TPZFMatrix<double> matrix(rows,cols);

    TPZVec<TPZFMatrix<double>> VecOfMat(layers,matrix);
    for (int layer = 0; layer < layers; ++layer) {
            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    if(col==1&row==3){
                        VecOfMat[layer](row,col) = 1;
                        VecOfMat[layer](row-2,col+2) = 1;
                        VecOfMat[layer](row-2,col) = 1;

                    }


                    else{
                        VecOfMat[layer](row,col) = 0;}
                }
            }
        }
    for (int row = 0; row < rows; ++row) {
            
        VecOfMat[10](row,10) = 1;
        VecOfMat[10](row,12) = 1;
    }
        
    
        
        
    VecOfMat[1](4,4) = 1;
    VecOfMat[4](1,2) = 1;
    VecOfMat[15](15,15) = 1;
    VecOfMat[10](12,12) = 1;

    
    
    std::cout<<VecOfMat<<std::endl;
    Image3D input("pixeldata",VecOfMat);
    Image3D output("objects",input.Depth(), input.Width(), input.Height());
    Image3D ordered("ordered",input.Depth(), input.Width(), input.Height());
    Image3D inputplane("objects",input.Depth(), input.Width(), input.Height());

    int count = input.identifyObjects(output);
    output.orderObjectsBySize(ordered, count);

    std::string outfilename ="RAWTest.vtk";
    std::string outfilename1 ="RAWTestoutput.vtk";
    std::string outfilename2 ="RAWTestordered.vtk";
    VisualMatrix3DVTK(input, outfilename);
    std::cout << "adding object data to " << outfilename1 << "\n";
    AddDataVTK(output, outfilename1);
    std::cout << "adding ordered data to " << outfilename2 << "\n";
    AddDataVTK(ordered, outfilename2);
    //Cubo encerado
    //TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
    //setPixel layer,col,row
    //inputbranco.setPixel(1, 2, 2, 1);
    std::string outfilenamebr="ArquivoBranco.vtk";
    std::cout << "adding ordered data to " << outfilenamebr << "\n";
    std::string outfilenamebr1="ArquivoBrancoComObjeto1.vtk";
    AddDataVTK(inputplane, outfilenamebr);
    
    output.Objects3DinPlane(output,inputplane,1);
    //output.highlightObject(output,inputplane,1);
    //output.highlightObject(output,inputplane,2);
    //output.highlightObject(output,inputplane,3);
    //output.highlightObject(output,inputplane,4);


    
    AddDataVTK(inputplane, outfilenamebr1);
}
int mainpython(){
//    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/geometrias.raw";
    //std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/Sample_Labels_3D_RAW1.raw";
    
//    int layers=20;//30
//    int rows=20;//676;
//    int cols=20;//616;
//
//    //TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
//    TPZFMatrix<double> matrix(rows,cols);
//
//    TPZVec<TPZFMatrix<double>> VecOfMat(layers,matrix);
//    for (int layer = 0; layer < layers; ++layer) {
//            for (int row = 0; row < rows; ++row) {
//                for (int col = 0; col < cols; ++col) {
//                    if(col==1&row==3){
//                        VecOfMat[layer](row,col) = 1;
//                        VecOfMat[layer](row-2,col+2) = 1;
//                        VecOfMat[layer](row-2,col) = 1;
//
//                    }
//
//
//                    else{
//                        VecOfMat[layer](row,col) = 0;}
//                }
//            }
//        }
//    for (int row = 0; row < rows; ++row) {
//
//        VecOfMat[10](row,10) = 1;
//        VecOfMat[10](row,12) = 1;
//    }
//
//
//
//
//    VecOfMat[1](4,4) = 1;
//    VecOfMat[4](1,2) = 1;
//    VecOfMat[15](15,15) = 1;
//    VecOfMat[10](12,12) = 1;
    int layers=150;
    int rows=676;//676;
    int cols=616;
    std::string ImagenRaw="/Users/victorvillegassalabarria/Downloads/Sample_Labels_3D_RAW1.raw";
    
    TPZVec<TPZFMatrix<double>> VecOfMat=create_raw_Vecmatrix(ImagenRaw, rows, cols, layers);
    
    TPZFMatrix<double> matrix(rows,cols);
    
//    TPZVec<TPZFMatrix<double>> VecOfMat(layers,matrix);
//    for (int layer = 0; layer < layers; ++layer) {
//            for (int row = 0; row < rows; ++row) {
//                for (int col = 0; col < cols; ++col) {
//                    if(col==2&row==2){
//                        VecOfMat[1](row,col) = 1;}
//
//                    else{
//                        VecOfMat[layer](row,col) = 0;}
//                }
//            }
//        }
    
    
    //std::cout<<VecOfMat<<std::endl;
    std::cout << "Vecofmat created " << std::endl;
    Image3D input("pixeldata",VecOfMat);
    Image3D output("objects",input.Depth(), input.Width(), input.Height());
    Image3D ordered("ordered",input.Depth(), input.Width(), input.Height());
    Image3D Segm("ordered",input.Depth(), input.Width(), input.Height());
    int count = input.identifyObjects(output);
    output.orderObjectsBySize(ordered, count);
    std::cout << count << "\n";
    std::string outfilename ="RAWTest.vtk";
    std::string outfilename1 ="RAWTestoutput.vtk";
    std::string outfilename2 ="RAWTestordered.vtk";
    std::string filen="Voxel_cords.txt";
//    ordered.SegmentVugFracture(Segm,30,filen);
    VisualMatrix3DVTK(input, outfilename);
    std::cout << "adding object data to " << outfilename1 << "\n";
    AddDataVTK(output, outfilename1);
    std::cout << "adding ordered data to " << outfilename2 << "\n";
    AddDataVTK(ordered, outfilename2);
    ordered.getTxtPixelsInObject(5,filen);

}
