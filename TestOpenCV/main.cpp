#include <iostream>
#include <filesystem>
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
#include <gmsh.h>
#include <fstream>
#include "TPZMaterial.h"
#include "TPZDarcyFlow.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZNullMaterialCS.h"
#include "TPZNullMaterial.h"
//#include "TPZAnalysis.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
//#include "TPZStepSolver.h"
//using std::cout;
//using std::endl;
//using std::cin;
//hola
using namespace cv;
using namespace std;
//Function to generate a mesh using gmsh library
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine);
void generarMalla(std::string nombreArchivo, double L, double lc) {
    gmsh::initialize();
    gmsh::model::add(nombreArchivo);
    
    // Se crean los puntos
    int p1 = gmsh::model::geo::addPoint(0, 0, 0, lc);
    int p2 = gmsh::model::geo::addPoint(L, 0, 0, lc);
    int p3 = gmsh::model::geo::addPoint(L, L, 0, lc);
    int p4 = gmsh::model::geo::addPoint(0, L, 0, lc);
    
    // Se crean las lineas
    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p4);
    int l4 = gmsh::model::geo::addLine(p4, p1);
    
    // Se crean las curvas
    int cl = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
    
    // Se crean las superficies
    int pl = gmsh::model::geo::addPlaneSurface({cl});
    gmsh::model::geo::synchronize();
    
    // Configurar para generar malla de cuadrados
    // gmsh::option::setNumber("Mesh.RecombineAll", 1);
    
    // Se genera la malla 2D
    gmsh::model::mesh::generate(2);
    gmsh::write(nombreArchivo + ".msh");
    gmsh::finalize();
}
//Function to create a txt with the coordinates of the contourns of binary images
void procesarImagen(std::string rutaImagen) {
    
    std::ofstream file("Coordenadas.txt");
    // Encuentra la última barra y el punto en la cadena
    size_t ultimaBarra = rutaImagen.find_last_of("/");
    size_t punto = rutaImagen.find_last_of(".");

    // Obtiene el nombre del archivo sin la extensión
    std::string nombre = rutaImagen.substr(ultimaBarra + 1, punto - ultimaBarra - 1);

    file << "OpenCV VERSION " << CV_VERSION <<std::endl;
    file << "Imagen " <<nombre <<std::endl;
    Mat src = cv::imread(rutaImagen, cv::IMREAD_GRAYSCALE);
    if (src.empty()) {
        file << "Error: la imagen no se pudo cargar." << std::endl;
        return;
    }

    Mat dst = Mat::zeros(src.rows, src.cols, CV_8UC3);
    src = src > 1;
    namedWindow("Initial Image", 1);
    imshow("Initial Image", src);

    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    findContours(src, contours, hierarchy, RETR_CCOMP, CHAIN_APPROX_SIMPLE);

    int idx = 0;
    for(; idx >= 0; idx = hierarchy[idx][0]) {
        Scalar color(255, 255, 255);
        drawContours(dst, contours, idx, color, FILLED, 8, hierarchy);
        Scalar colorverde(0, 0, 255);
        drawContours(dst, contours, idx, colorverde, 1.5, 8, hierarchy);
    }

    for(size_t i = 0; i < contours.size(); i++) {
        for(size_t j = 0; j < contours[i].size(); j++) {
            file << "Contorno " << i << ", Punto " << j << ": " << contours[i][j] << std::endl;
        }
    }

    file.close();
}
//Function to generate a mesh using neopz library

TPZGeoMesh* crearMallaHomogenea(int nx, int ny, double L, std::string nombreSalida) {
    TPZVec<int> nels(3,0);
    nels[0]=nx;         //Elements over x
    nels[1]=ny;         //Elements over y

    TPZVec<REAL> x0(3,0.0);
    TPZVec<REAL> x1(3,0.0);
    x1[0]=L;
    x1[1]=L;

    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGenGrid2D gen(nels,x0,x1);

    gen.SetElementType(MMeshType::ETriangular);
    gen.Read(gmesh);
    gmesh->SetDimension(2);
    gmesh->BuildConnectivity();

    int Nnels = gmesh->NElements();
    for (int iel=0; iel<Nnels; iel++) {
        TPZGeoEl * gel = gmesh->Element(iel);
        int nsides = gel->NSides();
        int nodes = gel->NNodes();
        int firstside = nodes;
        for(int iside = firstside; iside<nsides-1; iside++){
            TPZGeoElSide gelside(gel,iside);
            int nneig = gelside.NNeighbours();
            if(nneig==0){

                int index1 = gelside.SideNodeIndex(0);
                int index2 = gelside.SideNodeIndex(1);
                auto node1 = gmesh->NodeVec()[index1];
                auto node2 = gmesh->NodeVec()[index2];
                double x_1 = node1.Coord(0);
                double y_1 = node1.Coord(1);

                double x_2 = node2.Coord(0);
                double y_2 = node2.Coord(1);
                if ((y_1 == y_2 && y_2 == 0.0) || (y_1 == y_2 && y_2 == L)) {
                    gel->CreateBCGeoEl(iside, 2);
                }
                if ((x_1 == x_2 && x_2 == 0.0)) {
                    gel->CreateBCGeoEl(iside, 3);
                }
                if ( (x_1 == x_2 && x_2 == L)) {
                    gel->CreateBCGeoEl(iside, 4);
                }
            }
        }
    }
    std::ofstream file(nombreSalida);
    gmesh->Print(file);
    return gmesh;
//    std::ofstream file(nombreSalida);
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
}

#include <fstream>
#include <opencv2/opencv.hpp>
//Function to see 0 and 1 data from .raw files
cv::Mat readRawFile(const std::string& filename, int width, int height) {
    // Create a matrix to hold the image data
    cv::Mat img(height, width, CV_8UC1);

    // Open the file
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Failed to open the file." << std::endl;
        return cv::Mat();
    }

    // Read the raw data
    file.read(reinterpret_cast<char*>(img.data), width * height);

    return img;
}


int main (){
  
    TPZVec<int> nx(2, 10);
    //What does it mean the line 46?
    
    nx[0]=3;
    nx[1]=1;
    const TPZVec<REAL> x0(3, 0.);
    const TPZVec<REAL> x1(3, 0.);
    x1[0]=1.0;
    x1[1]=1.0;
    auto msh = TPZGenGrid2D(nx, x0, x1);
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    
    msh.Read(gmesh);
    msh.SetElementType(MMeshType::EQuadrilateral);
    
    msh.SetBC(gmesh, 7, 2);
    msh.SetBC(gmesh, 5, 3);
    msh.SetBC(gmesh, 4, 4);
    msh.SetBC(gmesh, 6, 4);
    std::ofstream file2("TestGeoMesh2D.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file2);
    
    
    
    
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagCoarse(4);
    dim_name_and_physical_tagCoarse[2]["k11"] = 1;
    dim_name_and_physical_tagCoarse[1]["inlet"] = 2;
    dim_name_and_physical_tagCoarse[1]["outlet"] = 3;
    dim_name_and_physical_tagCoarse[1]["noflux"] = 4;
    
    std::string filename="/Users/victorvillegassalabarria/Documents/Mastria/FEM2024/Malla2Dtestpointslinesap.msh";
    gmesh = generateGMeshWithPhysTagVec(filename, dim_name_and_physical_tagCoarse);
   
    std::ofstream file3("TestGeoMesh2Dnew.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file3);
    //Create CompMesh
    TPZCompMesh *cmesh =  new TPZCompMesh(gmesh);
    
    //Create Materials
    int matId=1;
    int dim = 2;
    TPZDarcyFlow *matDarcy = new TPZDarcyFlow(matId, dim);
    
    cmesh->InsertMaterialObject(matDarcy);
    int bc_id=2;
    int bc_typeN = 1;
    int bc_typeD = 0;
    TPZFMatrix<STATE> val1(1,1,0.0);
    TPZVec<STATE> val2(1,0.0);
    
    int bcinletId = 2;
    int bcOutletId = 3;
    int bcNoFlux = 4;
    TPZBndCond * face2 = matDarcy->CreateBC(matDarcy,bcNoFlux,bc_typeN,val1,val2);
    cmesh->InsertMaterialObject(face2);
    
    val2[0]=100; // Valor a ser impuesto como presión en la entrada
    TPZBndCond * face = matDarcy->CreateBC(matDarcy,bcinletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face);
    
    val2[0]=10; // Valor a ser impuesto como presión en la salida
    TPZBndCond * face1 = matDarcy->CreateBC(matDarcy,bcOutletId,bc_typeD,val1,val2);
    cmesh->InsertMaterialObject(face1);
    
   
    cmesh->AutoBuild();
    //Esto hace que el espacio de aproxiación sea H1
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    
    //Inicializa el tamaño del vector solución
    cmesh->ExpandSolution();
    
  
    //CreateAnalisys
    TPZLinearAnalysis *Analisys = new TPZLinearAnalysis(cmesh);
    bool mustOptimizeBandwidth = false;
   
    //Carga la solución a la malla computacional
    Analisys->LoadSolution();
    
   // Selecciona el método numérico para resolver el problema algebraico
    TPZStepSolver<STATE> step;
    
//    TPZSSpStructMatrix<STATE> matrix(cmesh);
      step.SetDirect(ELDLt);
  
//    Analisys->SetStructuralMatrix(matrix);
    

    Analisys->SetSolver(step);
    
    //Ensamblaje de la matriz de rigidez y vector de carga
    Analisys->Assemble();

    //Resolución del sistema algebraico
    Analisys->Solve();

    
    //Definición de variables escalares y vectoriales a posprocesar
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    
    //Configuración del posprocesamiento
    int ref =0; // Permite refinar la malla con la solucion obtenida
    std::string file_reservoir("SolVictor.vtk");
    Analisys->DefineGraphMesh(dim,scalnames,vecnames,file_reservoir);
    //Posprocesamiento
    Analisys->PostProcess(ref, dim);
 
    return 0;
}
int mainfake(){
cv::Mat img = cv::imread("/Users/victorvillegassalabarria/Downloads/Sample908.tif", cv::IMREAD_COLOR);

   // Comprueba si la imagen se ha cargado correctamente
   if(img.empty()) {
       std::cout << "Error al abrir la imagen" << std::endl;
       return -1;
   }


   // Binariza la imagen
   cv::Mat img_binaria;
//       cv::adaptiveThreshold(img, img_binaria, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 11, 2);
   cv::cvtColor(img, img_binaria, cv::COLOR_BGR2GRAY);
//       cv::threshold(img_binaria, img_binaria, 128, 255, cv::THRESH_BINARY);
  cv::adaptiveThreshold(img_binaria, img_binaria, 255, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY, 3, 3);



   // Guarda la imagen binaria
   cv::imwrite("/Users/victorvillegassalabarria/Downloads/908imagen_binaria.tif", img_binaria);
auto rutaImagen="/Users/victorvillegassalabarria/Downloads/908imagen_binaria.tif";
//    auto rutaImagen="/Users/victorvillegassalabarria/Documents/Github/ProjetoVictor2_build/BinaryImage.png";
procesarImagen(rutaImagen);
return 0;
}
TPZGeoMesh* generateGMeshWithPhysTagVec(std::string& filename, TPZManVector<std::map<std::string,int>,4>& dim_name_and_physical_tagFine){
            
    // Creating gmsh reader
    TPZGmshReader  GeometryFine;
    TPZGeoMesh *gmeshFine;
    REAL l = 1.0;
    GeometryFine.SetCharacteristiclength(l);
    
    // Reading mesh
    GeometryFine.SetDimNamePhysical(dim_name_and_physical_tagFine);
    gmeshFine = GeometryFine.GeometricGmshMesh(filename,nullptr,false);
    return gmeshFine;
}
