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
#include <gmsh.h>
#include <fstream>
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "TPZDarcyFlow.h"
#include "TPZLinearAnalysis.h"


//using std::cout;
//using std::endl;
//using std::cin;
using namespace cv;
using namespace std;
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
void procesarImagen(std::string rutaImagen) {
    std::ofstream file("Coordenadas.txt");

    file << "OpenCV VERSION " << CV_VERSION << std::endl;

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
//function to make a darcy test
//void DarcyTest(){
//    //Creating a geometric mesh
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
//    TPZManVector<int,3> nel(2,1);
//    nel[0] = 2;
//    nel[1] = 2;
//    TPZGenGrid gengrid(nel,x0,x1);
//    gengrid.SetElementType(MMeshType::ETriangular);
//    gengrid.Read(gmesh);
//    gmesh->SetDimension(2);
//    gmesh->BuildConnectivity();
//    //Creating a computational mesh
//    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
//    //Creating a material
//    TPZDarcyFlow *material = new TPZDarcyFlow(1);
//    cmesh->InsertMaterialObject(material);
//    //Creating a boundary condition
//    TPZFMatrix<STATE> val1(1,1,0.),val2(1,1,0.);
//    TPZBndCond *bc = material->CreateBC(material, -1, 0, val1, val2);
//    cmesh->InsertMaterialObject(bc);
//    //Creating a computational element
//    cmesh->AutoBuild();
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
//    //Solving the system
//    TPZAnalysis an(cmesh);
//    TPZSkylineStructMatrix skyl(cmesh);
//    an.SetStructuralMatrix(skyl);
//    TPZStepSolver<STATE> step;
//    step.SetDirect(ELDLt);
//    an.SetSolver(step);
//    an.Run();
//    //Post processing
//    TPZStack<std::string> scalnames,vecnames;
//    scalnames.Push("Pressure");
//    vecnames.Push("Flux");
//    std::string plotfile("DarcyTest.vtk");
//    an.DefineGraphMesh(2,scalnames,vecnames,plotfile);
//    an.PostProcess(0);
//}
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
        }
        
    }
}
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
    
    //std::ofstream out("MatrixTest1.txt");
    //matrix.Print(out);
    //std::string outfilename2 = "outputMatrixPZ2.vtk"; // Nombre del archivo de salida
    // Devolver la matriz binaria
    
 //   std::cout<<matrix<<std::endl;
    return matrix;
}

//
//EJEMPLO SENCILLO PARA COMPROBAR VISUALVTK3D
//
////Primera Matriz
//TPZFMatrix<double> matrix(3, 3,0.); // Crear una matriz 3x3
//std::cout<<"Primera Matriz \n"<<std::endl;
//matrix(0,2)=1.0;
//std::cout<<matrix<<std::endl;
////Segunda Matriz
//TPZFMatrix<double> matrix2(3, 3,0.); // Crear una matriz 3x3
//matrix2(1,1)=1.0;
//std::cout<<"Segunda Matriz \n"<<std::endl;
//std::cout<<matrix2<<std::endl;
////Tercera Matriz
//TPZFMatrix<double> matrix3(3, 3,0.); // Crear una matriz 3x3
//matrix3(2,2)=1.0;
//std::cout<<"Tercera Matriz \n"<<std::endl;
//std::cout<<matrix3<<std::endl;
//TPZVec<TPZFMatrix<double>> Vecmatrix(3);
//Vecmatrix[0]=matrix;
//Vecmatrix[1]=matrix2;
//Vecmatrix[2]=matrix3;
//std::string outfilename = "Matrix3D.vtk"; // Nombre del archivo de salida
//VisualMatrix3DVTK(Vecmatrix, outfilename); // Llamar a la función
struct AllSimulationData
{
    // geometric element index of the interface element
//    std::vector<int> sim_order;
    // computational element index of the interface element
    std::vector<std::string> TIFnames;

};
TPZFMatrix<double> Case10TIF(AllSimulationData alldata, int idata){
    std::string filename=alldata.TIFnames[idata];
    filename=filename;
    cv::Mat img=cv::imread(filename, cv::IMREAD_GRAYSCALE);
    TPZFMatrix<double> matrix=create_binary_matrix(img);
    return matrix;
}
int main (){
    AllSimulationData test;
    
    std::string common_name="/Users/victorvillegassalabarria/Downloads/Dados_TIF/";
    
//    std::string name_1 ="Sample000.tif";
//    std::string name_2 ="Sample001.tif";
//    std::string name_3 ="Sample002.tif";
//    std::string name_4 ="Sample003.tif";
//    std::string name_5 ="Sample004.tif";
//    std::string name_6 ="Sample005.tif";
//    std::string name_7 ="Sample006.tif";
//    std::string name_8 ="Sample007.tif";
//    std::string name_9 ="Sample008.tif";
//    std::string name_10="Sample009.tif";
//    test.TIFnames.push_back(common_name + name_1);
//    test.TIFnames.push_back(common_name + name_2);
//    test.TIFnames.push_back(common_name + name_3);
//    test.TIFnames.push_back(common_name + name_4);
//    test.TIFnames.push_back(common_name + name_5);
//    test.TIFnames.push_back(common_name + name_6);
//    test.TIFnames.push_back(common_name + name_7);
//    test.TIFnames.push_back(common_name + name_8);
//    test.TIFnames.push_back(common_name + name_9);
//    test.TIFnames.push_back(common_name + name_10);
    std::vector<std::string> TIFnames;
    int nsim=910;
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
//        std::cout<<"iteracion: "<<i<<std::endl;
//        std::cout<<Vecmat10[i].Rows()<<std::endl;
    }
    std::string outfilename = "Matrix3D10.vtk"; // Nombre del archivo de salida
    VisualMatrix3DVTK(Vecmat10, outfilename);
    return 0;
}
    
//    std::string nombreSalida = "malla.msh";
//    std::string file_name="TriangleVictor.msh";
//        int nx=10;
//        int ny=10;
//        double L =1;
//        TPZVec<int> nels(3,0);
//        nels[0]=nx;         //Elements over x
//        nels[1]=ny;         //Elements over y
//
//        TPZVec<REAL> x0(3,0.0);
//        TPZVec<REAL> x1(3,0.0);
//        x1[0]=L;
//        x1[1]=L;
//
//        TPZGeoMesh *gmesh = new TPZGeoMesh;
//        TPZGenGrid2D gen(nels,x0,x1);
//
//        gen.SetElementType(MMeshType::ETriangular);
//        gen.Read(gmesh);
//        gmesh->SetDimension(2);
//        gmesh->BuildConnectivity();
//
//        int Nnels = gmesh->NElements();
//        for (int iel=0; iel<Nnels; iel++) {
//            TPZGeoEl * gel = gmesh->Element(iel);
//            int nsides = gel->NSides();
//            int nodes = gel->NNodes();
//            int firstside = nodes;
//            for(int iside = firstside; iside<nsides-1; iside++){
//                TPZGeoElSide gelside(gel,iside);
//                int nneig = gelside.NNeighbours();
//                if(nneig==0){
//
//                    int index1 = gelside.SideNodeIndex(0);
//                    int index2 = gelside.SideNodeIndex(1);
//                    auto node1 = gmesh->NodeVec()[index1];
//                    auto node2 = gmesh->NodeVec()[index2];
//                    double x_1 = node1.Coord(0);
//                    double y_1 = node1.Coord(1);
//
//                    double x_2 = node2.Coord(0);
//                    double y_2 = node2.Coord(1);
//                    if ((y_1 == y_2 && y_2 == 0.0) || (y_1 == y_2 && y_2 == L)) {
//                        gel->CreateBCGeoEl(iside, 2);
//                    }
//                    if ((x_1 == x_2 && x_2 == 0.0)) {
//                        gel->CreateBCGeoEl(iside, 3);
//                    }
//                    if ( (x_1 == x_2 && x_2 == L)) {
//                        gel->CreateBCGeoEl(iside, 4);
//                    }
//                }
//            }
//
//
//        }
//
//
//    gmesh->Print(std::cout);
//    std::string nombreSalida = "testVictor1.vtk";
//    std::ofstream file(nombreSalida);
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
//    int dim=2;
//    int order=1;
//    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
//    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
//    TPZDarcyFlow *material= new TPZDarcyFlow(1,2);
//    TPZMaterial *mat(material);
////    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
//
////    TPZDarcyFlow(mat,dim) ;
//    cmesh->InsertMaterialObject(mat);
//    cmesh->SetDimModel(dim);
//    cmesh->SetDefaultOrder(order);
//    cmesh->AutoBuild();
//    TPZLinearAnalysis AN;
//    TPZEigenSolver<STATE> step;
//
//    AN.SetSolver(step);
//    cmesh->Print();
//    std::ofstream file1("cmesh_h.txt");
//    cmesh->Print(file1);
//    TPZVTKGeoMesh
    //    TPZLinearAnalysis
//    TPZLinearAnalysis cmsh;
//    TPZPrint
//    return 0;


//generarMalla("Malha_basica", 1.0, 0.1);
//    procesarImagen("/Users/victorvillegassalabarria/Documents/Github/ProjetoVictor2/Aranha.png");
//---------------------------------
// MALHA HOMOGENEA USANDO NEOPZ
//    int nx=10;
//    int ny=10;
//    double L =1;
//    TPZVec<int> nels(3,0);
//    nels[0]=nx;         //Elements over x
//    nels[1]=ny;         //Elements over y
//
//    TPZVec<REAL> x0(3,0.0);
//    TPZVec<REAL> x1(3,0.0);
//    x1[0]=L;
//    x1[1]=L;
//
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    TPZGenGrid2D gen(nels,x0,x1);
//
//    gen.SetElementType(MMeshType::ETriangular);
//    gen.Read(gmesh);
//    gmesh->SetDimension(2);
//    gmesh->BuildConnectivity();
//
//    int Nnels = gmesh->NElements();
//    for (int iel=0; iel<Nnels; iel++) {
//        TPZGeoEl * gel = gmesh->Element(iel);
//        int nsides = gel->NSides();
//        int nodes = gel->NNodes();
//        int firstside = nodes;
//        for(int iside = firstside; iside<nsides-1; iside++){
//            TPZGeoElSide gelside(gel,iside);
//            int nneig = gelside.NNeighbours();
//            if(nneig==0){
//
//                int index1 = gelside.SideNodeIndex(0);
//                int index2 = gelside.SideNodeIndex(1);
//                auto node1 = gmesh->NodeVec()[index1];
//                auto node2 = gmesh->NodeVec()[index2];
//                double x_1 = node1.Coord(0);
//                double y_1 = node1.Coord(1);
//
//                double x_2 = node2.Coord(0);
//                double y_2 = node2.Coord(1);
//                if ((y_1 == y_2 && y_2 == 0.0) || (y_1 == y_2 && y_2 == L)) {
//                    gel->CreateBCGeoEl(iside, 2);
//                }
//                if ((x_1 == x_2 && x_2 == 0.0)) {
//                    gel->CreateBCGeoEl(iside, 3);
//                }
//                if ( (x_1 == x_2 && x_2 == L)) {
//                    gel->CreateBCGeoEl(iside, 4);
//                }
//            }
//        }
//
//
//    }
//
//
//    gmesh->Print(std::cout);
//
//
//    std::ofstream file("victorTest2.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
