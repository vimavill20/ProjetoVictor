#include <iostream>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "pzmanvector.h"
#include "TPZGeoMeshTools.h"
#include "TPZRefPattern.h"
#include "TPZGenGrid2D.h"
#include "TPZVTKGeoMesh.h"
//using std::cout;
//using std::endl;
//using std::cin;
using namespace cv;
using namespace std;
int main (){
    std::cout<<"OpenCV VERSION "<< CV_VERSION <<std::endl;
    //Hola
    Mat src =cv::imread("/Users/victorvillegassalabarria/Documents/Github/ProjetoVictor2/SpiderTest2.jpeg",cv::IMREAD_GRAYSCALE);// cv::imread("/Users/victorvillegassalabarria/Documents/Github/ProjetoVictor/blackAndWhiteTEST.png",cv::IMREAD_GRAYSCALE);
    if (src.empty()) {
            std::cout << "Error: la imagen no se pudo cargar." << std::endl;
            return -1;
        }
//    cvtColor(src, src, COLOR_BGR2GRAY);
    Mat dst = Mat::zeros(src.rows, src.cols, CV_8UC3);
        src = src > 1;
        namedWindow( "Initial Image", 1 );
        imshow( "Initial Image", src );
        vector<vector<Point> > contours;
        vector<Vec4i> hierarchy;
        findContours( src, contours, hierarchy,
            RETR_CCOMP, CHAIN_APPROX_SIMPLE );
        // iterate through all the top-level contours,
        // draw each connected component with its own random color
        int idx = 0;
        for( ; idx >= 0; idx = hierarchy[idx][0] )
        {
            Scalar color(255,255,255 );//rand()&255, rand()&255, rand()&255 );
            drawContours( dst, contours, idx, color, FILLED, 8, hierarchy );
            Scalar colorverde(0,0,255 );
            drawContours( dst, contours, idx, colorverde, 1.5, 8, hierarchy );
        }
        namedWindow( "Recognized Image", 1 );
        imshow( "Recognized Image", dst );
        waitKey(0);

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
   

    return 0;
}
