#include "EditImage.h"

void Image3D::InsertSurroundingPixels(const Image3D &input, Image3D &output, int x, int y, int z, int value) {
    using Pixel = std::tuple<int, int, int>;
    if (input.getPixel(x, y, z) == 0) DebugStop();
    output.setPixel(x, y, z, value);
    std::set<Pixel> toVisit, visited;
    toVisit.insert({x, y, z});
    while (!toVisit.empty()) {
        auto [x, y, z] = *toVisit.begin();
        toVisit.erase(toVisit.begin());
        if (visited.count({x, y, z})) {
            continue;
        }
        visited.insert({x, y, z});
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                for (int k = -1; k <= 1; k++) {
                    if (i == 0 && j == 0 && k == 0) {
                        continue;
                    }
                    if (x + i >= 0 && x + i < input.depth && y + j >= 0 && y + j < input.width && z + k >= 0 && z + k < input.height) {
                        if(visited.count({x + i, y + j, z + k})) {
                            continue;
                        }
                        if (input.getPixel(x + i, y + j, z + k) > 0) {
                            output.setPixel(x+i, y+j, z+k, value);
                            toVisit.insert({x + i, y + j, z + k});
                        }
                        else {
                            visited.insert({x + i, y + j, z + k});
                        }
                    }
                }
            }
        }
    }
}
//int Image3D::identifyObjects(Image3D& output)const {
//    // Assume this method identifies objects in the 3D image and writes the result to the output image.
//    int count = 1;
//    // Initialize the output image.
//    for (auto &mat : output.data) {
//        mat.Zero();
//    }
//    for (int i = 0; i < depth; i++) {
//        for (int j = 0; j < width; j++) {
//            for (int k = 0; k < height; k++) {
//                if (getPixel(i,j,k) > 0 && output.getPixel(i, j, k) == 0) {
//                    InsertSurroundingPixels(*this, output, i, j, k, count);
//                    count++;
//                }
//            }
//        }
//    }
//    return count;
//}
int Image3D::identifyObjects(Image3D& output) const {
    int count = 1;
    for (auto &mat : output.data) {
        mat.Zero();
    }

    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < height; k++) {
                if (getPixel(i, j, k) > 0 && output.getPixel(i, j, k) == 0) {
//                    std::cout << "Labeling object " << count << " at position (" << i << ", " << j << ", " << k << ")" << std::endl;
                    InsertSurroundingPixels(*this, output, i, j, k, count);
                    count++;
                }
            }
        }
    }

    return count;
}
/// order the objects by size
void Image3D::orderObjectsBySize(Image3D& output, int numcolors)
{
    std::vector<int> count(numcolors,0);
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < height; k++) {
                if (getPixel(i,j,k) > 0) {
                    count[getPixel(i,j,k)]++;
                }
            }
        }
    }
    std::vector<int> index(numcolors);
    for (int i = 0; i < numcolors; i++) {
        index[i] = i;
    }
    std::sort(index.begin(),index.end(),[&count](int i, int j){return count[i] > count[j];});
    std::vector<int> newindex(numcolors);
    for (int i = 0; i < numcolors; i++) {
        newindex[index[i]] = i;
    }
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < height; k++) {
                if (getPixel(i,j,k) > 0) {
                    output.setPixel(i,j,k,newindex[getPixel(i,j,k)]);
                }
            }
        }
    }
    std::map<int,int> objectsinrange;
    for (int i = 0; i < numcolors; i++) {
        int numpixels = count[i];
        int j = 1;
        int val = numpixels>>2*j;
        while(val > 0) {
            j++;
            val = numpixels>>2*j;
        }
        objectsinrange[j]++;
    }
    
    std::cout << "objectsinrange: " << std::endl;
    for (const auto& pair : objectsinrange) {
        //std::cout << "number of objects smaller than " << (pair.first) << ": " << pair.second << std::endl;
        std::cout << "number of objects smaller than " << (2<<2*pair.first) << ": " << pair.second << std::endl;

    }
}
int Image3D::getPixelsInObject(int label) const {
    int pixelsCount = 0;
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < height; k++) {
                if (getPixel(i, j, k) == label) {
                    pixelsCount++;
                }
            }
        }
    }
    return pixelsCount;
}
//TPZVector

TPZVec<double> Image3D::obtenerObjetosYPixeles(const Image3D& output, int colors) {
    TPZVec<double> objetosYPixeles(colors);

    for (int i = 0; i < colors; i++) {
        int numPixeles = output.getPixelsInObject(i);
        objetosYPixeles[i]=numPixeles;
    }

    return objetosYPixeles;
}
void Image3D::countFacesByObject(std::map<int, int>& facesCount, const Image3D& Image) const {
    Image3D output("objects", Image.Depth(), Image.Width(), Image.Height());
    int count = Image.identifyObjects(output);

    // Initialize face count for each object
    for (int i = 1; i <= count; i++) {
        facesCount[i] = 0;
    }

    // Count faces for each object
    for (int i = 0; i < Image.Depth(); i++) {
        for (int j = 0; j < Image.Width(); j++) {
            for (int k = 0; k < Image.Height(); k++) {
                int currentLabel = output.getPixel(i, j, k);  // Use labels from the output image
                if (currentLabel > 0) {
                    for (int di = -1; di <= 1; di++) {
                        for (int dj = -1; dj <= 1; dj++) {
                            for (int dk = -1; dk <= 1; dk++) {
                                if (di == 0 && dj == 0 && dk == 0) {
                                    continue;
                                }
                                int ni = i + di;
                                int nj = j + dj;
                                int nk = k + dk;
                                if (ni >= 0 && ni < Image.Depth() && nj >= 0 && nj < Image.Width() && nk >= 0 && nk < Image.Height()) {
                                    int neighborLabel = output.getPixel(ni, nj, nk);  // Use labels from the output image
                                    if (neighborLabel != currentLabel) {
                                        facesCount[currentLabel]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Print faces count for each object
    for (int i = 1; i <= count; i++) {
        std::cout << "Object " << i << " has " << facesCount[i] << " faces." << std::endl;
    }
}
void Image3D::highlightObject(const Image3D& input, Image3D& output, int objectToHighlight) const {
    if (input.depth != output.depth || input.width != output.width || input.height != output.height) {
        
        return;
    }

    for (int i = 0; i < input.depth; i++) {
        for (int j = 0; j < input.width; j++) {
            for (int k = 0; k < input.height; k++) {
                if (input.getPixel(i, j, k) == objectToHighlight) {
                    //std::cout<<objectToHighlight<<std::endl;

                    output.setPixel(i, j, k, objectToHighlight);
                } //else {
//                    output.setPixel(i, j, k, 0);  // Otros píxeles quedan en 0 (blanco)
//                }
            }
        }
    }
}
void Image3D::Objects3DinPlane(/*const*/ Image3D& input, Image3D& output, int plano) const {
    if (input.depth != output.depth || input.width != output.width || input.height != output.height) {
        return; // Las dimensiones de entrada y salida no coinciden
    }
    TPZVec<int> objectvalues; /*TPZVec<double>*/
    for (int depthIndex = 0; depthIndex < input.depth; depthIndex++) {
        for (int widthIndex = 0; widthIndex < input.width; widthIndex++) {
            for (int heightIndex = 0; heightIndex < input.height; heightIndex++) {
                if (depthIndex == plano) {
                    
                    if (input.getPixel(depthIndex, widthIndex, heightIndex) != 0) {
                        int obj=input.getPixel(depthIndex, widthIndex, heightIndex);
                        //output.setPixel(depthIndex, widthIndex, heightIndex, 50); // Píxeles en el plano y no nulos
                        
                        //input.highlightObject(input,output,pix);
                        //std::cout<<obj<<std::endl;
                        objectvalues.push_back(obj);
                    } else {
                        output.setPixel(depthIndex, widthIndex, heightIndex, 1); // Píxeles en el plano y nulos
                    }
                } else {
                    output.setPixel(depthIndex, widthIndex, heightIndex, 0);  // Otros píxeles quedan en 0 (blanco)
                }
            }
        }
    }
    //std::cout<<objectvalues[0]<<std::endl;
    //std::cout<<objectvalues[1]<<std::endl;
    //std::cout<<objectvalues[2]<<std::endl;
    //std::cout<<objectvalues.size()<<std::endl;
    //std::cout<<objectvalues[4]<<std::endl;
    std::cout<<"Objetos que cortan el plano: "<<objectvalues.size()<<std::endl;
    for(int objeto=0;objeto<objectvalues.size();objeto++){
        input.highlightObject(input,output,objectvalues[objeto]);
        std::cout<<"Objeto: "<<objeto+1<<" Label: "<<objectvalues[objeto]<< " Numero de pixels: "<<input.getPixelsInObject(objeto)<<std::endl;
    }
}

//void Image3D::Objects3DinPlane( Image3D& input, Image3D& output, int plano) const {
//    if (input.depth != output.depth || input.width != output.width || input.height != output.height) {
//
//        return;
//    }
//    //Image3D input1("objects",input.Depth(), input.Width(), input.Height());
//
//    for (int i = 0; i < input.depth; i++) {
//        for (int j = 0; j < input.width; j++) {
//            for (int k = 0; k < input.height; k++) {
//                if(i==plano && input.getPixel(i, j, k) == 0 ){
//                    output.setPixel(i,j,k,10000);
//                }
//                if(i==plano && input.getPixel(i, j, k) != 0 ){
//                    output.setPixel(i,j,k,50);
//                    //input.highlightObject(input,output,input.getPixel(i, j, k));
//
//                }
//                else {
//                    output.setPixel(i, j, k, 0);  // Otros píxeles quedan en 0 (blanco)
//                }
//            }
//        }
//    }
////    for (int i = 0; i < input.depth; i++) {
////        if (i != plano) {
////            highlightObject(input, output, i);
////        }
////    }
//
//}

//std::vector<std::pair<int, int>> obtenerObjetosYPixeles(const Image3D& ordered, int numColors) {
//    std::vector<std::pair<int, int>> objetoYPixeles;
//
//    for (int i = 1; i < numColors; i++) {
//        int numPixeles = ordered.getPixelsInObject(i);
//        objetoYPixeles.push_back(std::make_pair(i, numPixeles));
//    }
//
//    return objetoYPixeles;
//}
