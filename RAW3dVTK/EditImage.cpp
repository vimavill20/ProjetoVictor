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
int Image3D::identifyObjects(Image3D& output) {
    // Assume this method identifies objects in the 3D image and writes the result to the output image.
    int count = 1;
    // Initialize the output image.
    for (auto &mat : output.data) {
        mat.Zero();
    }
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < height; k++) {
                if (getPixel(i,j,k) > 0 && output.getPixel(i, j, k) == 0) {
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
        std::cout << "number of objects smaller than " << (2 << 2*pair.first) << ": " << pair.second << std::endl;

    }
}
