#include "EditImage.h"

void Image3D::InsertSurroundingPixels(const Image3D &input, Image3D &output, int x, int y, int z, int value) {
    using Pixel = std::tuple<int, int, int>;
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
                    if (x + i >= 0 && x + i < input.width && y + j >= 0 && y + j < input.height && z + k >= 0 && z + k < input.depth) {
                        if(visited.count({x + i, y + j, z + k})) {
                            continue;
                        }
                        if (input.getPixel(x + i, y + j, z + k) > 0) {
                            output.setPixel(x+i, y+j, z+k, value);
                            toVisit.insert({x + i, y + j, z + k});
                        }
                        visited.insert({x + i, y + j, z + k});
                    }
                }
            }
        }
    }
}
void Image3D::identifyObjects(Image3D& output) {
    // Assume this method identifies objects in the 3D image and writes the result to the output image.
    int count = 1;
    // Initialize the output image.
    for (auto &mat : output.data) {
        mat.Zero();
    }
    for (int i = 0; i < depth; i++) {
        for (int j = 0; i < width; j++) {
            for (int k = 0; k < height; k++) {
                if (getPixel(i,j,k) > 0 && output.getPixel(i, j, k) == 0) {
                    InsertSurroundingPixels(*this, output, i, j, k, count);
                }
            }
        }
        count++;
    }
}