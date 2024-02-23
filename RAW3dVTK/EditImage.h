#ifndef EDITIMAGE_H
#define EDITIMAGE_H

#include "pzvec.h"
// Assume Matrix is a class that represents a 2D image.
// It should have methods like getPixel(x, y) and setPixel(x, y, value).
#include "pzfmatrix.h"

class Image3D {
public:
    Image3D(int depth, int width, int height)
        : depth(depth), width(width), height(height), data(depth, TPZFMatrix<REAL>(width, height)) {}
    // Overload the assignment operator.
    Image3D& operator=(const Image3D& other) {
        if (this != &other) {
            depth = other.depth;
            width = other.width;
            height = other.height;
            data = other.data;  // Assumes Matrix class has a working copy assignment operator.
        }
        return *this;
    }

    int getPixel(int x, int y, int z) const {
        return data[z](x, y);
    }

    void setPixel(int x, int y, int z, int value) {
        data[z].PutVal(x, y, value);
    }

    void applyFilter(const TPZFMatrix<REAL>& filter) {
        for (TPZFMatrix<REAL>& slice : data) {
            slice = filter;
        }
    }

    void identifyObjects(Image3D& output);

private:
    int depth;
    int width;
    int height;
    TPZVec<TPZFMatrix<REAL>> data;

    static void InsertSurroundingPixels(const Image3D& input, Image3D& output, int x, int y, int z, int value);
};

#endif  // EDITIMAGE_H