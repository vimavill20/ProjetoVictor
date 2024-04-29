#ifndef EDITIMAGE_H
#define EDITIMAGE_H

#include "pzvec.h"
// Assume Matrix is a class that represents a 2D image.
// It should have methods like getPixel(x, y) and setPixel(x, y, value).
#include "pzfmatrix.h"

class Image3D {
public:
    Image3D(const std::string &varname, int depth, int width, int height)
        : depth(depth), width(width), height(height), varname(varname), data(depth, TPZFMatrix<REAL>(width, height)) {}

    Image3D(const Image3D& other) {
        varname = other.varname;
        depth = other.depth;
        width = other.width;
        height = other.height;
        data = other.data;  // Assumes Matrix class has a working copy constructor.
    }

    Image3D(const std::string &varname, TPZVec<TPZFMatrix<REAL> > &input) : varname(varname){
        depth = input.size();
        width = input[0].Rows();
        height = input[0].Cols();
        data = input;
    }
    // Overload the assignment operator.
    Image3D& operator=(const Image3D& other) {
        if (this != &other) {
            varname = other.varname;
            depth = other.depth;
            width = other.width;
            height = other.height;
            data = other.data;  // Assumes Matrix class has a working copy assignment operator.
        }
        return *this;
    }

    int Depth() const { return depth; }
    int Width() const { return width; }
    int Height() const { return height; }

    std::string VarName() const { return varname; }

    void SetVarName(const std::string &varname) {
        this->varname = varname;
    }

    int getPixel(int x, int y, int z) const {
        return data[x](y, z);
    }

    void setPixel(int x, int y, int z, int value) {
        data[x].PutVal(y, z, value);
    }

    void applyFilter(const TPZFMatrix<REAL>& filter) {
        for (TPZFMatrix<REAL>& slice : data) {
            slice = filter;
        }
    }

    // returns the number of identified objects.
    int identifyObjects(Image3D& output)const;

    /// order the objects by size
    void orderObjectsBySize(Image3D& output, int numcolors);
    void SegmentVugFracture(Image3D& output, int numcolors,const std::string& filename);
    
    void countFacesByObject(std::map<int, int>& facesCount) const;
    int getPixelsInObject(int label)const;
    TPZVec<double> obtenerObjetosYPixeles(const Image3D& ordered, const int numColors);
    void countFacesByObject(std::map<int, int>& facesCount, const Image3D& Image) const;
    void highlightObject(const Image3D& input, Image3D& output, int objectToHighlight) const;
    void Objects3DinPlane( Image3D& input, Image3D& output, int plano) const;

private:
    int depth;
    int width;
    int height;
    std::string varname;
    TPZVec<TPZFMatrix<REAL>> data;
    
    static void InsertSurroundingPixels(const Image3D& input, Image3D& output, int x, int y, int z, int value);
};
#endif  // EDITIMAGE_H
