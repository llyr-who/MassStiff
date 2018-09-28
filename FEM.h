#ifndef FEM_H
#define FEM_H
#include"Mesh.h"
#include "Matrix.hpp"

class FEM
{
protected:
    Matrix elementMass(double x1,double x2, double x3,double y1,double y2, double y3);
    Matrix elementStiff(double x1,double x2, double x3,double y1,double y2, double y3);
    int nodeCount;
    int elementCount;
    Matrix GlobalMatrix;
    Matrix GlobalMass;
    Matrix GlobalStiff;
    Vector U;
    Vector b;
};


#endif
