#ifndef MESH_H
#define MESH_H
#include<vector>
#include"antmath.h"
#include"Vector.hpp"
typedef unsigned int UINT;
// Mesh class
class Mesh
{
public:
    std::vector<Vertex> mVertices;
    std::vector<AV3INT> mTriangles;
    std::vector<Edge> mEdges;
    std::vector<int> mTriSwapInds;
    std::vector<double> interpolated;
    double mArea;
private:
    // utility functions for triangulate
    int isEdgeIntersect(size_t n, const std::vector<Vertex> &p);
    int isVertexInterior(size_t n, const std::vector<Vertex> &p);
    int isEar(size_t n, const std::vector<Vertex> &p);
    double areaOfTriangle(AV3INT& triangle);
    double smallestAngle(AV3INT& triangle);
    void fixTri(AV3INT& t);
    void scanTriangles();
public:
    void Adapt(Vector&solution,double thres);
    Mesh AdaptNew(Vector & solution,double thres,Vector & solutionold,Vector & velocityold);
    void Improve();
    void ReadVertices(const char* filename);
    void PrintVertices();
    void Triangulate();
    void Refine(int refineCount);
};
#endif 
