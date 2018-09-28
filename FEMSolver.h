#ifndef FEMSOLVE_H
#define FEMSOLVE_H
#include"FEM.h"
class FEMSolver : public FEM
{
private:
    Vector V;
    double theta;
    double alpha;
    double intervalEnd;
    double dt;
    void initialSetup(Mesh&mesh);
public:
    std::vector<Vector> solutions;
    std::vector<Vector> velocities;
    std::vector<Vector> adaptedSolutions;
    std::vector<Vector> adaptedVelocities;
    std::vector<Mesh> adaptedMeshs;
    int numberOfTimeSteps;
    int timestep;
    void init(Mesh&mesh);
    void solve(Mesh&mesh);
    void adaptMeshToSolutions(Mesh&mesh);
    Vector getSolution();
    Mesh &getAdaptMesh();
    Vector getAdaptSolution();
};
#endif
