#include<iostream>
#include"FEMDriver.h"
int main()
{   
    // Instantiates a FEMDriver object
    FEMDriver femd;
    // Reads the vertices in from a file passed in as an argument
    femd.mesh.ReadVertices("VTEST3.txt");
    // Initial triangulation.
    femd.mesh.Triangulate();
    femd.mesh.Refine(1); femd.mesh.Improve();
    femd.mesh.Refine(1); femd.mesh.Improve();

    // Forms mass and stiffness matrix and sets intial condition
    femd.pfem.init(femd.mesh);

    // To be included in every MassStiff project
    femd.pfem.solve(femd.mesh);


    femd.pfem.adaptMeshToSolutions(femd.mesh);
    if( !femd.Init() )
        return 0;
    return femd.Run();
}


