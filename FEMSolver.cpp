#include"FEMSolver.h"
void FEMSolver::init(Mesh& mesh)
{  
    // This constructor deals with
    // memory allocation.
    nodeCount = mesh.mVertices.size();
    elementCount = mesh.mTriangles.size();
    GlobalMatrix.setMatrix(Matrix(nodeCount,nodeCount));
    GlobalMass.setMatrix(Matrix(nodeCount,nodeCount));
    GlobalStiff.setMatrix(Matrix(nodeCount,nodeCount));
    b.setVector(Vector(nodeCount));
    U.setVector(Vector(nodeCount));
    
    // Formation of global mass and stiffness matrices
    for(int l=0;l<elementCount;l++)
    {
        AV3INT tempTri  = mesh.mTriangles[l];
        double x1 = mesh.mVertices[tempTri[0]].Position.x;
        double y1 = mesh.mVertices[tempTri[0]].Position.y;
        double x2 = mesh.mVertices[tempTri[1]].Position.x;
        double y2 = mesh.mVertices[tempTri[1]].Position.y;
        double x3 = mesh.mVertices[tempTri[2]].Position.x;
        double y3 = mesh.mVertices[tempTri[2]].Position.y;

        Matrix emm = elementMass(x1,x2,x3,y1,y2,y3);
        Matrix esm = elementStiff(x1,x2,x3,y1,y2,y3);

        for(int i=0;i<3;i++)
        {
            int i1=tempTri[i];
            for(int j=0;j<3;j++)
            {
                int j1 = tempTri[j];
                GlobalStiff(i1,j1) = GlobalStiff(i1,j1) + esm(i,j);
                GlobalMass(i1,j1) = GlobalMass(i1,j1) + emm(i,j);
            }    
        }
    }
    initialSetup(mesh);
}

// THE USER OF THE CLASS WOULD BE ADVISED TO
// MODIFY THESE FUNCTIONS.

void FEMSolver::initialSetup(Mesh&mesh)
{
    timestep = 0;
    intervalEnd = 3;
    numberOfTimeSteps = 300;
    dt = (1.0*intervalEnd/numberOfTimeSteps);

    // Used to store the derivative of the solution.
    V.setVector(Vector(nodeCount));

    // User defined constants.
    theta = 1; /* Controls the numerical scheme */
    alpha = 1; /* Wave equation parameter */

    // Initial condition for the wave equation:
    //      1 on the boundary, zero everywhere else.
    for(int i=0;i<mesh.mVertices.size();i++)
        if(sqrt(mesh.mVertices[i].Position.x*mesh.mVertices[i].Position.x + mesh.mVertices[i].Position.x*mesh.mVertices[i].Position.x) < 0.2)
            U(i)=0.5;
}

// Calling solve will essentially solve the problem.a
void FEMSolver::solve(Mesh&mesh)
{
    double a = -1.0*dt*dt*theta*(1-theta)/alpha;
    // We solve the wave equation using the Crank-Nicolson scheme.
    for(int i=0;i<numberOfTimeSteps;i++){
    Vector temp = (GlobalMass+a*GlobalStiff)*U + theta*dt*GlobalMass*V;
    Vector Unew = temp/(GlobalMass+((dt*theta*dt*theta)/alpha)*GlobalStiff);
    temp = GlobalStiff*(theta*Unew+(1-theta)*U);
    Vector Vnew = ((GlobalMass*V) - (dt/alpha)*temp)/GlobalMass;
    U = Unew;
    V = Vnew;
    solutions.push_back(U);
    velocities.push_back(V);
    std::cout << "Solution number " << i << " obtained" << std::endl;
    }  
}

void FEMSolver::adaptMeshToSolutions(Mesh&mesh)
{
    for(int h=1;h<solutions.size();h++)
    {
    Vector oldsol = solutions[h-1];
        Vector oldvel = velocities[h-1];
    Mesh tMesh = mesh.AdaptNew(solutions[h],0.25,oldsol,oldvel);

    nodeCount = tMesh.mVertices.size();
    elementCount = tMesh.mTriangles.size();

    Matrix gmm(nodeCount,nodeCount);
    Matrix gsm(nodeCount,nodeCount);

    // Formation of global mass and stiffness matrices
    for(int l=0;l<elementCount;l++)
    {
        AV3INT tempTri  = tMesh.mTriangles[l];
        double x1 = tMesh.mVertices[tempTri[0]].Position.x;
        double y1 = tMesh.mVertices[tempTri[0]].Position.y;
        double x2 = tMesh.mVertices[tempTri[1]].Position.x;
        double y2 = tMesh.mVertices[tempTri[1]].Position.y;
        double x3 = tMesh.mVertices[tempTri[2]].Position.x;
        double y3 = tMesh.mVertices[tempTri[2]].Position.y;

        Matrix emm = elementMass(x1,x2,x3,y1,y2,y3);
        Matrix esm = elementStiff(x1,x2,x3,y1,y2,y3);

        for(int i=0;i<3;i++)
        {
            int i1=tempTri[i];
            for(int j=0;j<3;j++)
            {
                int j1 = tempTri[j];
                gsm(i1,j1) = gsm(i1,j1) + esm(i,j);
                gmm(i1,j1) = gmm(i1,j1) + emm(i,j);
            }    
        }
    }



    Vector temp = gmm*oldsol + dt*gmm*oldvel;
    Vector Unew = temp/(gmm+((dt*dt)/alpha)*gsm);
    temp = gsm*Unew;
    Vector Vnew = ((gmm*oldvel) - (dt/alpha)*temp)/gmm;

    if(norm(Unew) < 1.5*norm(oldsol)) 
    {
    oldsol = Unew;
        oldvel = Vnew;
        adaptedSolutions.push_back(oldsol);
        adaptedVelocities.push_back(oldvel);
        adaptedMeshs.push_back(tMesh);
    }
        else
    {
            adaptedSolutions.push_back(solutions[h]);
            adaptedVelocities.push_back(velocities[h]);
            adaptedMeshs.push_back(mesh);
    }


      std::cout << "adaptive stage: " << h << std::endl;



    }
    std::cout << adaptedSolutions.size() << " " << adaptedMeshs.size() << std::endl;
}
Vector FEMSolver::getSolution()
{
    if(timestep == numberOfTimeSteps-1)
    return *(solutions.end()-1);
    return solutions[timestep];
}

Mesh& FEMSolver::getAdaptMesh()
{
    if(timestep == numberOfTimeSteps-1)
    return *(adaptedMeshs.end()-1);
    return adaptedMeshs[timestep];
}
Vector FEMSolver::getAdaptSolution()
{
    if(timestep == numberOfTimeSteps-1)
    return *(adaptedSolutions.end()-1);
    return adaptedSolutions[timestep];
}



