#include"GLnixAPP.h"
#include"GeometryGen.h"
#include"MathHelper.h"
#include"FEMSolver.h"
// Instead of typing std::vector<Type> everytime
// we use a STL vector, we can just write
// vector<Type> after using the following line 
// of code:
using std::vector;

class FEMDriver : public GLnixAPP
{
public:
	FEMDriver(); /* Constructor */
	bool Init(); /* Initialisation routine */
	void UpdateScene(float dt); /* Calls Solve() from FEMSolver */
	void RedrawTheWindow();
        
        // Functions used to overide mouse input
	void OnMouseDown(XButtonEvent btn,int x, int y);
	void OnMouseUp(XButtonEvent btn,int x, int y);
	void OnMouseMove(int x, int y);
private:
        // Functions used to package data ready to
        // send over to the GPU
    	void BuildLandGridBuffers();
	void BuildMeshBuffer();
public:
        // Used in the grpahics pipeline
    	UINT VAOmesh;
    	UINT VAOLAND;
    	UINT meshIndexCount;
	UINT gridIndexCount;
	UINT BUFFERS[2];
	
        // Used to display the solution
    	GeometryGenerator::MeshData meshObj;
    	GeometryGenerator geoGen;

    	Mesh mesh;
    	FEMSolver pfem;
    
	int mousex;
	int mousey;
	int but;

        // Projection matrices for displaying
        // 3D geometry
	AV4X4FLOAT viewModelMatrix;
	AV4X4FLOAT projMatrix;
        
        // Used in mouse momement overrides
	float theta;
	float phi;
	float radius;
};

