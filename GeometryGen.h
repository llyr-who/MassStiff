#ifndef GEOMETRYGENERATOR_H
#define GEOMETRYGENERATOR_H
#include "antmath.h"
#include "Mesh.h"
#include "Vector.hpp"
typedef unsigned int UINT;

class GeometryGenerator
{
public:
	struct Vertex
	{

	    AV3FLOAT Position;
		AV3FLOAT Normal;
		AV3FLOAT TangentU;
		AV2FLOAT TexC;
		Vertex(){}
		Vertex(const AV3FLOAT& p, const AV3FLOAT& n, const AV3FLOAT& t, const AV2FLOAT& uv)
			: Position(p), Normal(n), TangentU(t), TexC(uv){}
		Vertex(
			float px, float py, float pz,
			float nx, float ny, float nz,
			float tx, float ty, float tz,
			float u, float v)
			: Position(px,py,pz), Normal(nx,ny,nz),
			  TangentU(tx, ty, tz), TexC(u,v){}
	};

	struct MeshData
	{
		std::vector<Vertex> Vertices;
		std::vector<UINT> Indices;
	};
	void CreateBox(float width, float height, float depth, MeshData& meshData);
	void CreateSphere(float radius, UINT sliceCount, UINT stackCount, MeshData& meshData);
	void CreateGrid(float width, float depth, UINT m, UINT n, MeshData& meshData);
	
	void CreateMeshFromObject(Mesh& mesh, MeshData& meshData);
    void CreateSolutionSurfaceFromObjects(Mesh& mesh, Vector solution, MeshData& meshData);
	void Subdivide(MeshData& meshData);
};
#endif // GEOMETRYGENERATOR_H
