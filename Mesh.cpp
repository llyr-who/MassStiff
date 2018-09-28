#include<string>
#include<fstream>
#include<iostream>
#include "Mesh.h"
#include<algorithm>

/*
 * ReadVertices: Takes in a .txt file and stores
 * the vertices in mVertices.
 *
 * Parameters: filename (.txt file)
 */
int Mesh::isEar(size_t n, const std::vector<Vertex> &p)
{
    return (isVertexInterior(n, p) && !isEdgeIntersect(n, p));
}

int Mesh::isEdgeIntersect(size_t n, const std::vector<Vertex> &p)
{
    std::vector<Vertex> a;

    for (size_t i = 0; i < p.size(); i++)
        if (i != n)
            a.push_back(p[i]);

    size_t c = 0, cnt = a.size(), prev = (cnt + (n - 1)) % cnt, next = n % cnt;

    Vertex v1 = a[prev], v2 = a[next];

    for (size_t i = 0, j = cnt - 1; i < cnt; j = i++)
    {
        if (prev == i || prev == j || next == i || next == j)
            continue;

        Vertex v4 = a[j], v3 = a[i];

        float d = ((v4.Position.y - v3.Position.y) * (v2.Position.x - v1.Position.x)) - ((v4.Position.x - v3.Position.x) * (v2.Position.y - v1.Position.y));

        if (!d)
            continue;

        float ua = (((v4.Position.x - v3.Position.x) * (v1.Position.y - v3.Position.y)) - ((v4.Position.y - v3.Position.y) * (v1.Position.x - v3.Position.x))) / d;
        float ub = (((v2.Position.x - v1.Position.x) * (v1.Position.y - v3.Position.y)) - ((v2.Position.y - v1.Position.y) * (v1.Position.x - v3.Position.x))) / d;

        if (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1)
        {
            c = 1;
            break;
        }
    }

    return c;
}

int Mesh::isVertexInterior(size_t n, const std::vector<Vertex> &p)
{
    Vertex v = p[n];
    std::vector<Vertex> a;

    for (size_t i = 0; i < p.size(); i++)
        if (i != n)
            a.push_back(p[i]);

    int c = 1;

    for (size_t i = 0, j = a.size() - 1; i < a.size(); j = i++)
    {
        if ((((a[i].Position.y <= v.Position.y) && (v.Position.y < a[j].Position.y)) || ((a[j].Position.y <= v.Position.y) && (v.Position.y < a[i].Position.y))) && (v.Position.x >(a[j].Position.x - a[i].Position.x) * (v.Position.y - a[i].Position.y) / (a[j].Position.y - a[i].Position.y) + a[i].Position.x))
            c = !c;
    }

    return c;
}

void Mesh::ReadVertices(const char* filename)
{
    std::string vertexstring;
    AV2FLOAT tempPosition;
    int vnum = 0;
    std::ifstream VFILE (filename);
    if(VFILE.is_open())
        {
        while(VFILE >>tempPosition)
        {
            mVertices.push_back(Vertex(vnum++,tempPosition,0,1));
        }
        
        VFILE.close();
        }
    else
        {
        std::cout << "Unable to open file";
        }

    for(int i=0;i<mVertices.size();i++)
        {
        mEdges.push_back(Edge(i,i+1,1));
        }
        mEdges[mVertices.size()-1].vnum2 = 0;
}

void Mesh::PrintVertices()
{

    std::vector<Vertex>::iterator id = mVertices.begin();
    while (id != mVertices.end())
    {
        std::cout << id->VertexNumber << std::endl;
        std::cout << id->Position.x << " " << id->Position.y << " " <<  std::endl;
        id++;
    }
    std::vector<AV3INT>::iterator i = mTriangles.begin();
    while (i != mTriangles.end())
    {
        std::cout << i->x << " " << i->y << " " <<  i->z <<  std::endl;
        i++;
    }

    std::vector<Edge>::iterator ie = mEdges.begin();
    while (ie != mEdges.end())
    {
        std::cout << ie->vnum1 << " " << ie->vnum2<<  std::endl;
        ie++;
    }
}

double Mesh::areaOfTriangle(AV3INT& triangle)
{
    return (mVertices[triangle.y].Position.x-mVertices[triangle.x].Position.x)*(mVertices[triangle.z].Position.y-mVertices[triangle.x].Position.y) - (mVertices[triangle.z].Position.x-mVertices[triangle.x].Position.x)*(mVertices[triangle.y].Position.y-mVertices[triangle.x].Position.y);
}

double Mesh::smallestAngle(AV3INT& triangle)
{
    float min = 180;
    float temp3 = 170;
 
    for(int i=4;i>1;i--)
    {
        float x1 = mVertices[triangle[(i+1)%3]].Position.x - mVertices[triangle[i%3]].Position.x;
        float x2 = mVertices[triangle[(i-1)%3]].Position.x - mVertices[triangle[i%3]].Position.x;
        float y1 = mVertices[triangle[(i+1)%3]].Position.y - mVertices[triangle[i%3]].Position.y;
        float y2 = mVertices[triangle[(i-1)%3]].Position.y - mVertices[triangle[i%3]].Position.y;
   
    
        float dot = x1*x2 + y1*y2;
        float det = x1*y2 - y1*x2;    
        temp3 = atan2(det, dot)*(180.0/AL_PI);
  
        if(temp3<min)min=temp3;
    }
  
    return min;

}

void Mesh::fixTri(AV3INT& t)
{
    float area= areaOfTriangle(t);
    if(area<0)
    {
        int tempind = t.x;
        t.x= t.y;
        t.y = tempind;
    }
}


void Mesh::Adapt(Vector&solution,double thres)
{
    int vertNum = mVertices.size();
    std::vector<Vertex> toBeAdded;
    int trilen = mTriangles.size();
    for(int i=0;i<trilen;i++)
    {

        AV3INT triangle = mTriangles[i];
        AV3INT tri1,tri2,tri3,tri4;
        std::vector<double> solutionTri;
        solutionTri.push_back(solution(triangle.x));
        solutionTri.push_back(solution(triangle.y));
        solutionTri.push_back(solution(triangle.z));
        std::sort(solutionTri.begin(),solutionTri.end());
        if(fabs(solutionTri[0]-solutionTri[2])>thres)
        {
           Vertex v0 = mVertices[ triangle[0] ];
           Vertex v1 = mVertices[ triangle[1] ];
           Vertex v2 = mVertices[ triangle[2] ];
           Vertex m0, m1, m2;
       
           m0.Position = AV2FLOAT(
                   0.5f*(v0.Position.x + v1.Position.x),
                   0.5f*(v0.Position.y + v1.Position.y));
          

           m1.Position = AV2FLOAT(
                   0.5f*(v1.Position.x + v2.Position.x),
                   0.5f*(v1.Position.y + v2.Position.y));
           

           m2.Position = AV2FLOAT(
                   0.5f*(v0.Position.x + v2.Position.x),
                   0.5f*(v0.Position.y + v2.Position.y));
           

           
           int m0add=1;
           int m1add=1;
           int m2add=1;
           int size = toBeAdded.size();
           if(size>0)
           {
                for(int i =0;i<size;i++)
                {
                    if(m0.distance(toBeAdded[i])<0.01)
                    {
                        m0add=0;
                        m0.VertexNumber = toBeAdded[i].VertexNumber;
                        continue;
                    }

                    if(m1.distance(toBeAdded[i])<0.01)
                    {
                        m1add=0;
                        m1.VertexNumber = toBeAdded[i].VertexNumber;
                        continue;
                    }

                    if(m2.distance(toBeAdded[i])<0.01)
                    {
                        m2add=0;
                        m2.VertexNumber = toBeAdded[i].VertexNumber;
                        continue;
                    }

                }
           } 

           if(m0add) {m0.VertexNumber=vertNum++;toBeAdded.push_back(m0);}
           if(m1add) {m1.VertexNumber=vertNum++;toBeAdded.push_back(m1);}
           if(m2add) {m2.VertexNumber=vertNum++;toBeAdded.push_back(m2);}

           tri4.setVertex(0,m0.VertexNumber);
           tri4.setVertex(1,m1.VertexNumber);
           tri4.setVertex(2,m2.VertexNumber);


           tri1.setVertex(0,tri4[0]);
           tri1.setVertex(1,tri4[2]);
           tri1.setVertex(2,triangle[0]);
           
           tri2.setVertex(0,tri4[0]);
           tri2.setVertex(1,triangle[1]);
           tri2.setVertex(2,tri4[1]);

           tri3.setVertex(0,tri4[2]);
           tri3.setVertex(1,tri4[1]);
           tri3.setVertex(2,triangle[2]);
 
           mTriangles[i]= tri1;
           mTriangles.push_back(tri4);mTriangles.push_back(tri2);mTriangles.push_back(tri3);

        }
    }
    
    for(int i=0;i<toBeAdded.size();i++)
        mVertices.push_back(toBeAdded[i]);
    for(int i=0;i<mTriangles.size();i++)
        fixTri(mTriangles[i]);
}

Mesh Mesh::AdaptNew(Vector&solution,double thres,Vector& oldsolution, Vector& velocityold)
{
    Mesh retmesh;
    retmesh.mVertices = mVertices;
    retmesh.mTriangles = mTriangles;
    retmesh.mEdges     = mEdges;
    int vertNum = retmesh.mVertices.size();
    std::vector<Vertex> toBeAdded;
    std::vector<double> interp;
    std::vector<double> interpV;
    int trilen = retmesh.mTriangles.size();
    for(int i=0;i<trilen;i++)
    {

        AV3INT triangle = retmesh.mTriangles[i];
        AV3INT tri1,tri2,tri3,tri4;
        std::vector<double> solutionTri;
        solutionTri.push_back(solution(triangle.x));
        solutionTri.push_back(solution(triangle.y));
        solutionTri.push_back(solution(triangle.z));
        std::sort(solutionTri.begin(),solutionTri.end());
        if(fabs(solutionTri[0]-solutionTri[2])>thres)
        {
           Vertex v0 = retmesh.mVertices[ triangle[0] ];
           Vertex v1 = retmesh.mVertices[ triangle[1] ];
           Vertex v2 = retmesh.mVertices[ triangle[2] ];
           Vertex m0, m1, m2;
       
           m0.Position = AV2FLOAT(
                   0.5f*(v0.Position.x + v1.Position.x),
                   0.5f*(v0.Position.y + v1.Position.y));
          

           m1.Position = AV2FLOAT(
                   0.5f*(v1.Position.x + v2.Position.x),
                   0.5f*(v1.Position.y + v2.Position.y));
           

           m2.Position = AV2FLOAT(
                   0.5f*(v0.Position.x + v2.Position.x),
                   0.5f*(v0.Position.y + v2.Position.y));
           

           
           int m0add=1;
           int m1add=1;
           int m2add=1;
           int size = toBeAdded.size();
           if(size>0)
           {
                for(int i =0;i<size;i++)
                {
                    if(m0.distance(toBeAdded[i])<0.01)
                    {
                        m0add=0;
                        m0.VertexNumber = toBeAdded[i].VertexNumber;
                        continue;
                    }

                    if(m1.distance(toBeAdded[i])<0.01)
                    {
                        m1add=0;
                        m1.VertexNumber = toBeAdded[i].VertexNumber;
                        continue;
                    }

                    if(m2.distance(toBeAdded[i])<0.01)
                    {
                        m2add=0;
                        m2.VertexNumber = toBeAdded[i].VertexNumber;
                        continue;
                    }

                }
           } 

           if(m0add) {m0.VertexNumber=vertNum++;
            toBeAdded.push_back(m0);
            interp.push_back(0.5*(oldsolution(v0.VertexNumber)+oldsolution(v1.VertexNumber)));
            interpV.push_back(0.5*(velocityold(v0.VertexNumber)+velocityold(v1.VertexNumber)));}
           if(m1add) {m1.VertexNumber=vertNum++;
            toBeAdded.push_back(m1);
            interp.push_back(0.5*(oldsolution(v1.VertexNumber)+oldsolution(v2.VertexNumber)));
            interpV.push_back(0.5*(velocityold(v1.VertexNumber)+velocityold(v2.VertexNumber)));}
           if(m2add) {m2.VertexNumber=vertNum++;
            toBeAdded.push_back(m2);
            interp.push_back(0.5*(oldsolution(v2.VertexNumber)+oldsolution(v0.VertexNumber)));
            interpV.push_back(0.5*(velocityold(v2.VertexNumber)+velocityold(v0.VertexNumber)));}

           tri4.setVertex(0,m0.VertexNumber);
           tri4.setVertex(1,m1.VertexNumber);
           tri4.setVertex(2,m2.VertexNumber);


           tri1.setVertex(0,tri4[0]);
           tri1.setVertex(1,tri4[2]);
           tri1.setVertex(2,triangle[0]);
           
           tri2.setVertex(0,tri4[0]);
           tri2.setVertex(1,triangle[1]);
           tri2.setVertex(2,tri4[1]);

           tri3.setVertex(0,tri4[2]);
           tri3.setVertex(1,tri4[1]);
           tri3.setVertex(2,triangle[2]);
           retmesh.mTriangles[i]= tri1;
           retmesh.mTriangles.push_back(tri4);retmesh.mTriangles.push_back(tri2);retmesh.mTriangles.push_back(tri3);

        }
    }
    
    for(int i=0;i<toBeAdded.size();i++)
        retmesh.mVertices.push_back(toBeAdded[i]);
    for(int i=0;i<mTriangles.size();i++)
        fixTri(retmesh.mTriangles[i]);

    Vector NEW(length(oldsolution)+interp.size());
    Vector NEWV(length(oldsolution)+interp.size());
    for(int i=0;i< length(oldsolution);i++)
    {
        NEW(i) = oldsolution(i);
        NEWV(i) = velocityold(i);
    }
    for(int i=length(oldsolution);i< interp.size();i++)
    {
        NEW(i) = interp[i-length(oldsolution)];
        NEWV(i) = interpV[i-length(oldsolution)];
    }
    oldsolution.setVector(NEW);
    velocityold.setVector(NEWV);
    return retmesh;
}



void Mesh::scanTriangles()
{
    std::vector<AV3INT> temptris;
    mTriSwapInds.clear();
    float oldSmallestAngle1;
    float oldArea1;
    for(int i=0;i<mTriangles.size();i++)
    {
        AV3INT triangle1 = mTriangles[i];
        oldSmallestAngle1 = smallestAngle(triangle1);
        oldArea1 = areaOfTriangle(triangle1);
        
        if(oldSmallestAngle1<25)
        {
            mTriSwapInds.push_back(i);
            temptris.push_back(triangle1);
        }
    }

}

void Mesh::Improve()
{

    scanTriangles();

    int contflag = 1;
    std::vector<int> swapped;
    std::vector<int>::iterator it = mTriSwapInds.begin();
    
    while(it!=mTriSwapInds.end())
    {     
        int k = *it;
        it++;
        AV3INT triangle1 = mTriangles[k];
        Edge edgearr[3];
        edgearr[0].setEndPoints(triangle1[0],triangle1[1]);
        edgearr[1].setEndPoints(triangle1[1],triangle1[2]);
        edgearr[2].setEndPoints(triangle1[2],triangle1[0]);
        for(int i=0;i<mTriangles.size();i++)
        {
            contflag =1;
            if(i==k) continue;
            if(std::find(swapped.begin(),swapped.end(),i)!=swapped.end()) continue;
            AV3INT triangle2 = mTriangles[i];
            Edge edges[3];
            edges[0].setEndPoints(triangle2[0],triangle2[1]);
            edges[1].setEndPoints(triangle2[1],triangle2[2]);
            edges[2].setEndPoints(triangle2[2],triangle2[0]);
            for(int j = 0;j<3;j++)
            {
                for(int l=0;l<3;l++)
                { 
                    if(edgearr[j].vnum1==edges[l].vnum2 && edgearr[j].vnum2==edges[l].vnum1 && contflag)
                    {
                       fixTri(triangle1);
                       fixTri(triangle2);
                       float sa1 = smallestAngle(triangle1);
                       float sa2 = smallestAngle(triangle2);
                       contflag=0;
                       float oldArea1 = areaOfTriangle(triangle1);
                       float oldArea2 = areaOfTriangle(triangle2);

                       triangle1.x = edgearr[j].vnum2;
                       triangle1.y = edgearr[(j+1)%3].vnum2;
                       triangle1.z = edges[(l+1)%3].vnum2;

                       triangle2.x = edgearr[j].vnum1;
                       triangle2.z = edgearr[(j+1)%3].vnum2;
                       triangle2.y = edges[(l+1)%3].vnum2;
                       float sa1n = smallestAngle(triangle1);
                       float sa2n = smallestAngle(triangle2);
                       float newArea1 = areaOfTriangle(triangle1);
                       float newArea2 = areaOfTriangle(triangle2);
                       fixTri(triangle1);fixTri(triangle2);
                       if(fabs(sa1n)<0.001||fabs(sa2n)<0.001||fabs(sa1n)>170||fabs(sa2n)>170) continue;
                       if(newArea1>2*oldArea1 || newArea2>2*oldArea2 || newArea1 < 0.1*oldArea1 || newArea2 < 0.1*oldArea2)continue;
                       if(0.5*(fabs(sa1+sa2))<0.5*(fabs(sa1n+sa2n)))
                       {
                            mTriangles[k] = triangle1;
                            mTriangles[i] = triangle2;
                            swapped.push_back(k);
                       }

                    }
                }
            }
        }
        
    }
}

//
void Mesh::Triangulate()
{
    std::vector<Vertex> a;
    AV3INT T;
    for (size_t i = 0; i < mVertices.size(); i++)
    {
        a.push_back(mVertices[i]);
    }
    for (size_t t = a.size() - 1, i = 0, j = 1; i < a.size(); t = i++, j = (i + 1) % a.size())
    {
        if (a.size() == 3)
        {
            T.x = a[0].VertexNumber; T.y = a[1].VertexNumber; T.z = a[2].VertexNumber;
            mTriangles.push_back(T);
            break;
        }

        if (isEar(i, a))
        {
            T.x = a[t].VertexNumber; T.y = a[i].VertexNumber; T.z = a[j].VertexNumber;

            //here we have the diagonals
            Edge d(a[j].VertexNumber,a[t].VertexNumber,0);
            mEdges.push_back(d);
            mTriangles.push_back(T);
            a.erase(a.begin() + i, a.begin() + i + 1);
            t = a.size() - 1;
            i = 0;
            j = 1;
        }
    }
    mArea=0;
    for(int i=0;i<mTriangles.size();i++)
        mArea+=areaOfTriangle(mTriangles[i]);
}




void Mesh::Refine(int refineCount)
{
    // Calculation of how many elements we are going
    // to have after triangulation. 
    int nElements = mTriangles.size();
    int nNewElements = nElements*pow(4,refineCount);

    int nVerts = mVertices.size();
    std::vector<Vertex> lVertices;
    std::vector<AV3INT> lTriangles;
    std::vector<Edge>  lEdges;
    for(int j=0;j<refineCount;j++)
    {
    //Save a copy of the input geometry
        lVertices = mVertices;
    lTriangles = mTriangles;
        lEdges = mEdges;
        mTriangles.clear();
        mVertices.clear();
        mEdges.clear();


        UINT numTris = lTriangles.size();
        for(UINT i = 0; i < numTris; ++i)
        {
                Vertex v0 = lVertices[ lTriangles[i][0] ];
                Vertex v1 = lVertices[ lTriangles[i][1] ];
                Vertex v2 = lVertices[ lTriangles[i][2] ];

                //
                // Generate the midpoints and new edges.
                // Edges are used to keep track of boundary data.

                Vertex m0, m1, m2;
                Edge e0,e1,e2,e3,e4,e5,e6,e7,e8;

                m0.Position = AV2FLOAT(
                        0.5f*(v0.Position.x + v1.Position.x),
                        0.5f*(v0.Position.y + v1.Position.y));
                m0.VertexNumber = nVerts++;

                // Here we use the diagonals to correctly allocate interior node
                // check if the segment v0v1 or the segment v1v0 is a diagonal

                Edge isD(v0.VertexNumber,v1.VertexNumber);
                auto fit = std::find(lEdges.begin(),lEdges.end(),isD);
                if(fit != lEdges.end()&& fit->boundary !=0 )
                   {
                       m0.Boundary = fit->boundary;
                   }
                else
                   {
                       m0.Boundary = 0;
                   }
                e0.setEndPoints(v0.VertexNumber,m0.VertexNumber); e0.boundary=m0.Boundary; 
                e1.setEndPoints(m0.VertexNumber,v1.VertexNumber); e1.boundary=m0.Boundary;
                                
                m1.VertexNumber = nVerts++;
                m1.Position = AV2FLOAT(
                        0.5f*(v1.Position.x + v2.Position.x),
                        0.5f*(v1.Position.y + v2.Position.y));

                isD.setEndPoints(v1.VertexNumber,v2.VertexNumber);
                fit = std::find(lEdges.begin(),lEdges.end(),isD);
                if(fit != lEdges.end()&& fit->boundary !=0 )
                   {
                       m1.Boundary = fit->boundary;
                   }
                else
                   {
                       m1.Boundary = 0;
                   }

                e2.setEndPoints(v1.VertexNumber,m1.VertexNumber); e2.boundary = m1.Boundary;
                e3.setEndPoints(m1.VertexNumber,v2.VertexNumber); e3.boundary = m1.Boundary;

                m2.VertexNumber = nVerts++;
                m2.Position = AV2FLOAT(
                        0.5f*(v0.Position.x + v2.Position.x),
                        0.5f*(v0.Position.y + v2.Position.y));

                isD.setEndPoints(v0.VertexNumber,v2.VertexNumber);
                fit = std::find(lEdges.begin(),lEdges.end(),isD);
                if(fit != lEdges.end()&& fit->boundary !=0 )
                   {
                       m2.Boundary = fit->boundary; 
                   }
                else
                   {
                       m2.Boundary = 0;
                   }
                e4.setEndPoints(v0.VertexNumber,m2.VertexNumber); e4.boundary=m2.Boundary;
                e5.setEndPoints(m2.VertexNumber,v2.VertexNumber); e5.boundary=m2.Boundary;               
                   
                e6.setEndPoints(m0.VertexNumber,m2.VertexNumber); 
                e7.setEndPoints(m2.VertexNumber,m1.VertexNumber); 
                e8.setEndPoints(m1.VertexNumber,m0.VertexNumber); 

                // Add new geometry.
                
                mVertices.push_back(v0);mVertices.push_back(v1); 
                mVertices.push_back(v2);mVertices.push_back(m0); 
                mVertices.push_back(m1); mVertices.push_back(m2); 

        mTriangles.push_back(AV3INT(i*6+0,i*6+3,i*6+5));mTriangles.push_back(AV3INT(i*6+3,i*6+4,i*6+5));
        mTriangles.push_back(AV3INT(i*6+5,i*6+4,i*6+2));mTriangles.push_back(AV3INT(i*6+3,i*6+1,i*6+4));

                mEdges.push_back(e0);mEdges.push_back(e1);
                mEdges.push_back(e2);mEdges.push_back(e3);
                mEdges.push_back(e4);mEdges.push_back(e5);
                mEdges.push_back(e6);mEdges.push_back(e7);
                mEdges.push_back(e8);
        }

    }
    //lets clear our temp buffers
    lTriangles.clear();
    lVertices.clear();
    lEdges.clear();

    //After our refinement routine we are going to have
    // repeated vertices. In order to solve this we  carry out
    // a routine that creates a container of unique vertices
    std::vector<int> vercount;
    lVertices.push_back(mVertices[mTriangles[0][0]]);
    lVertices.push_back(mVertices[mTriangles[0][1]]);
    lVertices.push_back(mVertices[mTriangles[0][2]]);
    vercount.push_back(1); vercount.push_back(1); vercount.push_back(1);
    lTriangles.push_back(AV3INT(0,1,2));
    lVertices[0].VertexNumber=0; 
    lVertices[1].VertexNumber=1;
    lVertices[2].VertexNumber=2;
    int nn=2;
    int newnode;

    for(int k=1;k<nNewElements;k++)
    {
        for(int l=0;l<3;l++)
        {
        newnode = 0;
        for(int i=0;i<nn+1;i++)
        {

            if(mVertices[mTriangles[k][l]].distance(lVertices[i])<0.001)
            {
                newnode = 1;
                vercount[i]++;
                if(k > lTriangles.size()-1)
                {
                    AV3INT appendtri;
                    appendtri.setVertex(l, i);
                    lTriangles.push_back(appendtri);
                }
                else
                {
                    lTriangles[k].setVertex(l,i);
                }
            }
        }
        if(newnode==0)
        {
            nn=nn+1;
            lVertices.push_back(mVertices[mTriangles[k][l]]);
            lVertices[nn].VertexNumber = nn;
            vercount.push_back(1);
            if(k > lTriangles.size()-1)
            {
                    AV3INT appendtri(0,0,0);
                    appendtri.setVertex(l, nn);
                    lTriangles.push_back(appendtri);
            }
            else
            {
                    lTriangles[k].setVertex(l,nn);
            }
        }

        }

    }
    mTriangles.swap(lTriangles);
    mVertices.swap(lVertices);
    mTriSwapInds.clear();

    for(int i=0;i<mTriangles.size();i++)
    {
    for(int j=0;j<3;j++)
        {
         int bdy=0;
             Edge e(mTriangles[i][j],mTriangles[i][(j+1)%3],bdy);
         int count = std::count(mEdges.begin(),mEdges.end(),e);
         if((count < 2) && (mVertices[e.vnum1].Boundary == 1 && mVertices[e.vnum2].Boundary == 1))
             {
                    e.boundary = 1;
             }
         lEdges.push_back(e);
        }
    }     
    mEdges.swap(lEdges);  
}

