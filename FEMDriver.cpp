#include"FEMDriver.h"
#include"myextloader.h"


PFNGLBINDBUFFERPROC GLnix_glBindBuffer;
PFNGLGENBUFFERSPROC GLnix_glGenBuffers;
PFNGLDELETEBUFFERSPROC GLnix_glDeleteBuffer;
PFNGLBUFFERDATAPROC GLnix_glBufferData;
PFNGLBUFFERSUBDATAPROC GLnix_glBufferSubData;
PFNGLMAPBUFFERPROC GLnix_glMapBuffer;
PFNGLUNMAPBUFFERPROC GLnix_glUnmapBuffer;
PFNGLBINDVERTEXARRAYPROC GLnix_glBindVertexArray;
PFNGLGENVERTEXARRAYSPROC GLnix_glGenVertexArrays;

//
FEMDriver::FEMDriver()
: GLnixAPP(), phi(1.5f*MathHelper::Pi),theta(1.5f*MathHelper::Pi), radius(10)
{
    load_extension_function_pointers();
    mousex = 0;
    mousey = 0;
    AV4FLOAT r(1,1,1,1);
    AV4X4FLOAT I;
    I.diag (r);

    projMatrix = I;
    viewModelMatrix = I;

}



bool FEMDriver::Init()
{
    

    if(!GLnixAPP::Init())
        return false;

        BuildLandGridBuffers();
        BuildMeshBuffer();
   
    
    projMatrix =formProjMatrix(0.25f*MathHelper::Pi, AspectRatio(), 1.0f, 1000.0f);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    return true;
}


void FEMDriver::UpdateScene(float dt)
{


    float x = radius*sinf(phi)*cosf(theta);
    float z = radius*sinf(phi)*sinf(theta);
    float y = radius*cosf(phi);


    AV4FLOAT position(x,y,z,1.0);
    
    AV4FLOAT target(0.0,0.0,0.0,0.0);
    AV4FLOAT up(0.0,1.0,0.0,0.0);

    viewModelMatrix = formViewModelMatrix(position,target,up);

    BuildMeshBuffer();
    if(pfem.timestep == pfem.numberOfTimeSteps-1)
    {
        
    }
    else
    {
        pfem.timestep++;
    } 
    
}



void FEMDriver::RedrawTheWindow()
{
    float const aspect = AspectRatio();

    float x = radius*sinf(phi)*cosf(theta);
    float z = radius*sinf(phi)*sinf(theta);
    float y = radius*cosf(phi);



    glDrawBuffer(GL_BACK);

    glViewport(0, 0, width, height);
    glClearColor(0, 0, 0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glLoadMatrixf(projMatrix.m);

    



    


    glMatrixMode(GL_MODELVIEW);
   
    glCullFace(GL_BACK);
    glLoadMatrixf(viewModelMatrix.m);
    glColor4f(1,1,1, 1);
    GLnix_glBindVertexArray(VAOmesh);
    glDrawElements(GL_TRIANGLES, meshIndexCount, GL_UNSIGNED_INT, 0);
    


    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    glLoadMatrixf(viewModelMatrix.m);
    glColor4f(0.1,0.5,0.5,1);
    GLnix_glBindVertexArray(VAOLAND);
    glDrawElements(GL_TRIANGLES, gridIndexCount, GL_UNSIGNED_INT, 0);

    glCullFace(GL_FRONT);

    
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
    glLoadMatrixf(viewModelMatrix.m);
    glColor4f(1,1,1, 1);
    GLnix_glBindVertexArray(VAOmesh);
    glDrawElements(GL_TRIANGLES, meshIndexCount, GL_UNSIGNED_INT, 0);
    
    glXSwapBuffers(Xdisplay, glX_window_handle);
}


void FEMDriver::BuildMeshBuffer()
{

    geoGen.CreateSolutionSurfaceFromObjects(pfem.getAdaptMesh(),pfem.getAdaptSolution(),meshObj);




    meshIndexCount = meshObj.Indices.size();
    GLnix_glGenVertexArrays(1,&VAOmesh);
    GLnix_glGenBuffers(2,BUFFERS);
    GLnix_glBindVertexArray(VAOmesh); 
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);


    GLnix_glBindBuffer(GL_ARRAY_BUFFER,BUFFERS[0]);
    GLnix_glBufferData(GL_ARRAY_BUFFER,
                                  meshObj.Vertices.size()*sizeof(GLfloat)*11,
                                  &meshObj.Vertices.front(), GL_STATIC_DRAW);

    
    GLnix_glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BUFFERS[1]);
    GLnix_glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                                   meshObj.Indices.size() * sizeof(UINT), 
                                   &meshObj.Indices.front(), GL_STATIC_DRAW);
    
    glVertexPointer(3, GL_FLOAT,sizeof(GLfloat)*11, 0); 
    glNormalPointer(GL_FLOAT,sizeof(GLfloat)*11,(GLvoid*)(3*sizeof(GLfloat)));
}




void FEMDriver::BuildLandGridBuffers()
{

        GeometryGenerator::MeshData grid;
        GeometryGenerator geoGen;
        geoGen.CreateGrid(20.0f, 20.0f, 20, 20, grid);

        gridIndexCount = grid.Indices.size();
        std::vector<GeometryGenerator::Vertex> vertices(grid.Vertices.size());
        for(size_t i = 0; i < grid.Vertices.size(); i++)
        {
                AV3FLOAT p = grid.Vertices[i].Position;
                
                p.y = -2; 
                vertices[i].Position   = p;
        }

        //now we are going to sort of the normals, because we have our land
        //as an explicit function we can analytically calculate the normal vecs!

        for(size_t i = 0; i < grid.Vertices.size(); i++)
        {
                AV3FLOAT p = grid.Vertices[i].Position;
                p.x = 0;
                p.y = 1;
                p.z = 0;
                p.normalize();
                vertices[i].Normal   = p;
        }

        GLnix_glGenVertexArrays(1,&VAOLAND);
        GLnix_glGenBuffers(2,BUFFERS);
        GLnix_glBindVertexArray(VAOLAND);

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);


        GLnix_glBindBuffer(GL_ARRAY_BUFFER,BUFFERS[0]);
        GLnix_glBufferData(GL_ARRAY_BUFFER,
                                       vertices.size() * sizeof(GLfloat) * 11,
                                       &vertices.front(), GL_STATIC_DRAW);


        GLnix_glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BUFFERS[1]);
        GLnix_glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                                   gridIndexCount * sizeof(UINT), 
                                   &grid.Indices.front(), GL_STATIC_DRAW);

        glVertexPointer(3, GL_FLOAT,sizeof(GLfloat) * 11, 0);
        glNormalPointer(GL_FLOAT, sizeof(GLfloat) * 11,  (GLvoid*)(3*sizeof(GLfloat)) );


}


void FEMDriver::OnMouseDown(XButtonEvent btn,int x, int y)
{
    mousex = x;
    mousey = y;
    but = btn.button;
}
void FEMDriver::OnMouseUp(XButtonEvent btn,int x, int y)
{

}
void FEMDriver::OnMouseMove(int x, int y)
{
    if(but == 1)
    {
        // Make each pixel correspond to a quarter of a degree.
        float dx = ANTMATHConvertToRadians(0.25f*static_cast<float>(x - mousex));
        float dy = ANTMATHConvertToRadians(0.25f*static_cast<float>(y -mousey));

        // Update angles based on input to orbit camera around box.
        theta += dx;
        phi   += dy;

        // Restrict the angle mPhi.
        phi = MathHelper::Clamp(phi, 0.1f, MathHelper::Pi-0.1f);
        std::cout << "Left Button Pressed" << std::endl;
    }
    else if (but == 3)
    {
        // Make each pixel correspond to 0.2 unit in the scene.
        float dx = 0.2f*static_cast<float>(x - mousex);
        float dy = 0.2f*static_cast<float>(y - mousey);

        // Update the camera radius based on input.
        radius += dx - dy;

        // Restrict the radius.
        radius = MathHelper::Clamp(radius, 1.0f, 200.0f);
        std::cout << "Right Button Pressed" << std::endl;
    }

    else if(but ==2)
    {
        std::cout << "Middle Button Pressed" << std::endl;
        pfem.timestep = 0;
    }
    

    mousex = x;
    mousey = y;
}
