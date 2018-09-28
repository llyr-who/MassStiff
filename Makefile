all: MassStiff

#For debugging
OPT=-g -Wall
LDFLAGS=-lm -lX11 -lGL -lGLU -lXext -lXrender
OBJECTS=FEMDriver.o FEMSolver.o Mesh.o GLnixAPP.o myextloader.h GeometryGen.o MathHelper.o antmath.o Timer.o Matrix.o Vector.o FEM.o Exception.o
		
#For optimistaion
#OPT=-O

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
MassStiff:	MassStiff.cpp ${OBJECTS}
		g++ ${OPT} ${OBJECTS} -o MassStiff MassStiff.cpp ${LDFLAGS}
clean:
	rm -f *.o *~ MassStiff
