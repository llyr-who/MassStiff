#include<iostream>
#include "FEM.h"

//

//Private methods

Matrix FEM::elementMass(double x1, double x2, double x3, double y1, double y2, double y3)
{
   double area = 0.5*((x3-x1)*(y1-y2)-(y3-y1)*(x1-x2));

   return (area/12)*(eye(3)+ones(3));
}
//
Matrix FEM::elementStiff(double x1, double x2, double x3, double y1, double y2, double y3)
{
   double area = 2*((x3-x1)*(y1-y2)-(y3-y1)*(x1-x2));

   Matrix X(2,3);
   X(0,0) = (y2-y3); X(0,1)=(y3-y1); X(0,2)=(y1-y2);
   X(1,0) = -(x2-x3); X(1,1)=-(x3-x1); X(1,2)=-(x1-x2);
   Matrix Xt = X.transpose();
   Matrix R= (1/area)*Xt*X;
   return R;
}


