#include "vertex.h"


Vertex::Vertex()
{
    points[0]=points[1]=points[2]=0;
}

Vertex::Vertex(const double &x, const double &y, const double &z)
{
    points[0]= x;
    points[1]= y;
    points[2]= z;
}
