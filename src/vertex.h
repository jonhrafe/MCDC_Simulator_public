//!  Auxiliary class. Implements basic vertices. ====================================/
/*!
*   Helper class to strore and handle vertices.
*   \author Jonathan Rafael
*   \date   July 2016
*   \version   0.2
*===================================================================================*/

#ifndef VERTEX_H
#define VERTEX_H

/*! \class Vertex
 *  \brief Vertex of a 3d poly
 */
class Vertex{
public:
    Vertex();
    Vertex(const double& x,const double& y,const double& z);
    unsigned index;
    double   points[3];
    double operator ()(unsigned i){ return points[i];}
};



#endif // VERTEX_H
