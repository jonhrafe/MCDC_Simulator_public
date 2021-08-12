//!  Auxiliary class. Implements trangular barriers. ===================================/
/*!
*   Helper class to strore and handle trangular barriers.
*   \author Jonathan Rafael
*   \date   July 2016
*   \version   0.2
*=====================================================================================*/
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vertex.h"
#include "collision.h"
#include "walker.h"
#include "Eigen/Core"

class Triangle{
public:
    unsigned        index;
    Vertex          *vertices;
    Eigen::Array3i  indexes;
    Eigen::Vector3d normal;

    //Collision sphere
    Eigen::Vector3d center;
    double radius;

    Triangle();
    Triangle(Vertex* vertices, unsigned index);

    void getVertex(const unsigned i, Eigen::Vector3d &v);
    void getNormal(Eigen::Vector3d &normal);
    void saveNormalAndAuxInfo();
    bool rayIntersects(const Eigen::Vector3d &ray_origin, const Eigen::Vector3d &step, double &t);
    void stepIntersects_MT(Walker& walker,const Eigen::Vector3d &step, const double &max_length, Collision &colision);
    void stepIntersects_MT_limits(const Eigen::Vector3d &ray_origin,const Eigen::Vector3d &step, const double &max_length, Collision &colision,
    const Eigen::Vector3d &limits_mod, double limit_x,double limit_y,double limit_z);

    bool rayIntersects_MT(const Eigen::Vector3d & ray_origin, const Eigen::Vector3d &step, double &u, double &v, double &t);

    //Returns the minimum distance from the point p to the triangle
    double minDistance(const Eigen::Vector3d p);
    //double minDistancePrecise(const Eigen::Vector3d p);

};

#endif // TRIANGLE_H
