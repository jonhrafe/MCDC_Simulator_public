#include "cylinder.h"
#include "constants.h"
#include "Eigen/Dense"
#include <iostream>

using namespace Eigen;

int Cylinder::count = 0;
Cylinder::Cylinder()
{
    id = count++;
}

Cylinder::~Cylinder()
{
    count--;
}

Cylinder::Cylinder(const Cylinder &cyl)
{

    D = cyl.D;
    Q = cyl.Q;
    P = cyl.P;
    radius = cyl.radius;
    id = count++;

}

bool Cylinder::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    //Origin of the ray
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d m = O - P;

    //minimum distance to the cylinder axis.
    double distance_to_cilinder = (D.cross(-m)).norm();
    double d_ = distance_to_cilinder - radius;

    //If the minimum distance from the walker to the cylinder is more than
    // the actual step size, we can discard this collision.
    if(d_> EPS_VAL){
        if(d_ > step_lenght+barrier_tickness){
            return false;
        }
    }

    double md = m.dot(D);
    double nd = step.dot(D);
    double nn = 1.0;
    double mm = m.dot(m);
    double a  = nn - nd*nd;
    double k  = mm - radius*radius;
    double c  = k  - md*md;


    //Parallel trajectory // WARNING: Check this stuff
    if(fabs(a) < 1e-5 && fabs(c)<barrier_tickness){
        colision.type = Collision::near;
        colision.rn = c;
        colision.obstacle_ind = id;
        return true;
    }

    double mn = m.dot(step);
    double b = mn - nd*md;
    double discr = b*b - a*c;

    //No real roots
    if(discr < 0.0){
        colision.type = Collision::null;
        return false;
    }

    //if we arrived here we need to compute the quadratic equation.
    return handleCollition(walker,colision,step,a,b,c,discr,step_lenght);

}

inline bool Cylinder::handleCollition(Walker& walker, Collision &colision, Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length){

    double t1 = (-b - sqrt(discr))/a;

    double t2 = (-b + sqrt(discr))/a;


    //if we are completely sure that no collision happened
    if( ( (t1 < 0.0) || (t1 > step_length+barrier_tickness) ) && ( (t2 < 0.0) || (t2 > step_length+barrier_tickness)) ){
        colision.type = Collision::null;
        return false;
    }

    //WARNING: Cuidar este patch
    // Implementa Percolacion
    if(percolation>0.0){
        double _percolation_ ((double)rand()/RAND_MAX);

        if( percolation - _percolation_ > EPS_VAL ){
            count_perc_crossings++;
            return false;
        }
    }

    // a spin that's bouncing ignores collision at 0 (is in a wall)
    if(walker.status == Walker::bouncing){

        //if the collision are too close or negative.
        if( ( (t1 < EPS_VAL) || (t1 > step_length+barrier_tickness)) && (( t2 < EPS_VAL) || (t2 > step_length+barrier_tickness)) ){
            colision.type = Collision::null;
            return false;
        }

        if( t1 >= EPS_VAL && t1 < t2)
            colision.t = fmin(t1,step_length);
        else
            colision.t = fmin(t2,step_length);
    }
    else{
        if( t1>0.0 && t1 <t2)
            colision.t = fmin(t1,step_length);
        else
            colision.t = fmin(t2,step_length);
    }

    colision.type = Collision::hit;
    colision.obstacle_ind = id;

    if(c<-1e-10){
        colision.col_location = Collision::inside;
        walker.in_obj_index = id;
    }
    else if(c>1e-10){
        colision.col_location = Collision::outside;
    }
    else{
        colision.col_location = Collision::unknown;
    }

    colision.rn = c;

    colision.colision_point = walker.pos_v + colision.t*step;

    if (fabs(a) < EPS_VAL){
        colision.col_location = Collision::on_edge;
        colision.bounced_direction = -step;
    }
    else{
        Eigen::Vector3d V = colision.colision_point - P;
        double v = V.dot(D);
        Eigen::Vector3d axis_point = P + v*D;
        //Normal point
        Eigen::Vector3d normal = (colision.colision_point-axis_point).normalized();

        Eigen::Vector3d temp_step = step;
        elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);

        colision.bounced_direction = temp_step.normalized();

    }

    return true;

}

double Cylinder::minDistance(Walker &w){

    //Origin of the ray
    Vector3d O;
    w.getVoxelPosition(O);
    Vector3d m = O - P;
    // minimum distance to the cylinder axis.
    double distance_to_cylinder = (D.cross(-m)).norm();

    //Minimum distance to the cylinders wall.
    double d_ = (distance_to_cylinder - radius);
   // return d_>0.0?d_:0.0;
    return d_;
}


