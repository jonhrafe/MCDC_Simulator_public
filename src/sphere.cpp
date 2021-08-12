#include "sphere.h"
#include "constants.h"
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;

int Sphere::count = 0;

Sphere::Sphere(const Sphere &sph)
{
    center = sph.center;
    radius = sph.radius;
    id = count++;
}

bool Sphere::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{

    //Origin of the ray
    Vector3d O;
    walker.getVoxelPosition(O);
    Vector3d m = O - this->center;

    // total distance
    double distance_to_sphere = m.norm();
    // collision distance
    double d_ = distance_to_sphere - radius;

    //If the minimum distance from the walker to the cylinder is more than
    // the actual step size, we can discard this collision.
    if(d_> EPS_VAL){
        if(d_ > step_lenght+barrier_tickness){
            return false;
        }
    }

    double a = 1;
    double b = m.dot(step);
    double c = m.dot(m) - radius*radius;
    if(b > EPS_VAL && c > EPS_VAL)
        return false;


    double discr = b*b - a*c;
    // no real roots
    if(discr <= 0.0){
        colision.type = Collision::null;
        return false;
    }

    //if we arrived here we need to compute the quadratic equation.
    return handleCollition(walker,colision,step,a,b,c,discr,step_lenght);

}


inline bool Sphere::handleCollition(Walker& walker, Collision &colision, Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length){

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
        double _percolation_ (double(rand())/RAND_MAX);

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
    colision.obstacle_ind = -1;

    if(c<-1e-10){
        colision.col_location = Collision::inside;
        walker.in_obj_index = -1;
    }
    else if(c>1e-10){
        colision.col_location = Collision::outside;
    }
    else{
        colision.col_location = Collision::unknown;
    }

    colision.rn = c;
    colision.colision_point = walker.pos_v + colision.t*step;

    //Normal point
    Eigen::Vector3d normal = (colision.colision_point-this->center).normalized();
    Eigen::Vector3d temp_step = step;
    elasticBounceAgainsPlane(walker.pos_v,normal,colision.t,temp_step);

    colision.bounced_direction = temp_step.normalized();

    return true;

}

double Sphere::minDistance(Walker &w){

    //Origin of the ray
    Vector3d O;
    w.getVoxelPosition(O);
    Vector3d m = O - this->center;
    // minimum distance to the cylinder axis.
    double distance_to_cylinder = m.norm();

    //Minimum distance to the cylinders wall.
    double d_ = (distance_to_cylinder - radius);
   // return d_>0.0?d_:0.0;
    return d_;
}
