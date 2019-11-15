#include "triangle.h"
#include "collision.h"
#include <limits>       // std::numeric_limits
#include <math.h>
#include <cstddef>
#include <Eigen/Dense>
#include <iostream>
#include "constants.h"

using namespace std;

Triangle::Triangle()
{
    vertices = NULL;
    normal[0]=0;
    normal[1]=0;
    normal[2]=0;
}

void Triangle::getVertex(const unsigned i, Eigen::Vector3d &v){
    v[0] = vertices[indexes[i]].points[0];
    v[1] = vertices[indexes[i]].points[1];
    v[2] = vertices[indexes[i]].points[2];
}

void Triangle::getNormal(Eigen::Vector3d &normal)
{
    if( (this->normal[0] == 0.0) && (this->normal[1] == 0.0) && (this->normal[2] == 0.0))
      saveNormalAndAuxInfo();

    normal= this->normal;
}

void Triangle::saveNormalAndAuxInfo()
{
    Eigen::Vector3d a,b,c,A,B;
    getVertex(0,a);
    getVertex(1,b);
    getVertex(2,c);
    A = a-b;
    B = a-c;

    this->normal = (A.cross(B)).normalized();

    this->center = (a+b+c)/3.0;

    double d1 = (center-a).squaredNorm();
    double d2 = (center-b).squaredNorm();
    double d3 = (center-c).squaredNorm();
    this->radius =sqrt(fmax(d1,fmax(d2,d3)));
}

bool Triangle::rayIntersects(const Eigen::Vector3d &ray_origin, const Eigen::Vector3d &step, double &t)
{
    Eigen::Vector3d e1,e2,pvec,tvec,qvec;
    Eigen::Vector3d a,b,c;
    getVertex(0,a);
    getVertex(1,b);
    getVertex(2,c);
    e1 = b-a;
    e2 = c-a;
    pvec = step.cross(e2);

    double det = e1.dot(pvec);

    if (det > -triangle_eps && det < triangle_eps )
        return false;

    double invDet = 1.0/det;

    tvec = ray_origin - a;
    double u = tvec.dot(pvec)* invDet;

    if (u < -triangle_eps || u  > 1.0+triangle_eps)
        return false;

    qvec = tvec.cross(e1);

    double v = step.dot(qvec) * invDet;

    if (v < -triangle_eps  ||  u + v  > 1.0+triangle_eps)
        return false;

    t = e2.dot(qvec) * invDet;


    return true;
}

bool Triangle::rayIntersects_MT(const Eigen::Vector3d & ray_origin, const Eigen::Vector3d &step, double &u, double &v, double &t)
{
    double EPS = 1e-13;

    Eigen::Vector3d e1,e2,pvec,tvec,qvec;
    Eigen::Vector3d a,b,c;
    getVertex(0,a);
    getVertex(1,b);
    getVertex(2,c);

    e1 = b-a;
    e2 = c-a;

    pvec = step.cross(e2);

    double det = e1.dot(pvec);

    //if determinant is near zero, ray lies in plane of triangle or ray is parallel to plane of triangle
    if (det > -EPS && det < EPS )
        return false;

    double invDet = 1.0/det;

    tvec = ray_origin-a;
    u = tvec.dot(pvec) * invDet;

    if (u < -triangle_eps || u > 1.0+triangle_eps)
        return false;

    qvec = tvec.cross(e1);
    v = step.dot(qvec) * invDet;

    if (v < -triangle_eps  ||  u + v > 1+triangle_eps)
        return false;

    t = e2.dot(qvec) * invDet;

    return true;
}

#if PRECISE_T_MIN_D == 0
double Triangle::minDistance(const Eigen::Vector3d p){
        //    distance to sphere
        return fmax(0,(p-center).norm()-radius);
}

#else

double Triangle::minDistance(const Eigen::Vector3d p)
{
    double EPS = 1e-1;

    Eigen::Vector3d a,b,c;
    getVertex(0,a);
    getVertex(1,b);
    getVertex(2,c);
    Eigen::Vector3d ab = b - a ;
    Eigen::Vector3d ac = c - a ;
    Eigen::Vector3d ap = b - c;

    // Check if P in vertex region outside A
    double d1 = ab.dot(ap);
    double d2 = ac.dot(ap);

    if(d1 <= -EPS && d2 <= -EPS)
        return (a-p).norm(); // barycentric coordinates (1,0,0)

    // Check if P in vertex region outside B
    Eigen::Vector3d bp = p - b;
    double d3 = ab.dot(bp);
    double d4 = ac.dot(bp);

    if(d3 <= -EPS && d4 <= d3)
        return (b-p).norm();  // barycentric coordinates (0,1,0)


    double vc = d1*d4 - d3*d2;

    if (vc <= -EPS && d1 >= EPS && d3 <= -EPS) {
        double v = d1 / (d1 - d3);
        return ((a + v*ab) - p).norm();             // barycentric coordinates (1-v,v,0)
    }

    // Check if P in vertex region outside C
    Eigen::Vector3d cp = p - c;
    double d5 = ab.dot(cp);
    double d6 = ac.dot(cp);

    if (d6 >= EPS && d5 <= d6)
        return (c-p).norm(); // barycentric coordinates (0,0,1)


    // Check if P in edge region of AC, if so return projection of P onto AC
    double vb = d5*d2 - d1*d6;
    if (vb <= -EPS && d2 >= EPS && d6 <= -EPS) {
        double w = d2 / (d2 - d6);
        return ((a+ w*ac) - p).norm();              // barycentric coordinates (1-w,0,w)
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    double  va = d3*d6 - d5*d4;
    if (va <= -EPS && (d4 - d3) >= EPS && (d5 - d6) >= EPS){
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return ((b + w*(c-b)) - p).norm();          // barycentric coordinates (0,1-w,w)
    }

    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    double  denom = 1.0 / (va + vb + vc);
    double v = vb * denom;
    double w = vc * denom;

    return ((a + ab * v + ac * w)-p).norm(); // = u*a + v*b + w*c, u = va * denom = 1.0f - v - w

}

#endif

void Triangle::stepIntersects_MT(Walker& walker, const Eigen::Vector3d &step, const double &max_length, Collision &colision)
{
    double EPS = 1e-14;
    Eigen::Vector3d e1,e2,pvec,tvec,qvec;
    Eigen::Vector3d a,b,c;
    getVertex(0,a);
    getVertex(1,b);
    getVertex(2,c);

    //Computation of the determinant of the system.

    e1 = b-a;
    e2 = c-a;
    pvec = step.cross(e2);
    double det = e1.dot(pvec);
    double t = numeric_limits<double>::infinity();

    //if determinant is near zero, ray lies in the triangle's plane or the ray is parallel to the triangle's plane
    if (det > -EPS && det < EPS ){
        colision.type = Collision::near;
        return;
    }

    double invDet = 1.0/det;

    // Computation of t
    tvec = walker.pos_v - a;
    qvec = tvec.cross(e1);

    t = e2.dot(qvec) * invDet;

    //First case, step is in the other direction of the step.
    //This means that we are moving in the opposite direction
    // Or we are too far away from the triangles plane anyway.
    if( t < 0.0 || ( t > max_length+barrier_tickness)){
        colision.type = Collision::null;
        return;
    }

    // a spin that's bouncing ignores collision at 0 (is in a wall)
    if(walker.status == Walker::bouncing){

        //if the collision are too close.
        if( ( t < EPS_VAL) ) {
            colision.type = Collision::null;
            return;
        }

        colision.t = fmin(t,max_length);
    }
    else{
        //if we are not bouncing, all collisions counts.
        colision.t = fmin(t,max_length);
    }

    // Computation of u
    colision.u = tvec.dot(pvec) * invDet;

    if (colision.u < -EPS_VAL || colision.u > 1.0+EPS_VAL){
        colision.type = Collision::null;
        return;
    }

    // Computation of c and u+v
    colision.v = step.dot(qvec) * invDet;

    if (colision.v < -EPS_VAL || colision.u + colision.v > 1.0+EPS_VAL){
        colision.type = Collision::null;
        return;
    }

    colision.type = Collision::hit;
}

void Triangle::stepIntersects_MT_limits(const Eigen::Vector3d &ray_origin, const Eigen::Vector3d &step, const double &max_length, Collision &colision,
                                        const Eigen::Vector3d &limits_mod, double limit_x, double limit_y, double limit_z)
{
    double EPS = 1e-15;
    Eigen::Vector3d e1,e2,pvec,tvec,qvec;
    Eigen::Vector3d a,b,c;
    getVertex(0,a);
    getVertex(1,b);
    getVertex(2,c);


    a[0]+= int(limits_mod[0])*limit_x; b[0]+= int(limits_mod[0])*limit_x; c[0]+= int(limits_mod[0])*limit_x;
    a[1]+= int(limits_mod[1])*limit_y; b[1]+= int(limits_mod[1])*limit_y; c[1]+= int(limits_mod[1])*limit_y;
    a[2]+= int(limits_mod[2])*limit_z; b[2]+= int(limits_mod[2])*limit_z; c[2]+= int(limits_mod[2])*limit_z;

    //Computation of the determinat of the system.
    e1 = b-a;
    e2 = c-a;
    pvec = step.cross(e2);
    double det = e1.dot(pvec);
    colision.t = numeric_limits<double>::infinity();

    //First case, det = 0;
    if (det > -EPS && det < EPS ){
        colision.type = Collision::near;
        return;
    }

    double invDet = 1.0/det;

    // Computation of t
    tvec = ray_origin - a;
    qvec = tvec.cross(e1);

    colision.t = e2.dot(qvec) * invDet;

    //First case, step is in the other direction of the step.
    //This means that we are moving in the opposite direction
    if(colision.t < 0){
        colision.type = Collision::null;
        return;
    }

    // Is the triangle plane on the distance of the step
    if(colision.t >= EPS && colision.t - max_length > EPS){
        colision.type = Collision::near;
        return;
    }

    // Computation of u
    colision.u = tvec.dot(pvec) * invDet;
    if (colision.u < EPS*2 || colision.u-1 > -EPS*2){
        colision.type = Collision::near;
        return;
    }

    // Computation of c and u+v
    colision.v = step.dot(qvec) * invDet;
    if (colision.v < EPS*2 || colision.v - 1 > -EPS*2  || colision.u + colision.v -1 > EPS*2){
        colision.type = Collision::near;
        return;
    }

    colision.type = Collision::hit;

}
