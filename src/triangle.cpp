#include "triangle.h"
#include "collision.h"
#include <limits>       // std::numeric_limits
#include <math.h>
#include <cstddef>
#include "Eigen/Dense"
#include <iostream>
#include "constants.h"

using namespace std;

Triangle::Triangle()
{

    indexes[0]=0;
    indexes[1]=1;
    indexes[2]=2;
    vertices = NULL;
    normal[0]=0;
    normal[1]=0;
    normal[2]=0;
}

void Triangle::getVertex(const unsigned i, Eigen::Vector3d &v) {
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

#if PRECISE_T_MIN_D == 1
double Triangle::minDistance(const Eigen::Vector3d p) {
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

static const double INTERSECTION_EPS = 1e-12;

/**
 * @brief Check if two 3D points are close to each other (same vertex).
 */
inline bool isClose(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, double eps = 1e-12) {
    return (p1 - p2).squaredNorm() < eps*eps;
}

/**
 * @brief Check if a point P is inside the triangle formed by A, B, C (using barycentric coordinates).
 * @param P - the point to test
 * @param A,B,C - the triangle vertices
 * @return true if inside (or on edge), false otherwise
 */
bool pointInTriangle(const Eigen::Vector3d &P,
                     const Eigen::Vector3d &A,
                     const Eigen::Vector3d &B,
                     const Eigen::Vector3d &C)
{
    // Compute vectors
    Eigen::Vector3d v0 = C - A;
    Eigen::Vector3d v1 = B - A;
    Eigen::Vector3d v2 = P - A;

    // Compute dot products
    double dot00 = v0.dot(v0);
    double dot01 = v0.dot(v1);
    double dot02 = v0.dot(v2);
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);

    // Compute barycentric coordinates
    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= -INTERSECTION_EPS) && (v >= -INTERSECTION_EPS) &&
           (u + v <= 1.0 + INTERSECTION_EPS);
}

/**
 * @brief Check intersection of a line segment [P0,P1] with triangle [A,B,C].
 *        This is a standard segment-triangle intersection using a Möller–Trumbore style test.
 * @param P0,P1 the endpoints of the segment
 * @param A,B,C the triangle’s vertices
 * @param tOut  the intersection parameter (0 <= t <= 1 for segment intersection)
 * @return true if they intersect, false otherwise
 */
bool segmentIntersectsTriangle(const Eigen::Vector3d &P0,
                               const Eigen::Vector3d &P1,
                               const Eigen::Vector3d &A,
                               const Eigen::Vector3d &B,
                               const Eigen::Vector3d &C,
                               double &tOut)
{
    const double EPS = 1e-10;

    Eigen::Vector3d e1 = B - A;
    Eigen::Vector3d e2 = C - A;
    Eigen::Vector3d d  = P1 - P0;  // segment direction

    Eigen::Vector3d p  = d.cross(e2);
    double det = e1.dot(p);

    // If det is near zero, there is no intersection (or the line is parallel).
    if (fabs(det) < EPS) {
        return false;
    }

    double invDet = 1.0 / det;
    Eigen::Vector3d T = P0 - A;

    // Calculate U parameter
    double u = T.dot(p) * invDet;
    if (u < 0.0 - EPS || u > 1.0 + EPS) {
        return false;
    }

    // Calculate V parameter
    Eigen::Vector3d q = T.cross(e1);
    double v = d.dot(q) * invDet;
    if (v < 0.0 - EPS || (u + v) > 1.0 + EPS) {
        return false;
    }

    // Calculate t to find out where intersection happens on the line
    double t = e2.dot(q) * invDet;
    if (t < 0.0 - EPS || t > 1.0 + EPS) {
        // For pure line intersection, we wouldn't check 0..1, 
        // but for segment we need t in [0,1].
        return false;
    }

    tOut = t;
    return true;
}

/**
 * @brief Check if any edge of triangle1 intersects triangle2
 */
bool anyEdgeIntersectsTriangle(const Eigen::Vector3d &A1,
                               const Eigen::Vector3d &B1,
                               const Eigen::Vector3d &C1,
                               const Eigen::Vector3d &A2,
                               const Eigen::Vector3d &B2,
                               const Eigen::Vector3d &C2)
{
    double tIgnore;
    // Check the 3 edges of the first triangle
    if (segmentIntersectsTriangle(A1, B1, A2, B2, C2, tIgnore)) return true;
    if (segmentIntersectsTriangle(B1, C1, A2, B2, C2, tIgnore)) return true;
    if (segmentIntersectsTriangle(C1, A1, A2, B2, C2, tIgnore)) return true;
    return false;
}

/**
 * @brief The main function to check if two triangles intersect.
 *        - Skips intersection if the triangles appear to share any vertex (adjacent).
 *        - Otherwise performs edge-edge tests and "vertex-in-triangle" tests.
 *
 * @param other the other triangle
 * @return true if they intersect, false otherwise
 */
bool Triangle::triangleIntersects(Triangle &other)
{
    // 1) Check if they share a vertex (adjacency check).
    //    If yes, skip intersection test and return false (or do whatever is needed).
    Eigen::Vector3d A1, B1, C1;
    getVertex(0, A1);
    getVertex(1, B1);
    getVertex(2, C1);

    Eigen::Vector3d A2, B2, C2;
    other.getVertex(0, A2);
    other.getVertex(1, B2);
    other.getVertex(2, C2);



    // Epsilon for "almost the same" vertex
    const double adjacencyEps = 1e-12;
    // Compare all 3 vertices from "this" to all 3 from "other"
    Eigen::Vector3d tri1Vertices[3] = {A1, B1, C1};
    Eigen::Vector3d tri2Vertices[3] = {A2, B2, C2};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (isClose(tri1Vertices[i], tri2Vertices[j], adjacencyEps)) {
                // They share (approximately) the same vertex => adjacent (same mesh?), skip
                return false;
            }
        }
    }
    

    // 2) Check edges of one triangle vs. the other
    if (anyEdgeIntersectsTriangle(A1, B1, C1, A2, B2, C2)) return true;
    if (anyEdgeIntersectsTriangle(A2, B2, C2, A1, B1, C1)) return true;

    // If none of the above conditions is true, the triangles do not intersect.
    return false;
}