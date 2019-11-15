#include "voxel.h"
#include "constants.h"
#include <Eigen/Dense>

Voxel::Voxel()
{

}

Voxel::Voxel(Eigen::Vector3d min_limits_, Eigen::Vector3d max_limits_):min_limits(min_limits_),max_limits(max_limits_)
{
    //Cube vertex
    Eigen::Vector3d a = min_limits;
    Eigen::Vector3d b = a; b[0] = max_limits[0];
    Eigen::Vector3d c = b; c[1] = max_limits[1];
    Eigen::Vector3d d = c; d[0] = min_limits[0];

    Eigen::Vector3d e = a; e[2] = max_limits[2];
    Eigen::Vector3d f = b; f[2] = max_limits[2];
    Eigen::Vector3d g = c; g[2] = max_limits[2];
    Eigen::Vector3d h = d; h[2] = max_limits[2];

    walls[0] = Plane(a,b,c,d);
    walls[1] = Plane(b,f,g,c);
    walls[2] = Plane(e,f,g,h);
    walls[3] = Plane(a,e,h,d);
    walls[4] = Plane(a,b,f,e);
    walls[5] = Plane(d,c,g,h);
}

bool Voxel::CheckCollision(Walker &walker, Eigen::Vector3d& step, double& tmax, Collision &colision)
{
    Collision colision_temp = colision;

    for(int i = 0 ; i < 6; i++){
        walls[i].CheckCollision(walker,step,tmax,colision_temp);

        if(colision_temp.doIHaveMorePiorityThan(colision))
            colision = colision_temp;
    }

    return colision.type != Collision::null;
}

Plane::Plane(Eigen::Vector3d normal_, Eigen::Vector3d plane_center_, double d_):normal(normal_),plane_center(plane_center_),d(d_)
{

}

Plane::Plane(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& d_)
{
    normal = (b-a).cross(c-a);
    normal.normalize();
    plane_center = (a+b+c+d_)/4.0;
    d = normal.dot(a);
}

bool Plane::CheckCollision(Walker &walker, Eigen::Vector3d& step, double tmax, Collision& colision)
{
    colision.type = Collision::null;

    double t = (d - normal.dot(walker.pos_v))/(normal.dot(step));

    if(walker.status == Walker::on_voxel){
        if(t < EPS_VAL || t  > tmax)
            return false;
    }
    else{
        if(t < 0.0 || t >  tmax)
            return false;
    }

    colision.type = Collision::hit;
    colision.t = t;
    colision.bounced_direction = step;
    colision.col_location = Collision::voxel;
    colision.colision_point = walker.pos_v + t*step;

    return true;
}
