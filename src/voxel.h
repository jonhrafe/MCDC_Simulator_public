//!  Main class. Implements basic voxel limits and operations. =================================================/
/*! Class to handle and manage the voxels in the simulations. So far only one voxel at the time can
 *  be handled. To improve to several voxels, modifications shall be done.
 *
*   \author Jonathan Rafael
*   \date   July 2016
*   \version   0.2
*=============================================================================================================+*/

#ifndef VOXEL_H
#define VOXEL_H
#include <Eigen/Core>
#include <collision.h>
#include "walker.h"

/*! \class Plane
 *  \brief Auxiliary class to implements plane's interactions with particles.
 */
class Plane{
public:
    //Enough to save the parametric representation
    Eigen::Vector3d normal, plane_center;
    double d;

    //X-Y plane in the origin
    Plane(){
        normal = {0,0,1};
        plane_center = {0,0,0};
        d = 0;
    }

    Plane(Eigen::Vector3d normal_, Eigen::Vector3d plane_center_, double d_);

    Plane(Eigen::Vector3d &a, Eigen::Vector3d &b, Eigen::Vector3d &c, Eigen::Vector3d &d);

    bool CheckCollision(Walker& walker, Eigen::Vector3d &step, double tmax, Collision &colision);
};

/*! \class Voxel
 *  \brief //!  Main class. Implements basic voxel limits and operations. Class to handle and manage
 *  the voxels in the simulations. So far only one voxel at the time can be handled. To improve to
 *  several voxels, modifications shall be done.
 */
class Voxel
{
public:
    Eigen::Vector3d min_limits, max_limits;

    Voxel();

    Voxel(Eigen::Vector3d min_limits_,Eigen::Vector3d max_limits_);

    Plane walls[6];

    bool CheckCollision(Walker& walker, Eigen::Vector3d &step, double &tmax, Collision& colision);
};



#endif // VOXEL_H
