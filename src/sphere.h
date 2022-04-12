//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere class derived from an Obstacle. Defines sphere of radius R
*   \author    Remy Gardier
*   \date      January 2021
*   \version   0.0
=================================================================================================*/


#ifndef SPHERE_H
#define SPHERE_H

#include "obstacle.h"


class Sphere : public Obstacle
{
public:

    static int count;

    Eigen::Vector3d P;      /*!< Center of the sphere   */
    double radius;          /*!< Radius of the sphere   */

    double volume;

    /*!
     *  \brief Default constructor. Does nothing
     */
    Sphere();

    ~Sphere(); 

    /*!
     *  \param P_ Sphere origin
     *  \param radius_ sphere's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    Sphere(Eigen::Vector3d P_, double radius_, double scale = 1, double percolation_=0.0):P(P_*scale), radius(radius_*scale){
        percolation = percolation_;
        id = count++;
        volume = 4./3.*M_PI * (radius_*scale) *  (radius_*scale)  *  (radius_*scale);
    }

    /*!
     *  \param P_ Sphere origin
     *  \param radius_ sphere's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
   Sphere(Sphere const &sph);

    /*! \fn  checkCollision
     *  \param walker, Walker instance in the simulation.
     *  \param 3d step. Is assumed to be normalized.
     *  \param step_length, length used as the maximum step collision distance.
     *  \param collision, Collision instance to save the collision (if any) details.
     *  \return true only if there was a Collision::hit status. \see Collision.
     *  \brief Basic collision function. Returns the if there was any collision on against the obstacle.
     */
    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    /*! \fn  minDistance
     *  \param walker, Walker instance in the simulation.
     *  \brief Returns the minimum distance from the walker to the sphere. Used to set the reachable
     *  sphere that a given walker can reach.
     */
    double minDistance(Walker &w);

private:

    /*! \fn  handleCollition
     *  \param walker, Walker instance in the simulation.
     *  \param collision, Collision instance to save all the information.
     *  \param step, step vector where to move.
     *  \brief Returns true if it was any analytical collision to the infinite plane
     */
    inline bool handleCollition(Walker& walker, Collision &colision, Eigen::Vector3d& step,double& a,double& b, double& c,double& discr,double& step_length);

};

#endif // SPHERE_H
