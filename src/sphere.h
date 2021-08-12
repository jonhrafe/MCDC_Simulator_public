//!  Sphere Obstacle Derived Class =============================================================/
/*!
*   \details   Sphere obstacle class derived from an Obstacle. Defines a analyitical sphere of radius
*   r and centered in center
*   \author    Jonathan Rafael
*   \date      2020
*   \version   1.5
=================================================================================================*/

#ifndef SPHERE_H
#define SPHERE_H

#include "obstacle.h"

class Sphere : public Obstacle
{
public:

    static int count;

    Eigen::Vector3d center;    /*!< Cilinder Axis reference Points, P should be the "center"      */
    double radius;             /*!< Radius of the Sphere                                          */

    /*!
     *  \brief Default constructor. Does nothing
     */
    Sphere(){id = count++;}
    /*!
     *  \brief Default destructor. Does nothing
     */
    ~Sphere(){count--;}

    /*!
     *  \param center Sphere origin
     *  \param radius Sphere's radius
     *  \param scale  overall scale for when reading files.
     *  \brief Initialize everything.
     */
    Sphere(Eigen::Vector3d center_, double radius_,double scale =1):center(center_*scale),radius(radius_*scale){
        id = count++;
    }

    /*!
     *  \brief constrcutor by copy
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
     *  \brief Returns the minimum distance from the walker to the Sphere. Used to set the reachable
     *  Spheres that a given walker can reach.
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

#endif // Sphere_H
