//!  Cylinder Obstacle Derived Class =============================================================/
/*!
*   \details   Cylinder class derived from an Obstacle. Defines infinite long cylinders
*              in the direction set by P,Q.
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
=================================================================================================*/


#ifndef CYLINDER_H
#define CYLINDER_H

#include "obstacle.h"


class Cylinder : public Obstacle
{
public:

    static int count;

    Eigen::Vector3d P,Q;    /*!< Cilinder Axis reference Points, P should be the "center"       */
    Eigen::Vector3d D;      /*!< Pre-computed and normalized P - Q vector                       */
    double radius;          /*!< Radius of the cylinder                                         */

    /*!
     *  \brief Default constructor. Does nothing
     */
    Cylinder();

    ~Cylinder();

    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    Cylinder(Eigen::Vector3d P_, Eigen::Vector3d Q_, double radius_, double scale = 1):P(P_*scale),Q(Q_*scale),radius(radius_*scale){
        D  = (Q_-P_).normalized();
        Q = P+D;
        id = count++;
    }

    /*!
     *  \param P_ Cylinder origin
     *  \param Q_ cylinder direction.
     *  \param radius_ cylinder's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    Cylinder(Cylinder const &cyl);

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
     *  \brief Returns the minimum distance from the walker to the cylinder. Used to set the reachable
     *  cylinders that a given walker can reach.
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

#endif // CYLINDER_H
