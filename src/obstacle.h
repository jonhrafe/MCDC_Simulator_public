//!  Obstacle Base Class ==============================================================================/
/*!
*   \details   Father class to define the base of any other obstacle (wall or substrate)
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
 =====================================================================================================*/

#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "collision.h"
#include "walker.h"
#include <Eigen/Core>
#include <random>


class Obstacle
{
public:
    int id;                         /*!< Unique id of the simulation                                                */

    int count_perc_crossings;       /*!< Auxiliar value to count the number of percolatin crossings in a simulation */
    double percolation;             /*!< Percolation value between 0 and 1.                                         */
    double T2;                      /*!< T2 decay, not used by default                                              */
    double prob_cross_e_i;
    double prob_cross_i_e;
    double diffusivity_i;
    double diffusivity_e;
    

    /*! \fn  Obstacle
     *  \brief Default constructor. Does nothing.
     */
    Obstacle();

    /*! \fn  checkCollision
     *  \param walker, Walker instance in the simulation.
     *  \param 3d step. Is assumed to be normalized.
     *  \param step_lenght, length used as the maximum step collision distance.
     *  \param colilsion, Collision instance to save the collision (if any) details.
     *  \return true only if there was a Collision::hit status. \see Collision.
     *  \brief Basic collision function. Returns the if there was any collision on against the obstacle.
     */
    bool checkCollision(Walker& walker, Eigen::Array3d& step,const double& step_lenght, Collision& colision);

    /*! \fn     elasticBounceAgainsPlane
     */
    void elasticBounceAgainsPlane(Eigen::Vector3d& ray_origin, Eigen::Vector3d& normal, double& t, Eigen::Vector3d &step);

    /*!
     *  \param  walker to find the (closest) distance.
     *  \brief  Returns the minimum distance of collision.
     */
    double minDistance(Walker& w);

    void setPercolation(double& percolation_);

    void setDiffusion(double& diffusivity_i_, double& diffusivity_e_);


};

#endif // OBSTACLE_H
