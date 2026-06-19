//!  Obstacle Base Class ==============================================================================/
/*!
*   \details   Father class to define the base of any other obstacle (wall or substrate)
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.6
 =====================================================================================================*/

#ifndef OBSTACLE_H
#define OBSTACLE_H
#include "collision.h"
#include "walker.h"
#include "Eigen/Core"
class Obstacle
{
public:

    int id;                         /*!< Unique id of the simulation                                                */
    int count_perc_crossings;       /*!< Auxiliar value to count the number of percolatin crossings in a simulation */
    // Directional permeability counters (validation of the Powles model). A "hit"
    // is counted whenever the membrane is reached and a crossing draw is made
    // (covers both step and bouncing collisions); a "cross" when that draw succeeds.
    // Empirical p_hat = count_cross_* / count_hits_* should match prob_cross_*.
    // NOTE: obstacles are shared across processes, so these are exact only for
    // num_process 1 (like count_perc_crossings, they race under multithreading).
    unsigned long count_hits_i_e;   /*!< membrane reached from the intra side (draw made)                            */
    unsigned long count_hits_e_i;   /*!< membrane reached from the extra side (draw made)                            */
    unsigned long count_cross_i_e;  /*!< successful intra->extra crossings                                          */
    unsigned long count_cross_e_i;  /*!< successful extra->intra crossings                                          */
    double percolation;             /*!< Percolation value between 0 and 1.                                         */
    double permeability;            /*!< Membrane permeability kappa (velocity; m/s in the .conf, mm/ms internal). Physical input; per-encounter statistics derived from it. */
    double T2;                      /*!< T2 decay, not used by default                                              */
    double d_intra;                 /*!< Internal Diffusion coefficient                                             */
    double prob_cross_e_i;         /*!< Probability of crossing from the the exterior                               */
    double prob_cross_i_e;         /*!< Probability of crossing to   the the interior                               */

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

};

#endif // OBSTACLE_H
