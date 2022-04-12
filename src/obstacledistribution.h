//!  ObstacleDistribution Class =============================================================/
/*!
*   \details   Base class to construct a substrate made of obstacles from a distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      April 2022
*   \version   0.0
=================================================================================================*/

#ifndef OBSTACLEDISTRIBUTION_H
#define OBSTACLEDISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "obstacle.h" 
#include "statdistribution.h"

class ObstacleDistribution
{
public:

    double icvf;                                    /*!< Achieved intra-celular volum fraction in the substrate                     */

    Eigen::Vector3d min_limits;                     /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits;                     /*!< voxel max limits (if any)                                                  */
    std::vector<double> radiis;                     /*!< List of radii of the obstacles                                             */
    
    
    /*!
     *  \brief Default constructor. Does nothing
     */
    ObstacleDistribution(); 

    /*!
     *  \param icvf_ Intra-cellular volume fraction
     *  \param min_l voxel min limits
     *  \param max_l voxel min limits
     *  \brief Initialize everything.
     */
    ObstacleDistribution(double& ,Eigen::Vector3d &, Eigen::Vector3d &, std::vector<double> &);


    /*!
     *  \brief Samples and constructs a distribution
    */
    void createSubstrate();

    /*!
     *  \brief Prints the obstacle positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out); 

protected:

    /*!
     *  \brief Checks for collision between inside a voxel (with periodic boundaries)
     *  \param obs New obstacle
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param obstacles_list obstacle already added.
     *  \param min_distance that two spheres can be close to.
    */
    bool checkForCollition(Obstacle obs, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector<Obstacle>& obstacles_list, double &min_distance);

    /*!
     *  \brief Auxiliary function to check the BOundary collision
     *  \param obs  obstacle to check collision with.
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param obstacles_to_add spheres already added.
    */
    void checkBoundaryConditions(Obstacle obs, std::vector<Obstacle>& obstacles_to_add, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits);

    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added spheres.
     *  \param obstacles List of included spheres.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
    */
    double  computeICVF(std::vector<Obstacle> &obstacles, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);

    /*!
     *  \brief Estimate the minimal voxel side length from Intra Celular Volum Fraction and the radiis of the obstacles
     *  \param radiis List of radiis.
     *  \param icvf_ Intra-cellular volume fraction
     *  \param l voxel side length.
    */
    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);

};

#endif // OBSTACLEDISTRIBUTION_H
