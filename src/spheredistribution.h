//!  SphereDistribution Class =============================================================/
/*!
*   \details   Base class to construct a substrate made of spheres from a distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      April 2022
*   \version   0.0
=================================================================================================*/

#ifndef SPHEREDISTRIBUTION_H
#define SPHEREDISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "sphere.h"
#include "obstacledistribution.h"

class SphereDistribution : public ObstacleDistribution
{
public:

    std::vector<Sphere> spheres;                    /*!< Spheres vector                                                             */
    

    /*!
     *  \brief Default constructor. Does nothing
     */
    SphereDistribution(); 


    /*!
     *  \param icvf_ Intra-cellular volume fraction
     *  \param min_l voxel min limits
     *  \param max_l voxel min limits
     *  \brief Initialize everything.
     */
    SphereDistribution(double &, Eigen::Vector3d &, Eigen::Vector3d &, std::vector<double> &);

    /*!
     *  \brief Samples and constructs a distribution
    */
    void createSubstrate();

    /*!
     *  \brief Prints the sphere positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out); 


private:

    /*!
     *  \brief Checks for collision between inside a voxel (with periodic boundaries)
     *  \param sph sphere to check collision with
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param spheres_list spheres already added.
     *  \param min_distance that two spheres can be close to.
    */
    bool checkForCollition(Sphere sph, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector<Sphere>& spheres_list, double &min_distance);

    /*!
     *  \brief Auxiliary function to check the BOundary collision
     *  \param sph sphere to check collision with.
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param spheres_list spheres already added.
    */
    void checkBoundaryConditions(Sphere sph, std::vector<Sphere>& spheres_list, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits);

    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added spheres.
     *  \param spheres List of included spheres.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
    */
    double  computeICVF(std::vector<Sphere> &spheres, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);

    /*!
     *  \brief Estimate the minimal voxel side length from Intra Celular Volum Fraction and the radiis of the obstacles
     *  \param radiis List of radiis.
     *  \param icvf_ Intra-cellular volume fraction
     *  \param l voxel side length.
    */
    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);

};

#endif // SPHEREDISTRIBUTION_H
