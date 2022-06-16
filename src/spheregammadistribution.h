//!  SphereGammaDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from a Gamma distribution of radiis placed in
*              a single voxel structure.
*   \author    JR and Remy G
*   \date      January 2021
*   \version   1.5
=================================================================================================*/


#ifndef SPHEREGAMMADISTRIBUTION_H
#define SPHEREGAMMADISTRIBUTION_H

#include "cylindergammadistribution.h"
#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "sphere.h"


class SphereGammaDistribution : public CylinderGammaDistribution
{
public:

    std::vector<Sphere> spheres;                    /*!< Spheres vector                                                             */

    /*!
     *  \param P_ Sphere origin
     *  \param radius_ sphere's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    SphereGammaDistribution(unsigned ,double,double, double,Eigen::Vector3d &, Eigen::Vector3d &,double min_radii = 0.001);

    /*!
     *  \brief Samples and constructs a Gamma distribution
    */
    void createGammaSubstrate();

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

    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);


};

#endif // SPHEREGAMMADISTRIBUTION_H
