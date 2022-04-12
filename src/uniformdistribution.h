//!  UniformDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from Multi uniform distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      March 2022
*   \version   2.0
=================================================================================================*/

#ifndef UNIFORMDISTRIBUTION_H
#define UNIFORMDISTRIBUTION_H


#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "statdistribution.h"


class UniformDistribution : public StatDistribution
{
public:

    std::vector<double> radiis_sphere;                /*!< Radii of each population of spheres                                       */
    std::vector<unsigned> num_obstacles;              /*!< number of spheres fit inside the substrate                                */
    double radius_min;
    double radius_max;


    /*!
     *  \brief Default constructor. Does nothing
     */
    UniformDistribution(); 

    UniformDistribution(std::vector<unsigned> ,std::vector<double>, double r_min=EPS_VAL, double r_max=10e6);

    /*!
     *  \brief Shows a small histogram of the gamma distribution
    */
    void displayDistribution();

    std::vector<double> createRadiiList();

};

#endif // UNIFORMDISTRIBUTION_H