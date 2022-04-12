//!  GaussianDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from Multi Gaussian distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      March 2022
*   \version   2.0
=================================================================================================*/

#ifndef GAUSSIANDISTRIBUTION_H
#define GAUSSIANDISTRIBUTION_H


#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "statdistribution.h"


class GaussianDistribution : public StatDistribution
{
public:

    std::vector<double> mu_mean;                       /*!< Mean of the Gaussian distribution                                */
    std::vector<double> sigma_std;                    /*!< STD of the Gaussian distribution                                 */
    std::vector<unsigned> num_obstacles;            /*!< number of spheres fit inside the substrate                       */
    double radius_min;
    double radius_max;


    /*!
     *  \brief Default constructor. Does nothing
     */
    GaussianDistribution(); 

    GaussianDistribution(std::vector<unsigned> ,std::vector<double>, std::vector<double>, double r_min=EPS_VAL, double r_max=10e6);

    /*!
     *  \brief Shows a small histogram of the gamma distribution
    */
    void displayDistribution();

    std::vector<double> createRadiiList();

};

#endif // GAUSSIANDISTRIBUTION_H