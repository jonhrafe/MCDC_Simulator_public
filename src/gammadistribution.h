//!  GammaDistribution Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from Multi Gamma distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      March 2022
*   \version   2.0
=================================================================================================*/

#ifndef GAMMADISTRIBUTION_H
#define GAMMADISTRIBUTION_H


#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "statdistribution.h"


class GammaDistribution : public StatDistribution
{
public:

    std::vector<double> alpha;                                   /*!< alpha coefficient of the Gamma distribution                                */
    std::vector<double> beta;                                    /*!< beta coefficient of the gamma distribution                                 */
    std::vector<unsigned> num_obstacles;              /*!< number of spheres fit inside the substrate                                */
    double radius_min;
    double radius_max;


    /*!
     *  \brief Default constructor. Does nothing
     */
    GammaDistribution(); 

    GammaDistribution(std::vector<unsigned> ,std::vector<double>, std::vector<double>, double r_min=EPS_VAL, double r_max=10e6);

    /*!
     *  \brief Shows a small histogram of the gamma distribution
    */
    void displayDistribution();

    std::vector<double> createRadiiList();

};

#endif // GAMMADISTRIBUTION_H