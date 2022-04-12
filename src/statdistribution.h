#ifndef STATSDISTRIBUTION_H
#define STATSDISTRIBUTION_H


#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>

class StatDistribution
{
public:

    

    StatDistribution(); 

    /*!
     *  \brief Shows a small histogram of the distribution
    */
    void displayDistribution();

    std::vector<double> createRadiiList();

};

#endif // STATSDISTRIBUTION_H