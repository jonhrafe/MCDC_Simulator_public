#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include "constants.h"
#include "Eigen/Core"
#include <vector>


class Propagator
{
public:

    uint num_dirs =0;
    uint num_times = 0;
    Eigen::Matrix3Xf directions;

    std::vector<unsigned> log_times;

    std::vector<std::vector<float>> propagator_log;

    Propagator();

    void initPropagator();

};

#endif // PROPAGATOR_H
