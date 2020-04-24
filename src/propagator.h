#ifndef PROPAGATOR_H
#define PROPAGATOR_H

#include <Eigen/Core>
#include <vector>

#include "win_types.h"  // quick hack for compilitation on Windows

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
