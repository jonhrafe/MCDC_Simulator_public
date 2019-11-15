#include "propagator.h"

Propagator::Propagator()
{


}

void Propagator::initPropagator()
{
    for (uint i = 0 ; i < num_times; i++){
        std::vector<float> jkr(this->num_dirs,0.0);
        propagator_log.push_back(jkr);
    }
}
