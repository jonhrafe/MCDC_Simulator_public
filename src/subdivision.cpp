#include "subdivision.h"

Subdivision::Subdivision()
{
    density = 0 ;
    density_intra = density_extra = 0;
}

Subdivision::Subdivision(Eigen::Vector3f &min_, Eigen::Vector3f &max_)
{
    this->min_limits = min_;
    this->max_limits = max_;
    this->density = 0;
    this->density_intra=0;
}

bool Subdivision::isInside(Eigen::Vector3d& pos)
{
    bool flag =  (pos(0) >= min_limits(0)) && (pos(1) >= min_limits(1)) && (pos(2) >= min_limits(2)) ;
    flag     &=  (pos(0) <= max_limits(0)) && (pos(1) <= max_limits(1)) && (pos(2) <= max_limits(2)) ;

    return flag;
}
