#include "obstacledistribution.h"
#include <algorithm>    // std::sort
#include <random>

using namespace std;
using namespace Eigen;

ObstacleDistribution::ObstacleDistribution(){}

ObstacleDistribution::ObstacleDistribution(double &icvf_, Eigen::Vector3d &min_l, Eigen::Vector3d &max_l, std::vector<double> &radiis_)
{
    icvf = icvf_;
    min_limits = min_l;
    max_limits = max_l;
    radiis = radiis_;
}

void ObstacleDistribution::createSubstrate(){}

void ObstacleDistribution::printSubstrate(ostream &out){}

bool ObstacleDistribution::checkForCollition(Obstacle obs, Vector3d min_limits, Vector3d max_limits, std::vector<Obstacle>& obstacles_to_add,double &min_distance)
{
    return false;
}

void ObstacleDistribution::checkBoundaryConditions(Obstacle obs, std::vector<Obstacle>& obstacles_to_add, Vector3d min_limits, Vector3d max_limits){}

double ObstacleDistribution::computeICVF(std::vector<Obstacle>& obstacles, Vector3d& min_limits, Vector3d& max_limits,int& num_no_repeat)
{
    return 0.0;
}

void ObstacleDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_,Eigen::Vector3d& l){}
