#include "statdistribution.h"
#include <algorithm>    // std::sort
#include <random>

using namespace std;
using namespace Eigen;


StatDistribution::StatDistribution(){}



void StatDistribution::displayDistribution(){}


std::vector<double> StatDistribution::createRadiiList(){
    return  std::vector<double> (1, 0.0);
}
