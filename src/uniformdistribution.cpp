#include "uniformdistribution.h"
#include <algorithm>    // std::sort
#include <random>

using namespace std;
using namespace Eigen;


UniformDistribution::UniformDistribution(){}


UniformDistribution::UniformDistribution(std::vector<unsigned> num_obs, std::vector<double> radii_spheres_, double r_min, double r_max)
{
    num_obstacles = num_obs; 
    radius_min = r_min;
    radius_max =r_max;
    radiis_sphere = radii_spheres_;
}

void UniformDistribution::displayDistribution()
{
    int nrolls=0;  // number of experiments
    const int nstars=100;    // maximum number of stars to distribute
    string message;
    
    int p[11]={};

    for(int i=0;i<this->num_obstacles.size(); ++i){
        double curr_r   = this->radiis_sphere[i];
        if (curr_r<10){
            p[int(curr_r)] += this->num_obstacles[i];
        }
        else{
            p[10]+= this->num_obstacles[i];
        }
        nrolls+=this->num_obstacles[i];
    }

    for (int i=0; i<9; ++i) {
        message = std::to_string(i) + "-" + std::to_string(i+1) + ": " + std::string(p[i]*nstars/nrolls,'*');
        SimErrno::info(message,cout);
    }
    message = "9-10:" + std::string(p[9]*nstars/nrolls,'*') ;
    SimErrno::info(message,cout);
    message = ">10: " +  std::string(p[10]*nstars/nrolls,'*') + "\n" ;
    SimErrno::info(message,cout);
}


std::vector<double> UniformDistribution::createRadiiList(){
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);


    unsigned tot_num_obstacles = 0; 
    for(int i=0;i<this->num_obstacles.size(); ++i){
        tot_num_obstacles += this->num_obstacles[i];
    }
   
    std::vector<double> radiis_packing(tot_num_obstacles,0.0);

    unsigned counter    = 0; 
    double rad_curr     = 0.0;

    for(int i=0;i<this->num_obstacles.size(); ++i){
        unsigned curr_counter    = 0; 
        rad_curr = this->radiis_sphere[i]*1e-3 ;
                
        double curr_num_obstacle = this->num_obstacles[i];

        while(curr_counter < curr_num_obstacle){
            double random_rad = udist(gen) * 1e-10;

            if((rad_curr > radius_min) & (rad_curr < radius_max)){
                radiis_packing[counter] = rad_curr + random_rad;
                curr_counter++;
                counter++;
            }
    }
    }
    
    std::sort(radiis_packing.begin(),radiis_packing.end(),[](const double a, double  b) -> bool
    {
        return a> b;
    });

    return radiis_packing;

}