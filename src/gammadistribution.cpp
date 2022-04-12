#include "gammadistribution.h"
#include <algorithm>    // std::sort
#include <random>

using namespace std;
using namespace Eigen;


GammaDistribution::GammaDistribution(){}


GammaDistribution::GammaDistribution(std::vector<unsigned> num_obs, std::vector<double> a, std::vector<double> b, double r_min, double r_max)
{
    num_obstacles = num_obs; 
    radius_min = r_min;
    radius_max =r_max;
    alpha = a;
    beta  = b;
}

void GammaDistribution::displayDistribution()
{
    const int nrolls=10000;  // number of experiments
    const int nstars=100;    // maximum number of stars to distribute
    unsigned tot_num_obstacles = 0; 
    for(int i=0;i<num_obstacles.size(); ++i){
        tot_num_obstacles += num_obstacles[i];
    }

    string message;
    std::random_device rd;
    std::default_random_engine generator(rd());
    
    int p[11]={};  


    for(int i=0;i<num_obstacles.size(); ++i){
        unsigned curr_counter    = 0; 
        double curr_alpha   = alpha[i];
        double curr_beta    = beta[i];
        double curr_num_obstacle = num_obstacles[i];

        unsigned curr_nb_roll   = nrolls * num_obstacles[i]/tot_num_obstacles; 

        std::gamma_distribution<double> distribution(curr_alpha,curr_beta);


        for (int i=0; i<curr_nb_roll; ++i) {
            double number = distribution(generator);
            if (number<10) ++p[int(number)];
            else ++p[10];
        }
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


std::vector<double> GammaDistribution::createRadiiList(){
    
    std::random_device rd;
    std::default_random_engine generator(rd());


    unsigned tot_num_obstacles = 0; 
    for(int i=0;i<num_obstacles.size(); ++i){
        tot_num_obstacles += num_obstacles[i];
    }
   
    std::vector<double> radiis_packing(tot_num_obstacles,0.0);

    unsigned counter    = 0; 
    double rad_curr     = 0.0;

    for(int i=0;i<num_obstacles.size(); ++i){
        unsigned curr_counter    = 0; 
        double curr_alpha   = alpha[i];
        double curr_beta    = beta[i];
        double curr_num_obstacle = num_obstacles[i];

        std::gamma_distribution<double> distribution(curr_alpha,curr_beta);

        while(curr_counter < curr_num_obstacle){
            rad_curr = distribution(generator)*1e-3;
            if((rad_curr > radius_min) & (rad_curr < radius_max)){
                radiis_packing[counter] = rad_curr;
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