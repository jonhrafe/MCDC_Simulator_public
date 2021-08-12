#include "cylindergammadistribution.h"
#include <algorithm>    // std::sort
#include <random>
#include  "simerrno.h"

using namespace std;
using namespace Eigen;

CylinderGammaDistribution::CylinderGammaDistribution(unsigned num_cyl, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l, float min_radius)
{
    num_obstacles = num_cyl;
    alpha = a;
    beta  = b;
    icvf = icvf_;
    min_limits = min_l;
    max_limits = max_l;
    cylinders.clear();
    this->min_radius = min_radius;
}

void CylinderGammaDistribution::computeMinimalSize(std::vector<double> radiis, double icvf_,Eigen::Vector3d& l){

    /*A little heuristic for complicated ICVF: > 0.7*/
    if(icvf_>= 0.7 && icvf_ < 0.99){
        icvf_+=0.01;
    }

    double area = 0;

    for(uint i = 0 ; i < radiis.size();i++){
        area+= radiis[i]*radiis[i]*M_PI;
    }

    double l_ = sqrt(area/icvf_);

    l = {l_,l_,l_};
}


void CylinderGammaDistribution::displayGammaDistribution()
{
    const int nrolls=10000;  // number of experiments
    const int nstars=100;    // maximum number of stars to distribute
    string message;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha,beta);

    int p[11]={};

    for (int i=0; i<nrolls; ++i) {
        double number = distribution(generator);

        if (number<10) ++p[int(number)];
        else ++p[10];
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

void CylinderGammaDistribution::createGammaSubstrate()
{
    // generate the gamma distribution
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::gamma_distribution<double> distribution(alpha,beta);
    uint repetition = 40;
    uint max_adjustments = 5;
    double best_icvf = 0;
    vector<Cylinder> best_cylinders;
    Eigen::Vector3d best_max_limits;
    min_limits = {0.,0.,0.};

    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);
    std::vector<double> radiis(num_obstacles,0);

    bool achieved = false;

    int tried = 0;

    for (unsigned i=0; i< num_obstacles; ++i) {

        if(tried > 10000){
            string message = " Radii distribution cannot be sampled [Min. radius Error]\n";
            SimErrno::error(message,cout);
            assert(0);
        }
        double jkr =  distribution(generator);

        if(jkr< this->min_radius){
            i--;
            tried++;
            continue;
        }
        tried=0;

        radiis[i] = jkr*1e-3; //WE CONVERT FROM UM TO MM HERE
    }

    // using a lambda function:
    std::sort(radiis.begin(),radiis.end(),[](const double a, double  b) -> bool
    {
        return a> b;
    });


    uint adjustments = 0;
     // We increease 1% the total area. (Is prefered to fit all the cylinders than achieve a perfect ICVF.)
    double adj_increase = icvf*0.01;
    while(!achieved){

        double target_icvf = this->icvf+adjustments*adj_increase;
        computeMinimalSize(radiis,target_icvf,max_limits);

        for(uint t = 0 ;  t < repetition; t++){
            vector<Cylinder> cylinders_to_add;

            cylinders.clear();
            for(unsigned i = 0 ; i < num_obstacles; i++){
                unsigned stuck = 0;

                while(++stuck <= 1000){

                    double t = udist(gen);
                    double x = (t*max_limits[0]) + (1-t)*min_limits[0];
                    t = udist(gen);
                    double y = (t*max_limits[1] + (1-t)*min_limits[1]);
                    double z = 0;

                    Vector3d Q = {x,y,z};
                    Vector3d D = {x,y,z+1};
                    Cylinder cyl(Q,D,radiis[i]);


                    double min_distance;

                    bool collision = checkForCollition(cyl,min_limits,max_limits,cylinders_to_add,min_distance);

                    if(!collision){
                        for (unsigned j = 0; j < cylinders_to_add.size(); j++)
                            cylinders.push_back(cylinders_to_add[j]);
                        break;
                    }
                }

                int dummy;
                double icvf_current = computeICVF(cylinders,min_limits, max_limits,dummy);
                if(icvf_current > best_icvf ){
                    best_icvf = icvf_current;
                    best_cylinders.clear();
                    best_cylinders = cylinders;
                    best_max_limits = max_limits;
                }
            } // end for cylinders

            if(this->icvf - best_icvf  < 0.0005){
                achieved = true;
                break;
            }
        }
        cylinders.clear();
        adjustments++;
        cout << best_icvf << endl;
        if(adjustments > max_adjustments){
            break;
        }
    }

    cylinders = best_cylinders;
    max_limits = best_max_limits;

    //TODO cambiar a INFO
    int perc_;
    double icvf_current = computeICVF(cylinders,min_limits, max_limits,perc_);

     string  message = "Percentage of cylinders selected: "+ to_string(double(perc_)/radiis.size()*100.0)
            + "%,\nICVF achieved: " + to_string(icvf_current*100) + "  ("+ to_string( int((icvf_current/icvf*100))) + "% of the desired icvf)\n";
     SimErrno::info(message,cout);

}

void CylinderGammaDistribution::printSubstrate(ostream &out)
{
    out << 1e-3 << endl;
    for(unsigned i = 0; i < cylinders.size(); i++){

        out << cylinders[i].P[0]*1e3 << " " << cylinders[i].P[1]*1e3 << " " << cylinders[i].P[2]*1e3 << " "
                                     << cylinders[i].Q[0]*1e3 << " " << cylinders[i].Q[1]*1e3 << " " << cylinders[i].Q[2]*1e3 << " "
                                     << cylinders[i].radius*1e3 << endl;
    }
}

bool CylinderGammaDistribution::checkForCollition(  Cylinder cyl, Vector3d min_limits, Vector3d max_limits, std::vector<Cylinder>& cylinders_to_add,double &min_distance)
{
    cylinders_to_add.clear();

    checkBoundaryConditions(cyl,cylinders_to_add,min_limits,max_limits);

    min_distance = 1e10;

    bool collision = false;

    for(unsigned i = 0 ; i < cylinders.size(); i++){
        for(unsigned j = 0 ; j < cylinders_to_add.size(); j++){

            double distance = (cylinders[i].P - cylinders_to_add[j].P).norm();

            if(distance - (cylinders[i].radius + cylinders_to_add[j].radius) < 1e-15){
                min_distance = 0;
                collision = true;
                break;
            }
            if(distance < min_distance)
                min_distance = distance;
        }
    }

    // we need to check that the cylinders to add don't interse4ct each other (very small voxel sizes)

    for(unsigned i = 0 ; i < cylinders_to_add.size()-1; i++){
        for(unsigned j = i+1 ; j < cylinders_to_add.size(); j++){

            double distance = (cylinders_to_add[i].P - cylinders_to_add[j].P).norm();

            if(distance - (cylinders_to_add[i].radius + cylinders_to_add[j].radius) < 1e-15){
                min_distance = 0;
                collision = true;
                break;
            }
            if(distance < min_distance)
                min_distance = distance;
        }
    }

    return collision;

}

/*
WARNING: The way we discard repeated cylinders is using radius. Repreated radius (like really the same)
are considered the same. This becasuse we don't track wich cylinders had to be replciated to mantain the voxel
symmetry
*/

double CylinderGammaDistribution::computeICVF(std::vector<Cylinder>& cylinders, Vector3d& min_limits, Vector3d& max_limits,int& num_no_repeat)
{
    if (cylinders.size() == 0)
        return 0;

    double AreaV = (max_limits[0] - min_limits[0])*(max_limits[1] - min_limits[1]);

    double AreaC = 0;

    // using a lambda function:
    std::sort(cylinders.begin(),cylinders.end(),[](const Cylinder a, Cylinder b) -> bool
    {
        return a.radius > b.radius;
    });

    double rad_holder = -1;
    num_no_repeat = 0;
    for (uint i = 0; i < cylinders.size(); i++){

        if( fabs(rad_holder - cylinders[i].radius) < 1e-15 ){
            continue;
        }
        else{
            rad_holder = cylinders[i].radius;
        }

        double rad = cylinders[i].radius;
        AreaC += M_PI*rad*rad;
        num_no_repeat++;
    }
    return AreaC/AreaV;
}


void CylinderGammaDistribution::checkBoundaryConditions(Cylinder cyl, std::vector<Cylinder>& cylinders_to_add, Vector3d min_limits, Vector3d max_limits){
    vector<Cylinder> to_add;

    to_add.push_back(cyl);

    for(int i = 0 ; i < 2; i++){

        double rad = cyl.radius;

        if(cyl.P[i]+rad >= max_limits[i]){

            Cylinder tmp = cyl;
            tmp.P[i]+= min_limits[i] - max_limits[i];
            tmp.Q[i]+= min_limits[i] - max_limits[i];
            to_add.push_back(tmp);
        }

        if(cyl.P[i] - rad <= min_limits[i]){
            Cylinder tmp = cyl;
            tmp.P[i]+= max_limits[i] - min_limits[i];
            tmp.Q[i]+= max_limits[i] - min_limits[i];
            to_add.push_back(tmp);
        }
    }

    if(to_add.size() == 3)
        for(unsigned j = 1 ; j < 3; j++){
            Cylinder jkr(to_add[j]);
            for(int i = 0 ; i < 2; i++){
                double rad = cyl.radius;

                if(jkr.P[i]+rad >= max_limits[i]){
                    Cylinder tmp(jkr);
                    tmp.P[i]+= min_limits[i] - max_limits[i];
                    tmp.Q[i]+= min_limits[i] - max_limits[i];
                    to_add.push_back(tmp);
                }

                if(jkr.P[i] - rad <= min_limits[i]){
                    Cylinder tmp(jkr);
                    tmp.P[i]+= max_limits[i] - min_limits[i];
                    tmp.Q[i]+= max_limits[i] - min_limits[i];
                    to_add.push_back(tmp);
                }
            }
        }

    for(unsigned i = 0 ; i < to_add.size(); i++)
    {
        bool rep = false;
        for( unsigned j = 0; j < cylinders_to_add.size(); j++){
            if( ( fabs(to_add[i].P[0] - cylinders_to_add[j].P[0]) < 1e-12) && (fabs(to_add[i].P[1] - cylinders_to_add[j].P[1]) < 1e-12))
            {
                rep = true;
                break;
            }
        }

        if(rep == false){
            cylinders_to_add.push_back(to_add[i]);
        }
    }
}
