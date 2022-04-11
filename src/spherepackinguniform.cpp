#include "spherepackinguniform.h"
#include <algorithm>    // std::sort
#include <random>

using namespace std;
using namespace Eigen;

SpherePackingUniform::SpherePackingUniform(std::vector<unsigned> num_sph, std::vector<double> radii_spheres_, double icvf_,Eigen::Vector3d &min_l, Eigen::Vector3d &max_l)
{
    num_spheres = num_sph;
    nb_spheres=0;
    for (unsigned i =0; i<num_spheres.size(); i++){
        nb_spheres = nb_spheres +  num_spheres[i];
    }
    radii_spheres = radii_spheres_;
    icvf = icvf_;
    min_limits = min_l;
    max_limits = max_l;
    spheres.clear();
}


void SpherePackingUniform::computeMinimalSize(std::vector<double> radiis, double icvf_,Eigen::Vector3d& l){

   
    if(icvf_>= 0.7 && icvf_ < 0.99){
        icvf_+=0.01;
    }

    double vol = 0;

    for(uint i = 0 ; i < radiis.size();i++){
        vol+= 4./3.*radiis[i]*radiis[i]*radiis[i]*M_PI;
        
    }

    double l_ = fmax(pow(vol/icvf_, 1/3.), radiis[radiis.size()-1]*2/icvf_);
    l = {l_,l_,l_};
}


void SpherePackingUniform::createSubstrate()
{
    uint repetition = 40;
    uint max_adjustments = 5;
    double best_icvf = 0;
    vector<Sphere> best_spheres;
    Eigen::Vector3d best_max_limits;
    min_limits = {0.,0.,0.};

    std::random_device rd;

    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);
    std::vector<double> radiis(nb_spheres,0.0);

    bool achieved = false;
    unsigned nb_sphere_acc = 0;

    for (unsigned j=0; j<num_spheres.size(); j++){
        for(unsigned i=0; i< num_spheres[j]; i++){
            double random_rad = udist(gen) * 1e-10;
            radiis[nb_sphere_acc] = radii_spheres[j]*1e-3 + random_rad;
            nb_sphere_acc += 1;
        }

    }

    // using a lambda function:
    std::sort(radiis.begin(),radiis.end(),[](const double a, double  b) -> bool
    {
        return a> b;
    });


    uint adjustments = 0;
     // We increease 1% the total area. (Is prefered to fit all the spheres than achieve a perfect ICVF.)
    double adj_increase = icvf*0.01;
    while(!achieved){

        double target_icvf = this->icvf+adjustments*adj_increase;
        computeMinimalSize(radiis,target_icvf,max_limits);

        for(uint t = 0 ;  t < repetition; t++){
            vector<Sphere> spheres_to_add;

            spheres.clear();
            for(unsigned i = 0 ; i < nb_spheres; i++){
                unsigned stuck = 0;

                while(++stuck <= 1000){

                    double t = udist(gen);
                    double x = (t*max_limits[0]) + (1-t)*min_limits[0];
                    t = udist(gen);
                    double y = (t*max_limits[1] + (1-t)*min_limits[1]);
                    t = udist(gen);
                    double z = (t*max_limits[2] + (1-t)*min_limits[2]);

                    Vector3d P = {x,y,z};
                    Sphere sph(P,radiis[i]);


                    double min_distance;

                    bool collision = checkForCollition(sph,min_limits,max_limits,spheres_to_add,min_distance);

                    if(!collision){
                        for (unsigned j = 0; j < spheres_to_add.size(); j++)
                            spheres.push_back(spheres_to_add[j]);
                        break;
                    }
                }

                int dummy;
                double icvf_current = computeICVF(spheres,min_limits, max_limits,dummy);
                if(icvf_current > best_icvf ){
                    best_icvf = icvf_current;
                    best_spheres.clear();
                    best_spheres = spheres;
                    best_max_limits = max_limits;
                }
            } // end for spheres

            if(this->icvf - best_icvf  < 0.0005){
                achieved = true;
                break;
            }
        }

        
        spheres.clear();
        adjustments++;
        if(adjustments > max_adjustments){
            break;
        }
    }

    spheres = best_spheres;
    max_limits = best_max_limits;

    //TODO cambiar a INFO
    int perc_;
    double icvf_current = computeICVF(spheres,min_limits, max_limits,perc_);
    cout << "Percentage of spheres selected: "+ to_string(double(perc_)/radiis.size()*100.0)
            + "%,\nICVF achieved: " + to_string(icvf_current*100) + "  ("+ to_string( int((icvf_current/icvf*100))) + "% of the desired icvf)\n" << endl;
}

void SpherePackingUniform::printSubstrate(ostream &out)
{
    out << 1e-3 << endl;
    for(unsigned i = 0; i <spheres.size(); i++){

        out << spheres[i].P[0]*1e3 << " " << spheres[i].P[1]*1e3 << " " << spheres[i].P[2]*1e3 << " "
        << spheres[i].radius*1e3 << endl;;
    }
}

bool SpherePackingUniform::checkForCollition(Sphere sph, Vector3d min_limits, Vector3d max_limits, std::vector<Sphere>& spheres_to_add,double &min_distance)
{
    spheres_to_add.clear();

    checkBoundaryConditions(sph,spheres_to_add,min_limits,max_limits);

    min_distance = 1e10;

    bool collision = false;

    for(unsigned i = 0 ; i < spheres.size(); i++){
        for(unsigned j = 0 ; j < spheres_to_add.size(); j++){

            double distance = (spheres[i].P - spheres_to_add[j].P).norm();

            if(distance - (spheres[i].radius + spheres_to_add[j].radius) < 1e-15){
                min_distance = 0;
                collision = true;
                break;
            }
            if(distance < min_distance)
                min_distance = distance;
        }
    }

    // we need to check that the spheres to add don't interse4ct each other (very small voxel sizes)

    for(unsigned i = 0 ; i < spheres_to_add.size()-1; i++){
        for(unsigned j = i+1 ; j < spheres_to_add.size(); j++){

            double distance = (spheres_to_add[i].P - spheres_to_add[j].P).norm();

            if(distance - (spheres_to_add[i].radius + spheres_to_add[j].radius) < 1e-15){
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

double SpherePackingUniform::computeICVF(std::vector<Sphere>& spheres, Vector3d& min_limits, Vector3d& max_limits,int& num_no_repeat)
{
    if (spheres.size() == 0)
        return 0;

    double VolV = (max_limits[0] - min_limits[0])*(max_limits[1] - min_limits[1])*(max_limits[2] - min_limits[2]);
    double VolS = 0;

    // using a lambda function:
    std::sort(spheres.begin(),spheres.end(),[](const Sphere a, Sphere b) -> bool
    {
        return a.radius > b.radius;
    });

    double rad_holder = -1;
    num_no_repeat = 0;
    for (uint i = 0; i < spheres.size(); i++){

        if( fabs(rad_holder - spheres[i].radius) < 1e-15 ){
            continue;
        }
        else{
            rad_holder = spheres[i].radius;
        }

        double rad = spheres[i].radius;
        VolS += 4./3.*M_PI * rad * rad * rad;
        num_no_repeat++;
    }
    
    return VolS/VolV;
}


void SpherePackingUniform::checkBoundaryConditions(Sphere sph, std::vector<Sphere>& spheres_to_add, Vector3d min_limits, Vector3d max_limits){
    vector<Sphere> to_add;

    to_add.push_back(sph);

    for(int i = 0 ; i < 3; i++){

        double rad = sph.radius;

        if(sph.P[i]+rad >= max_limits[i]){

            Sphere tmp = sph;
            tmp.P[i]+= min_limits[i] - max_limits[i];
            to_add.push_back(tmp);
        }

        if(sph.P[i] - rad <= min_limits[i]){
            Sphere tmp = sph;
            tmp.P[i]+= max_limits[i] - min_limits[i];
            to_add.push_back(tmp);
        }
    }

    unsigned nb_s_first = to_add.size();
    if(nb_s_first >= 3){
        for(unsigned j = 1 ; j < nb_s_first; j++){
            Sphere jkr(to_add[j]);
            for(int i = 0 ; i < 3; i++){
                double rad = jkr.radius;

                if(jkr.P[i]+rad >= max_limits[i]){
                    Sphere tmp(jkr);
                    tmp.P[i]+= min_limits[i] - max_limits[i];
                    to_add.push_back(tmp);
                }

                if(jkr.P[i] - rad <= min_limits[i]){
                    Sphere tmp(jkr);
                    tmp.P[i]+= max_limits[i] - min_limits[i];
                    to_add.push_back(tmp);
                }
            }
        }
    }

    unsigned nb_s_second = to_add.size();
    if(nb_s_second > nb_s_first){
        for(unsigned j = nb_s_first-1 ; j < nb_s_second; j++){
            Sphere jkr(to_add[j]);
            for(int i = 0 ; i < 3; i++){
                double rad = jkr.radius;

                if(jkr.P[i]+rad >= max_limits[i]){
                    Sphere tmp(jkr);
                    tmp.P[i]+= min_limits[i] - max_limits[i];
                    to_add.push_back(tmp);
                }

                if(jkr.P[i] - rad <= min_limits[i]){
                    Sphere tmp(jkr);
                    tmp.P[i]+= max_limits[i] - min_limits[i];
                    to_add.push_back(tmp);
                }
            }
        }
    }

    unsigned nb_s_third = to_add.size();
    if(nb_s_third > nb_s_second){
        for(unsigned j = nb_s_second-1 ; j < nb_s_third; j++){
            Sphere jkr(to_add[j]);
            for(int i = 0 ; i < 3; i++){
                double rad = jkr.radius;

                if(jkr.P[i]+rad >= max_limits[i]){
                    Sphere tmp(jkr);
                    tmp.P[i]+= min_limits[i] - max_limits[i];
                    to_add.push_back(tmp);
                }

                if(jkr.P[i] - rad <= min_limits[i]){
                    Sphere tmp(jkr);
                    tmp.P[i]+= max_limits[i] - min_limits[i];
                    to_add.push_back(tmp);
                }
            }
        }
    }    

    
    for(unsigned i = 0 ; i < to_add.size(); i++)
    {
        bool rep = false;
        for( unsigned j = 0; j < spheres_to_add.size(); j++){
            if( (fabs(to_add[i].P[0] - spheres_to_add[j].P[0]) < 1e-12) && (fabs(to_add[i].P[1] - spheres_to_add[j].P[1]) < 1e-12)  && (fabs(to_add[i].P[2] - spheres_to_add[j].P[2]) < 1e-12) && (fabs(to_add[i].radius - spheres_to_add[j].radius) < 1e-14))
            {
                rep = true;
                break;
            }
        }

        if(rep == false){
            spheres_to_add.push_back(to_add[i]);
        }
    }
}
