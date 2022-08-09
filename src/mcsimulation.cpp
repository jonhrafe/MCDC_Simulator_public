#include "mcsimulation.h"
#include <Eigen/Dense>
#include "simerrno.h"
#include "pgsesequence.h"
#include "gradientwaveform.h"
#include <iostream>

int MCSimulation::count =0;

MCSimulation::MCSimulation()
{
    dynamicsEngine = NULL;
    dataSynth = NULL;
    id = count;
    count++;
}

/*DEPRECATED*/
MCSimulation::MCSimulation(std::string config_file)
{
    dynamicsEngine = NULL;
    dataSynth      = NULL;

    params.readSchemeFile(config_file);
    dynamicsEngine = new DynamicsSimulation(params);

    if(params.scheme_file.length() > 2){
        scheme.readSchemeFile(params.scheme_file,params.scale_from_stu);
    }

    if(scheme.type == "PGSE"){
        dataSynth = new PGSESequence(scheme);
        dataSynth->setNumberOfSteps(dynamicsEngine->params.num_steps);

        if(params.subdivision_flag){
            dataSynth->subdivision_flag = true;
            dataSynth->subdivisions = params.subdivisions;
            dataSynth->initializeSubdivisionSignals();
        }
    }

    dynamicsEngine->id = count;
    id = count;
    count++;
}

MCSimulation::MCSimulation(Parameters& params_)
{
    dynamicsEngine = NULL;
    dataSynth = NULL;

    params = params_;

    dynamicsEngine = new DynamicsSimulation(params);
  

    if(params.scheme_file.length() > 2){
        scheme.readSchemeFile(params.scheme_file,params.scale_from_stu);
    }

    if(scheme.type == "PGSE"){
        dataSynth = new PGSESequence(scheme);
    }
    if(scheme.type == "WAVEFORM"){
        dataSynth = new GradientWaveform(scheme);
    }



    if (dataSynth){
        dataSynth->setNumberOfSteps(dynamicsEngine->params.num_steps);
    }
    
    if(params.subdivision_flag){
        dataSynth->subdivision_flag = true;
        dataSynth->subdivisions = params.subdivisions;
        dataSynth->initializeSubdivisionSignals();
    }

    dynamicsEngine->id = count;
    id = count;
    count++;
}


void MCSimulation::startSimulation()
{

    iniObstacles();

    if(dataSynth != NULL){
        dynamicsEngine->startSimulation(dataSynth);
    }
    else{
        dynamicsEngine->startSimulation();
    }

}

double MCSimulation::getExpectedFreeeDecay(unsigned i)
{
    if(dataSynth){
        double b = dataSynth->getbValue(i);
        return exp(-b*params.diffusivity_extra);
    }

    return -1;
}


void MCSimulation::iniObstacles()
{
    addCylindersObstaclesFromFiles();

    addPLYObstaclesFromFiles();

    addVoxels();

    addCylindersConfigurations();
    //Used only if there's a voxel (deprecated)
    //addExtraObstacles();

    addSpheresObstaclesFromFiles();

}


//* Auxiliare method to split words in a line using the spaces*//
template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void MCSimulation::addCylindersObstaclesFromFiles()
{
    for(unsigned i = 0; i < params.cylinders_files.size(); i++){

        bool z_flag = false;
        std::ifstream in(params.cylinders_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}

            std::vector<std::string> jkr = split(line,' ');
            if (jkr.size() != 7){
                z_flag = true;
                //std::cout << "\033[1;33m[Warning]\033[0m Cylinder orientation was set towards the Z direction by default" << std::endl;
            }
            break;
        }
        in.close();

        // Permeability file - if any
        double perm_; 

        std::ifstream in_perm;
        if(params.cylinder_permeability_files.size() >0){
            in_perm.open(params.cylinder_permeability_files[i]);
        }

        // Diffusion coefficients
        double diff_i; 
        double diff_e;
        
        in.open(params.cylinders_files[i]);


        if(z_flag){
            double x,y,z,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> r)
            {
                Cylinder cyl(Cylinder(Eigen::Vector3d(x,y,z),Eigen::Vector3d(x,y,z+1.0),r,scale, perm_));

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
        
                cyl.setPercolation(perm_);

                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                cyl.setDiffusion(diff_i, diff_e);

                dynamicsEngine->cylinders_list.push_back(cyl);
            }
            in_perm.close();
            in.close();
        }
        else{
            double x,y,z,ox,oy,oz,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> ox >> oy >> oz >> r)
            {

                Cylinder cyl(Eigen::Vector3d(x,y,z),Eigen::Vector3d(ox,oy,oz),r,scale, perm_);

                // Local permeability - Different for each obstacle
                if(in_perm){
                    in_perm >> perm_;
                }
                // Global permeability - Same for all obstacle
                else{
                    perm_ = params.obstacle_permeability;
                }  
        
                cyl.setPercolation(perm_);

                // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
                diff_i = params.diffusivity_intra; 
                diff_e = params.diffusivity_extra;
                cyl.setDiffusion(diff_i, diff_e);

                dynamicsEngine->cylinders_list.push_back(cyl);
            }
            in_perm.close();
            in.close();
        }

    }
}


void MCSimulation::addPLYObstaclesFromFiles()
{
    for(unsigned i = 0; i < params.PLY_files.size(); i++){

        PLYObstacle ply_(params.PLY_files[i],params.PLY_scales[i]);

        // Permeability - Kept outside initialization to be consistent with cylinders and spheres. Easily moved to ply constructor. 
        double perm_; 
        perm_ = params.PLY_permeability[i];
        ply_.setPercolation(perm_);

        // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
        double diff_i; 
        double diff_e;
        
        diff_i = params.diffusivity_intra; 
        diff_e = params.diffusivity_extra;
        ply_.setDiffusion(diff_i, diff_e);

        // Add PLY to list
        dynamicsEngine->plyObstacles_list.push_back(ply_);
    }
}

void MCSimulation::addVoxels()
{
    for(unsigned i = 0 ; i < params.voxels_list.size(); i++){
        dynamicsEngine->voxels_list.push_back(Voxel(params.voxels_list[i].first,params.voxels_list[i].second));
    }
}

void MCSimulation::addCylindersConfigurations()
{

    if(params.hex_packing){
        double rad = params.hex_packing_radius,sep = params.hex_packing_separation;

        // h = sqrt(3)/2 * sep
        double h = 0.866025404*sep;

        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(0,0,0),Eigen::Vector3d(0,0,1.0),rad));
        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(sep,0,0),Eigen::Vector3d(sep,0,1.0),rad));

        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(0,2.0*h,0),Eigen::Vector3d(0,2.0*h,1.0),rad));
        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(sep,2.0*h,0),Eigen::Vector3d(sep,2.0*h,1.0),rad));

        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(0.5*sep,h,0),Eigen::Vector3d(0.5*sep,h,1.0),rad));

        // To avoid problems with the boundaries
        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(-0.5*sep,h,0),Eigen::Vector3d(-0.5*sep,h,1.0),rad));
        dynamicsEngine->cylinders_list.push_back(Cylinder(Eigen::Vector3d(1.5*sep,h,0),Eigen::Vector3d(1.5*sep,h,1.0),rad));

        if(dynamicsEngine->voxels_list.size()>0)
            dynamicsEngine->voxels_list.clear();

        dynamicsEngine->voxels_list.push_back(Voxel(Eigen::Vector3d(0,0,0),Eigen::Vector3d(sep,2.0*h,2.0*h)));

    }
}

void MCSimulation::addSpheresObstaclesFromFiles()
{
    for(unsigned i = 0; i < params.spheres_files.size(); i++){

        std::ifstream in(params.spheres_files[i]);

        if(!in){
            return;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {first-=1;continue;}
            break;
        }
        in.close();

        // Permeability file - if any
        double perm_; 

        std::ifstream in_perm;
        if(params.sphere_permeability_files.size() >0){
            in_perm.open(params.sphere_permeability_files[i]);
        }

        // Diffusion coefficients
        double diff_i; 
        double diff_e;
            
        in.open(params.spheres_files[i]);
        double x,y,z,r;
        double scale;
        in >> scale;

        while (in >> x >> y >> z >> r)
        {
            Sphere sph(Eigen::Vector3d(x,y,z),r,scale);

            // Local permeability - Different for each obstacle
            if(in_perm.is_open()){
                in_perm >> perm_;
            }
            // Global permeability - Same for all obstacle
            else{
                perm_ = params.obstacle_permeability;
            }            

            sph.setPercolation(perm_);

            // Diffusion coefficient - Useless now, to be implemented for obstacle specific Di
            diff_i = params.diffusivity_intra; 
            diff_e = params.diffusivity_extra;
            sph.setDiffusion(diff_i, diff_e);
            
            // Add sphere to list
            dynamicsEngine->spheres_list.push_back(sph);     
        }
        in.close();
        in_perm.close();
    }
}

bool cylinderIsCloseBoundery(Cylinder& cyl, Eigen::Vector3d min_limits,Eigen::Vector3d max_limits){

    //gap to the boundary
    double gap = 1e-6;
    //3 dimensional vector
    for (int i = 0 ; i < 3; i++)
        if( (cyl.P[i] - cyl.radius - gap < min_limits[i]) || (cyl.P[i] + cyl.radius + gap  > max_limits[i]) )
            return true;

    return false;
}

void MCSimulation::addExtraObstacles()
{
    if(dynamicsEngine->voxels_list.size() == 0)
        return;

    std::vector<Eigen::Vector3d> multipliers;

    Eigen::Vector3d gap = params.max_limits - params.min_limits;

    for(int i = -1  ;i <= 1; i++)
        for(int j = -1  ;j <= 1; j++)
            for(int k = -1 ;k <= 1; k++){
                Eigen::Vector3d jkr(i*gap[0],j*gap[1],k*gap[2]);
                multipliers.push_back(jkr);
            }


    unsigned long cylinders_num = dynamicsEngine->cylinders_list.size();

    for (unsigned c = 0; c < cylinders_num ;c++)
        for (unsigned i = 0 ; i < multipliers.size(); i++)
            if(multipliers[i][0]!=0.0 || multipliers[i][1]!=0.0 || multipliers[i][2]!=0.0)
            {
                Eigen::Vector3d P_ = dynamicsEngine->cylinders_list[c].P;
                Eigen::Vector3d Q_ = dynamicsEngine->cylinders_list[c].Q;
                P_[0]+= multipliers[i][0];P_[1]+=multipliers[i][1];P_[2]+=multipliers[i][2];
                Q_[0]+= multipliers[i][0];Q_[1]+=multipliers[i][1];Q_[2]+=multipliers[i][2];
                Cylinder tmp_cyl(P_,Q_,dynamicsEngine->cylinders_list[c].radius);


                //if the obstacle is close enough
                //if (cylinderIsCloseBoundery(tmp_cyl,params.min_limits,params.max_limits))
                    dynamicsEngine->cylinders_list.push_back(tmp_cyl);
            }

}


MCSimulation::~MCSimulation()
{
    if(dynamicsEngine != NULL)
        delete dynamicsEngine;

    if(dataSynth != NULL)
        delete dataSynth;
}


