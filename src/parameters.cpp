#include "parameters.h"
#include <fstream>
#include <iostream>
#include "constants.h"
#include "simerrno.h"
using namespace std;

Parameters::Parameters()
{
    //Dummy initial values;
    scale_from_stu      = false;
    seed                = -1;
    save_phase_shift    = false;
    write_traj          = false;
    write_txt           = false;
    write_bin           =  true;

    hex_cyl_packing    = false;
    hex_sphere_packing = false;
    hex_packing_radius      = 0;
    hex_packing_separation  = 0;

    gamma_cyl_packing = false;
    gamma_sph_packing = false;
    gamma_packing_alpha = 0;
    gamma_packing_beta  = 0;
    gamma_num_obstacles = 0;
    gamma_icvf          = 0;
    min_obstacle_radii  = 0;

    ini_walkers_file = "";
    num_proc    = 0;
    verbatim    = false;

    record_phase_times.clear();
    record_pos_times.clear();
    subdivisions_file = "";

    computeVolume = false;
    custom_sampling_area = false;
    separate_signals = false;
    img_signal = false;
    hex_sphere_packing = false;

    for (auto i= 0;i<3; i++)
        min_sampling_area[i] = max_sampling_area[i] = 0.0;
}

void Parameters::readSchemeFile(std::string conf_file_path)
{

    ifstream in(conf_file_path);

    if(!in){
        cout << "[ERROR] Configuration file" << endl;
        return;
    }

    string tmp="";
    while((in >> tmp) && (str_dist(tmp,"<END>") >= 2) ){

        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"n") == 0){
            in >> num_walkers;
        }
        else if(str_dist(tmp,"t") == 0){
            in >> num_steps;
        }
        else if(str_dist(tmp,"duration") <= 1){
            in >> sim_duration;
        }
        else if(str_dist(tmp,"scheme_file") <= 1){
            in >> scheme_file;
        }
        else if(str_dist(tmp,"diffusivity") <= 1){
            in >> diffusivity;
        }
        else if( (str_dist(tmp,"out_traj_file_index") <= 2) or (str_dist(tmp,"exp_prefix") <= 2)) {
            in >> traj_file;
            output_base_name = traj_file;
        }
        else if( str_dist(tmp,"ini_walkers_file") <= 3){
            in >> ini_walkers_file;
        }
        else if(str_dist(tmp,"write_txt") <= 1){
            in >> write_txt;
        }
        else if(str_dist(tmp,"write_bin") <= 1){
            in >> write_bin;
        }
        else if(str_dist(tmp,"write_traj_file") <= 2){
            in >> write_traj;
        }
        else if(str_dist(tmp,"scale_from_stu") <= 2){
            in >> scale_from_stu;
        }
        else if(str_dist(tmp,"seed") <= 1){
            in >> seed;
        }
        else if(str_dist(tmp,"<obstacle>") == 0){
            readObstacles(in);
        }
        else if ((str_dist(tmp,"<voxels>") == 0) or (str_dist(tmp,"<voxel>") == 0)){
            readVoxels(in);
        }
        else if(str_dist(tmp,"num_process") <= 1 || str_dist(tmp,"processors") <= 1){
            in >> num_proc;
        }
        else if(str_dist(tmp,"<log>") == 0){
            readInfoGatheringParams(in);
        }
        else if(str_dist(tmp,"<sampling_area>") ==0 || str_dist(tmp,"<spawning_area>") ==0){
            float tmp;
            in >> tmp; min_sampling_area[0] = tmp;
            in >> tmp; min_sampling_area[1] = tmp;
            in >> tmp; min_sampling_area[2] = tmp;
            in >> tmp; max_sampling_area[0] = tmp;
            in >> tmp; max_sampling_area[1] = tmp;
            in >> tmp; max_sampling_area[2] = tmp;

            custom_sampling_area = true;
        }
        else if(str_dist(tmp,"<delta>") == 0)
        {
            float tmp;
            //Save the three positions
            in >> tmp;ini_delta_pos.push_back(tmp);
            in >> tmp;ini_delta_pos.push_back(tmp);
            in >> tmp;ini_delta_pos.push_back(tmp);
        }
        else if(str_dist(tmp,"ini_walkers_pos") <= 2)
        {
            in >> ini_walker_flag;
            std::transform(ini_walker_flag.begin(), ini_walker_flag.end(), ini_walker_flag.begin(), ::tolower);

            if(str_dist(ini_walker_flag,"intra") > 1 && str_dist(ini_walker_flag,"extra") > 1 && str_dist(ini_walker_flag,"delta") > 1){
                SimErrno::error("Wrong walker's initialization position (Did you mean: 'walker_ini_file'?) ",cout);

                assert(0);
            }
        }
        else if(str_dist(tmp,"subdivisions_number") <= 2)
        {
            in >> number_subdivisions;
            subdivision_flag = (number_subdivisions>1)?true:false;

            if(number_subdivisions > 500 || number_subdivisions <=0 ){
                SimErrno::error("Unrealistic number of resulting subdivision voxels : " + std::to_string(number_subdivisions) + "^3",cout);

                assert(0);
            }
        }
        else if(str_dist(tmp,"subdivisions_file") <= 1)
        {
            in >> subdivisions_file;
            readSubdivisionFile();
            subdivision_flag |= (subdivisions.size()>0);
        }
        else if((str_dist(tmp,"permeability") <= 1) || (str_dist(tmp,"obstacle_permeability") <= 2))
        {
            in >> obstacle_permeability;
        }
        else if( str_dist(tmp,"coll_sphere_dist") <= 2 || str_dist(tmp,"colision_dist") <= 2 || str_dist(tmp,"sphere_size") <= 2)
        {
            in >> collision_sphere_distance;
        }
        else if( str_dist(tmp,"verbatim") <= 2 )
        {
           this->verbatim = true;
        }
        else if( str_dist(tmp,"log_opp") <= 1 )
        {
           this->log_opp = true;
        }
        else if( str_dist(tmp,"compute_volume") <= 1 )
        {
           this->computeVolume = true;
        }
        else if( str_dist(tmp,"separate_signals") <= 1 )
        {
           this->separate_signals = true;
           this->computeVolume = true;
        }
        else if( str_dist(tmp,"log_phase_shift") <= 2 )
        {
           this->log_phase_shift = true;
        }
        else if(str_dist(tmp,"max_sim_time") <= 2 || str_dist(tmp,"max_simulation_time") <= 3){
            in >> max_simulation_time;

            if(max_simulation_time<=0){
                SimErrno::error("Max simulation time must be greater than 1 second. ",cout);
                assert(0);
            }
        }
        else if( str_dist(tmp,"deportation") <= 2 )
        {
            in >> discard_illegals;
        }
        else if( str_dist(tmp,"discard_stucks") <= 2 )
        {
            in >> discard_stucks;
        }
        else{
            if( str_dist(tmp.substr(0,2),"</") > 0 )
                SimErrno::warning("Parameter: " + tmp + " Unknown",cout);
        }
    }

    if(scale_from_stu){
        //m^2/s to mm^2/ms
        diffusivity*=m2_to_mm2/s_to_ms;
        //seconds to ms
        sim_duration*=s_to_ms;
    }



    in.close();

    return;
}

//Set Methods

void Parameters::setNumWalkers(unsigned N)
{
    num_walkers = N;
}

void Parameters::setNumSteps(unsigned T)
{
    num_steps = T;
}

void Parameters::setDiffusivity(double D)
{
    diffusivity = D;
}

void Parameters::setSimDuration(double duration)
{
    sim_duration = duration;
}

void Parameters::setWriteTrajFlag(bool write_bin)
{
    write_traj = write_bin;
}

void Parameters::setWriteTextFlag(bool write_txt_)
{
    write_txt = write_txt_;
}

void Parameters::setMinLimits(Eigen::Vector3d min_limits_)
{
    min_limits = min_limits_;

}

void Parameters::setMaxLimits(Eigen::Vector3d max_limits_)
{
    max_limits = max_limits_;
}

void Parameters::setTrajFileName(std::string traj_file_)
{
    traj_file = traj_file_;
}

void Parameters::setOutputBaseFileName(std::string output_base_name_)
{
    output_base_name = output_base_name_;
}

void Parameters::iniWalkersFileName(std::string ini_walkers_file_)
{
    ini_walkers_file = ini_walkers_file_;
}

void Parameters::setSchemeFileName(std::string scheme_file_)
{
    scheme_file = scheme_file_;
}


//GET METHODS

unsigned Parameters::getNumWalkers()
{
    return num_walkers;
}

unsigned Parameters::getNumSteps()
{
    return num_steps;
}

double Parameters::getDiffusivity()
{
    return diffusivity;
}

bool Parameters::getWriteTrajFlag()
{
    return write_traj;
}

bool Parameters::getWriteTextFlag()
{
    return write_txt;
}

Eigen::Vector3d Parameters::getMinLimits()
{
    return min_limits;
}

Eigen::Vector3d Parameters::getMaxLimits()
{
    return max_limits;
}

std::string Parameters::getTrajFileName()
{
    return traj_file;
}

std::string Parameters::getOutputBaseFileName()
{
    return output_base_name;
}

std::string Parameters::getIniWalkersFileName()
{
    return ini_walkers_file;
}

std::string Parameters::getSchemeFileName()
{
    return scheme_file;
}


void Parameters::readObstacles(ifstream& in)
{
    string tmp="";
    unsigned num_obstacles = 0;

    while( !(str_dist(tmp,"</obstacle>") <= 2)){
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"cylinders_list") <= 2){
            string path;
            in >> path;
            cylinders_files.push_back(path);
            num_obstacles++;
        }
        if(str_dist(tmp,"spheres_list") <= 2){
            string path;
            in >> path;
            spheres_files.push_back(path);
            num_obstacles++;
        }
        if(str_dist(tmp,"oriented_cylinders_list") <= 2){
            string path;
            in >> path;
            cylinders_files.push_back(path);
            num_obstacles++;
        }
        if(str_dist(tmp,"ply") <= 1){
            string path;
            in >> path;
            PLY_files.push_back(path);
            PLY_percolation.push_back(0);
            num_obstacles++;
        }
        if(str_dist(tmp,"ply_scale") <= 1){
            double scale_;
            in >> scale_;
            PLY_scales.push_back(scale_);
        }
        if(str_dist(tmp,"ply_file_list") <= 2){
            string path;
            in >> path;
            readPLYFileList(path);
            num_obstacles++;
        }
        if(str_dist(tmp,"ply_file_list_scale_permeability") <= 3){
            string path;
            in >> path;
            readPLYFileListScalePercolation(path);
            num_obstacles++;
        }
        if(str_dist(tmp,"<cylinder_hex_packing>") <=1){
            this->hex_cyl_packing = true;
            readHexagonalParams(in);
            num_obstacles++;
        }
        if(str_dist(tmp,"<sphere_hex_packing>") <=1){
            this->hex_sphere_packing = true;
            readHexagonalParams(in);
            num_obstacles++;
        }
        if(str_dist(tmp,"<cylinder_gamma_packing>") <=1){
            gamma_cyl_packing = true;
            readGammaParams(in);
            num_obstacles++;
        }
        if(str_dist(tmp,"<sphere_gamma_packing>") <=1){
            gamma_sph_packing = true;
            readGammaParams(in);
            num_obstacles++;
        }
    }

    if(num_obstacles ==0){
        SimErrno::warning("<obstacle> tag initialized, but no valid obstacle tag found",cout);
    }

}

void Parameters::readVoxels(ifstream& in)
{
    string tmp="";
    double x,y,z;

    in >> x >> y >> z;
    min_limits = {x,y,z};
    in >> x >> y >> z;
    max_limits = {x,y,z};

    pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(min_limits,max_limits);
    voxels_list.push_back(voxel_);

    in >> tmp;
}

int Parameters::str_dist(string s, string t)
{
    ulong len_s = s.length();
    ulong len_t = t.length();

    /* base case: empty strings */
    if (len_s == 0) return int(len_t);
    if (len_t == 0) return int(len_s);

    if(len_s == 1 && len_t ==1)
        return s[0] != t[0];

    Eigen::MatrixXd costos(len_s,len_t);

    for(unsigned i = 0 ; i < s.size(); i++){
        for (unsigned j = 0 ; j < t.size(); j++){
            costos(i,j) = 0;
            costos(0,j) = j;
        }
        costos(i,0) = i;
    }

    int cost;

    for(unsigned i = 1 ; i < s.size(); i++){
        for (unsigned j = 1 ; j < t.size(); j++){
            /* test if last characters of the strings match */
            if (s[i] == t[j])
                cost = 0;
            else
                cost = 1;

            /* return minimum of delete char from s, delete char from t, and delete char from both */
            costos(i,j) =  min(min( costos(i-1,j) + 1,
                                    costos(i,j-1) + 1),
                               costos(i-1,j-1) + cost);
        }
    }

    return costos(s.length()-1,t.length()-1);
}

void Parameters::readInfoGatheringParams(ifstream& in)
{
    string tmp="";
    while(str_dist(tmp,"</log>"))
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"<positions>") <= 3)
        {
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

            while(str_dist(tmp,"</positions>"))
            {
                if(tmp.compare("t") == 0){
                    record_pos_times.push_back(num_steps);
                }
                else{
                    unsigned long time = stoul(tmp);
                    record_pos_times.push_back(unsigned(time));
                }

                in >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            }
        }
        else if(str_dist(tmp,"<phase>") <= 2)
        {
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

            while(str_dist(tmp,"</phase>"))
            {
                if(tmp.compare("t") == 0){
                    record_phase_times.push_back(num_steps);
                }
                else{
                    unsigned time = unsigned(stoul(tmp));
                    record_phase_times.push_back(time);
                }
                in >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
            }
        }
        else if(str_dist(tmp,"<propagator>")<= 2){
            in >> tmp;
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

            while(str_dist(tmp,"</propagator>"))
            {
                if(str_dist(tmp,"directions")<2){
                    std::string dir_path;
                    in >> dir_path;
                    readPropagatorDirections(dir_path);
                    // Read directions
                }
                else if(tmp.compare("t") == 0){
                    this->record_prop_times.push_back(num_steps);
                }
                else{
                    unsigned time = unsigned(stoul(tmp));
                    this->record_prop_times.push_back(time);
                }
                in >> tmp;
                std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

                sort(this->record_prop_times.begin(), this->record_prop_times.end());

                this->record_prop_times.erase(unique( this->record_prop_times.begin(), this->record_prop_times.end() ), this->record_prop_times.end() );
            }
            this->log_propagator = true;
        }
    }
}

void Parameters::readHexagonalParams(ifstream &in)
{
    string tmp="";

    while(true)
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(str_dist(tmp,"radius") <= 1){
            in >> hex_packing_radius;
        }
        if(str_dist(tmp,"separation") <= 1){
            in >> hex_packing_separation;
        }

        if(str_dist(tmp,"icvf") <= 1){
            in >> hex_packing_icvf;
        }

        if(str_dist(tmp,"</cylinder_hex_packing>") == 0 || str_dist(tmp,"</sphere_hex_packing>") == 0){
            break;
        }
    }
}


void Parameters::readGammaParams(ifstream &in)
{
    string tmp="";

    while(str_dist(tmp,"</cylinder_gamma_packing>") > 0 || str_dist(tmp,"</sphere_gamma_packing") > 0 )
    {
        in >> tmp;
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);
        if(str_dist(tmp,"output_conf") <= 1){
            string tst;

            in >> tst;
            in >> gamma_output_conf;
        }
        else if(str_dist(tmp,"alpha") <= 1 or str_dist(tmp,"shape") <= 1){
            in >> gamma_packing_alpha;
        }
        else if(str_dist(tmp,"beta") <= 1 or str_dist(tmp,"scale") <= 1){
            in >> gamma_packing_beta;
        }
        else if(str_dist(tmp,"num_cylinders") <= 1){
            in >> gamma_num_obstacles;
        }
        else if(str_dist(tmp,"num_spheres") <= 1){
            in >> gamma_num_obstacles;
        }
        else if(str_dist(tmp,"icvf") <= 1){
            in >> gamma_icvf;
        }
        else if(str_dist(tmp,"min_radius") <= 1){
            in >> min_obstacle_radii;
        }
        else if(str_dist(tmp,"") == 0){
            in.clear();
            //in.ignore();
        }
        else if(str_dist(tmp,"</cylinder_gamma_packing>") ==0 || str_dist(tmp,"</sphere_gamma_packing>") == 0  ){
            break;
        }

        tmp = "";
    }
}

void Parameters::readSubdivisionFile()
{
    ifstream in(subdivisions_file);

    Subdivision tmp;

    while( in >> tmp.min_limits[0]){
        in >> tmp.min_limits[1];
        in >> tmp.min_limits[2];
        in >> tmp.max_limits[0];
        in >> tmp.max_limits[1];
        in >> tmp.max_limits[2];
        subdivisions.push_back(tmp);
    }
}


void Parameters::addSubdivisions()
{

    if( (number_subdivisions > 0) && (voxels_list.size() <=0) ){
        SimErrno::error("subdivisions_number parameter passed without a defined voxel.",cout);
        assert(0);
    }
    float gap[3];
    for(uint i = 0; i < 3; i++){
        gap[i] =  float(this->voxels_list[0].second[i] - this->voxels_list[0].first[i])/float(number_subdivisions);
    }

    for(uint x_ind = 0 ; x_ind < number_subdivisions; x_ind++){
        for(uint y_ind = 0 ; y_ind < number_subdivisions; y_ind++){
            for(uint z_ind = 0; z_ind < number_subdivisions; z_ind++){

                Subdivision tmp;

                //x
                tmp.min_limits[0] = float(voxels_list[0].first[0]  + x_ind*gap[0]);
                //y
                tmp.min_limits[1] = float(voxels_list[0].first[1]  + y_ind*gap[1]);
                //z
                tmp.min_limits[2] = float(voxels_list[0].first[2]  + z_ind*gap[2]);

                //x
                tmp.max_limits[0] = tmp.min_limits[0] + gap[0];
                //y
                tmp.max_limits[1] = tmp.min_limits[1] + gap[1];
                //z
                tmp.max_limits[2] = tmp.min_limits[2] + gap[2];

                subdivisions.push_back(tmp);

                //cout << ' ' << tmp.min_limits[0] << ' ' << tmp.min_limits[1] <<  ' ' <<tmp.min_limits[2] << endl;
                //cout << ' ' << tmp.max_limits[0] << ' ' << tmp.max_limits[1] <<  ' ' <<tmp.max_limits[2] << endl;
            }
        }
    }

}

void Parameters::readPropagatorDirections(string dir_path)
{
    ifstream in(dir_path);

    if(in.fail()){

        cout << " WHY!!! " << endl;
        //File does not exist code here
    }

    Eigen::Vector3f direction;


    while( in >> direction[0]){
        in >> direction[1];
        in >> direction[2];

        if(direction.norm() >0)
            this->prop_dirs.push_back(direction.normalized());
    }

    in.close();
}

void Parameters::readPLYFileList(string path){

    ifstream in(path);

    if(in.fail()){
        SimErrno::error("PLY file list not found in:",cout);
        cout << path << endl;
        assert(0);
    }

    float scale;
    in >> scale;

    if(scale <=0.0){
        SimErrno::error("PLY scale must be a positive number",cout);
        assert(0);
    }

    if(scale >=1e6 || scale <= 1e-6){
        SimErrno::warning("PLY may be unsuitable for simulation.",cout);
        assert(0);
    }

    string ply_file;
    while( in >> ply_file){
        PLY_files.push_back(ply_file);
        PLY_scales.push_back(scale);
        PLY_percolation.push_back(0.0);
    }
    in.close();
}

void Parameters::readPLYFileListScalePercolation(string path)
{
    ifstream in(path);

    if(in.fail()){
        SimErrno::error("PLY file list not found in:",cout);
        cout << path << endl;
        assert(0);
    }

    float scale,percolation;
    string ply_file;
    while( in >> ply_file){
        in >> scale;
        in >> percolation;
        PLY_files.push_back(ply_file);
        PLY_scales.push_back(scale);
        PLY_percolation.push_back(percolation);
    }
    in.close();
}
