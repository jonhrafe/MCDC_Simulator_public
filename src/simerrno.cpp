#include "simerrno.h"
#include "iostream"
#include "constants.h"
#include <thread>
#include <iostream>
#include <sstream>
#include <string>
#include <queue>
#include<iomanip>
using namespace std;

SimErrno::SimErrno()
{

}

bool SimErrno::checkSimulationParameters(Parameters &params)
{
    cout << "\033[1;35m/***************   MC/DC Simulation parameters check out:  ***************/" << SH_DEFAULT << "\n";

    if(params.num_walkers > 1e9){
        error( " Maximum number of particles is fixed to 1e9.",cout);
        assert(0);
        return true;
    }

    if(params.num_walkers < 1){
        error( " Minimum number of particles is fixed to 1." ,cout);
        assert(0);
        return true;
    }

    if(params.num_steps > 1e7){
        error( " Maximum number of steps is fixed to 1e7.",cout);
        assert(0);
        return true;
    }

    if(params.num_steps < 1){
        error( " Minimum number of steps is fixed to 1.",cout);
        assert(0);
        return true;
    }

    if(params.sim_duration > 200.0 && params.scale_from_stu == 1){
        warning("Simulation duration might be unsuitable.", cout);
    }

    if(params.sim_duration > 1 && params.scale_from_stu == 0){
        warning("Simulation duration might be unsuitable.",cout);
    }

    if(params.sim_duration <= 0.0){
        error( " Simulation duration wrongly initialized.",cout);
        assert(0);
        return true;
    }

    if(params.diffusivity <= 0.0){
        error( " Paticle diffusivity wrongly initialized.",cout);
        assert(0);
        return true;
    }

    if (params.num_proc == 0){
         warning( "Number of processors is set by defult (1).",cout);
        params.num_proc = 1;
    }

    unsigned int nthreads = std::thread::hardware_concurrency();
    if (params.num_proc > nthreads){
         warning( "The number of processors to be used (" + to_string(params.num_proc) + ") is higher than the physical available processors (" + to_string(nthreads) + ").",cout);
    }

    if(params.scheme_file.size() >1){
        assert(checkSchemeFile(params));
    }

    if(params.PLY_files.size() > 0){
        info("Checking PLY format...",cout);
        assert(checkPLYFiles(params));
        info("Done...",cout);
    }

    if(params.cylinders_files.size()>0){
        info("Checking Cylinder list format...",cout);
        checkCylindersListFile(params);
        info("Done...",cout);
    }
    if(params.cylinders_files.size()>0){
        info("Checking Sphere list format...",cout);
        checkCylindersListFile(params);
        info("Done...",cout);
    }

    if(params.ini_walkers_file.size() > 2){
        info("Checking walker initial position list format...",cout);
        checkInitWalkerFile(params);
        info("Done...",cout);
    }

    if(params.voxels_list.size()>1){
        checkVoxelLimits(params);
    }

    if(params.hex_cyl_packing == true){
        if(params.hex_packing_radius<= 0){
            error( "Cylinder radius incoherent: " + to_string(params.hex_packing_radius) ,cout);
            assert(0);
            return true;
        }

        if(params.hex_packing_icvf > 0.90){
            error( "Max achievable ICVF is 0.9 ",cout);
            assert(0);
            return true;
        }

        if(params.hex_packing_icvf <= 0.0){
            error( "ICVF must be greater than 0.0 ",cout);
            assert(0);
            return true;
        }else{
            params.hex_packing_separation = sqrt( (2*M_PI*params.hex_packing_radius*params.hex_packing_radius)/(sqrt(3)*params.hex_packing_icvf));
        }


        if(params.hex_packing_separation - 2.0*params.hex_packing_radius < 0.0){
            error( "Cylinder separation can't be less that twice the radius (or epsilon close): " + to_string(params.hex_packing_separation) ,cout);
            assert(0);
            return true;
        }

        if(params.hex_packing_separation - 2.0*params.hex_packing_radius <= 1e-6){
            warning("Cylinder separation is too close (barrier collision): " + to_string(params.hex_packing_separation) ,cout);
        }

    }

    if(params.hex_sphere_packing == true){

        if(params.hex_packing_radius<= 0){
            error( "Spheres' radius incoherent: " + to_string(params.hex_packing_radius) ,cout);
            assert(0);
            return true;
        }

        if(params.hex_packing_icvf > 0.69){
            error( "Max achievable ICVF is 0.69 ",cout);
            assert(0);
            return true;
        }

        if(params.hex_packing_icvf <= 0.0){
            error( "ICVF must be greater than 0.0 ",cout);
            assert(0);
            return true;
        }else{
            params.hex_packing_separation = pow((4.*4./3.*M_PI*params.hex_packing_radius*params.hex_packing_radius*params.hex_packing_radius)/(3*params.hex_packing_icvf),1./3.);
        }


        if(params.hex_packing_separation - 2.0*params.hex_packing_radius < 0.0){
            error( "Cylinder separation can't be less that twice the radius (or epsilon close): " + to_string(params.hex_packing_separation) ,cout);
            assert(0);
            return true;
        }

        if(params.hex_packing_separation - 2.0*params.hex_packing_radius <= 1e-6){
            warning("Cylinder separation is too close (barrier collision): " + to_string(params.hex_packing_separation) ,cout);
        }

    }

    if(params.gamma_cyl_packing){
        checkGammaDistributionParamaters(params);
    }

    if(params.subdivision_flag){
        if(params.subdivisions.size() > 1000 ){
            warning("Huge number of sudivision voxels. A considerable amount of RAM will be needed for the ouput computation.",cout);
        }

        if(params.number_subdivisions > 100){
            error("Unrealistic number of resulting subdivision voxels : " + std::to_string(params.number_subdivisions) + "^3",cout);
            assert(0);
            return true;
        }

        if( (params.number_subdivisions > 0) && (params.voxels_list.size() <=0) && params.gamma_cyl_packing ==false){
            error("subdivisions_number parameter passed without a defined voxel.",cout);
            assert(0);
            return true;
        }

        if(params.subdivisions_file.size() > 2)
            checkSubdivisionsFile(params);
    }

    checkOuputPrefixAndWriteInfo(params);

    if(params.obstacle_permeability < 0.0 || params.obstacle_permeability > 1){
        error(" Permeability coefficient must be set in the range [0,1].",cout);
        assert(0);
        return true;
    }

    if(!(params.write_bin || params.write_txt)){
        error(" No output will be written; write_bin and write_txt flags are deactivated.",cout);
        assert(0);
        return true;

    }

    if(params.custom_sampling_area){

        for(auto j = 0; j < 3; j++){
            for (unsigned i=0; i < params.voxels_list.size();i++)
                if((params.voxels_list[i].first[j]  - params.min_sampling_area[j])>1e-8 || (params.max_sampling_area[j]-params.voxels_list[i].second[j])>1e-8)
                {
                    SimErrno::error("Custom sampling area cannot be outside the defined voxel\n",cout);
                    assert(0);
                }

            if(params.max_sampling_area[j]-params.min_sampling_area[j] <=  0 ){
                SimErrno::error("Custom sampling area wrongly defined (bad limits)\n",cout);
                assert(0);
            }
        }
    }

    if(params.computeVolume && params.voxels_list.size() <=0 && params.gamma_cyl_packing==false and params.hex_cyl_packing ==false and params.hex_sphere_packing ==false and params.gamma_sph_packing ==false){
        warning(" Flag: 'compute_volume' ignored, no voxel."  ,cout);
    }

    return false;
}

bool SimErrno::checkSchemeFile(Parameters &params)
{
    info("Checking Sequence Scheme file format...",cout);

    ifstream in(params.scheme_file.c_str());

    if(!in.is_open()){
        error( " Scheme file cannot be open: " + params.scheme_file ,cout);
        in.close();
        return false;
    }

    string header;
    in >> header;
    in >> header;

    std::size_t found = header.find("STEJSKALTANNER");
    std::size_t found_APGSE = header.find("APGSE");
    std::size_t found_waveForm = header.find("WAVEFORM");
    if(found!=std::string::npos){
        vector<double> sample_vector;

        unsigned counter = 0;
        double tmp = 0;
        while( in >> tmp){
            if (counter < 7)
                sample_vector.push_back(tmp);
            counter++;
        }

        if(params.scale_from_stu == 1){
            if(sample_vector[6] > 1.0 || sample_vector[3] > 1.0){
                warning("Scheme file might not be in standard units (meters, seconds, Tesla). Units Warning.", cout);
            }
        }
        else{
            if(sample_vector[6] < 1.0 || sample_vector[3] < 1.0){
                warning("Scheme file might be in standard units (meters, seconds, Tesla). Units Warning.",cout);
            }
        }

        if(counter%7 != 0){
            error("Scheme file has inconsistent format. PGSE Format ERROR.",cout);
            assert(0);
        }
    }
    else if(found_APGSE!=std::string::npos){
        vector<double> sample_vector;

        unsigned counter = 0;
        double tmp = 0;
        while( in >> tmp){
            if (counter < 9)
                sample_vector.push_back(tmp);
            counter++;
        }

        if(params.scale_from_stu == 1){
            if(sample_vector[8] > 1.0 || sample_vector[3] > 1.0){
                warning("Scheme file might not be in standard units (meters, seconds, Tesla). Units Warning.", cout);
            }

        }
        else{
            if(sample_vector[8] < 1.0 || sample_vector[3] < 1.0){
                warning("Scheme file might be in standard units (meters, seconds, Tesla). Units Warning.",cout);
            }
        }

        if(counter%9 != 0){
            error("Scheme file has inconsistent format. APGSE Format ERROR.",cout);
            assert(0);
        }

    }
    else if (found_waveForm!=std::string::npos){

        float wave_duration,wave_bins,num_rep;

        in >> wave_duration;
        in >> wave_bins;
        float holder;
        in >> holder;
        num_rep = uint(holder);

        if(params.scale_from_stu == 1){
            if(wave_duration > 1){
                warning("Scheme file might not be in standard units (meters, seconds, Tesla). Units Warning.", cout);
            }
        }
        else{
            if(wave_duration <= 1){
                warning("Scheme file might be in standard units (meters, seconds, Tesla). Units Warning.",cout);
            }
        }

        if(params.scale_from_stu == 1)
            wave_duration *= s_to_ms;

        if(params.sim_duration > double(wave_duration)+EPS_VAL){
            warning("Gradient waveform is shorter than the dynamic duration. The Waveform will be filled with 0's.",cout);
        }
        else if(params.sim_duration < double(wave_duration)-EPS_VAL){
            warning("Gradient waveform TE is larger than the dynamic duration.",cout);
        }

        unsigned counter = 0;
        double tmp = 0;
        while( in >> tmp){
            counter++;
        }

        if(counter != uint(wave_bins*num_rep*3)){
            error("Waveform Scheme file has inconsistent size. WAVEFORM Format ERROR.",cout);
            assert(0);
        }
    }
    else{
        error( "Scheme file version error. Valid Header was not found. Check documentation." ,cout);
        return false;
    }

    info("Done...",cout);
    return true;
}

bool SimErrno::checkPLYFiles(Parameters &params)
{

    bool degenerated = false;
    unsigned int degenerated_triangles = 0;
    for (unsigned i = 0 ;i < params.PLY_files.size(); i++)
    {
        std::ifstream in(params.PLY_files[i].c_str(),std::ifstream::in);

        if(!in){
            error( " PLY file cannot be opened: " ,cout);
            in.close();
            return false;
        }

        std::string first_word;
        in >> first_word;
       if(first_word.compare("ply")){
           error( " Input file is not a PLY mesh model: Missing  \"ply\"  header" ,cout);
           in.close();
           return false;
       }

        unsigned vert_number=0,face_number=0;
        std::string tmp = "";

        while(tmp.compare("end_header")){
            in >> tmp;
            //cout << tmp << endl;

            if(!tmp.compare("vertex")){
                in >> vert_number;
            }
            if(!tmp.compare("face")){
                in >> face_number;
            }

            if(!tmp.compare("nx")){
                in.close();
                error( " PLY file should not contain face normals, or vertex colors. PLY format error: " ,cout);
                in.close();
                return false;
            }
        }


        Eigen::Matrix3Xf vertices(3,vert_number);

        // We load all the vertices in a strucutre of size (3,num_of_vertices)
        for(unsigned v = 0 ; v < vert_number; v++){
            vector<float> tmp_v = {0,0,0};
            in >> tmp_v[0];
            in >> tmp_v[1];
            in >> tmp_v[2];
            vertices(0,v) = tmp_v[0];
            vertices(1,v) = tmp_v[1];
            vertices(2,v) = tmp_v[2];
        }

        double num;
        //for each face (index f) we check the distance between the 3 edges
        for(unsigned f = 0 ; f < face_number; f++)
        {
            in >> num;
            // if something is not a triangle then we throw an error
            if(num != 3.0){
                in.close();
                error( " PLY mesh should be completely triangulated. PLY format error: ",cout);
                in.close();
                return false;
            }

            //checks if the triangles are "degenerated"
            unsigned tmp_e[3]; // the indexes of its vertices (3 numbers)
            float tmp_t[3][3]; // (the actual coordinates (triangle))
            in >> tmp_e[0];
            in >> tmp_e[1];
            in >> tmp_e[2];
            for(int ii = 0 ; ii < 3; ii++)
                for (int jj = 0 ; jj < 3; jj++)
                    tmp_t[ii][jj] = vertices(jj,tmp_e[ii]);

            float edge_lengths[3]= {0,0,0};
            for(int ii = 0 ; ii < 3; ii++){
                unsigned vertex_a = ii;
                unsigned vertex_b = (ii+1)%3;

                for (unsigned jj=0; jj < 3; jj++)
                    edge_lengths[ii] += (tmp_t[vertex_a][jj] - tmp_t[vertex_b][jj])*(tmp_t[vertex_a][jj] - tmp_t[vertex_b][jj]) ;
            }

            for(int ii = 0 ; ii < 3; ii++)
                for (unsigned jj=ii; jj < 3; jj++)
                    if(edge_lengths[jj]>0)
                        if(edge_lengths[ii]/edge_lengths[jj] > 1000 or edge_lengths[ii]/edge_lengths[jj] < 1e-3 ){
                            degenerated=true;
                            degenerated_triangles++;
                        }

        }
        in.close();
    }

    if(degenerated){
        warning( "PLY contains ("+  std::to_string(degenerated_triangles) + ") highly irregular triangles. Possible numerical errors and optimization failures may occur.",cout);

    }


    if(params.PLY_files.size() > params.PLY_scales.size()){
        warning( "PLY scale is not set for all files. Scale will be set as default (1e-3). Substrate scale Warning.",cout);

        while(params.PLY_files.size() > params.PLY_scales.size())
            params.PLY_scales.push_back(1.0e-3);

    }

    return true;
}


//* Auxiliare method to split words in a line using the spaces*//
template<typename Out>
void split__(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split_(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split__(s, delim, std::back_inserter(elems));
    return elems;
}


bool SimErrno::checkCylindersListFile(Parameters &params)
{
    for(unsigned i = 0; i < params.cylinders_files.size(); i++){
        bool z_flag = false;
        ifstream in(params.cylinders_files[i]);

        if(!in){
            error( "Cylinder list file cannot be open." ,cout);
            assert(0);
            in.close();
            return true;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {
                std::vector<std::string> jkr = split_(line,' ');
                if (jkr.size()!= 1){
                    error( "First line must be only the overall scale factor: ",cout);
                    in.close();
                    assert(0);
                    return true;
                }
                first-=1;continue;
            }

            std::vector<std::string> jkr = split_(line,' ');

            if(jkr.size() != 7 && jkr.size() != 4){
                error( "Cylinder list file is not in the correct format." ,cout);
                in.close();
                assert(0);
                return true;
            }

            if (jkr.size() != 7){
                z_flag = true;
                warning("No cylinders orientation inlcluded. Cylinder orientation was set towards the Z direction by default for all cylinders.",cout);
            }
            break;
        }
        in.close();

        in.open(params.cylinders_files[i]);

        if(!z_flag){
            double x,y,z,ox,oy,oz,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> ox >> oy >> oz >> r)
            {
                if ((x - ox) == 0.0 && (z - oz) == 0.0 && (y - oy) == 0.0){
                    error( "Cylinder list has wrongly defined cylinders. Invalid orientation: ",cout);
                    in.close();
                    assert(0);
                    return true;
                }
            }
            in.close();
        }
    }

    return true;
}


bool SimErrno::checkSphereListFile(Parameters &params)
{
    for(unsigned i = 0; i < params.spheres_files.size(); i++){
        ifstream in(params.cylinders_files[i]);

        if(!in){
            error( "Spheres list file cannot be open." ,cout);
            assert(0);
            in.close();
            return true;
        }

        bool first=true;
        for( std::string line; getline( in, line ); )
        {
            if(first) {
                std::vector<std::string> jkr = split_(line,' ');
                if (jkr.size()!= 1){
                    error( "First line must be only the overall scale factor: ",cout);
                    in.close();
                    assert(0);
                    return true;
                }
                first-=1;continue;
            }

            std::vector<std::string> jkr = split_(line,' ');

            if(jkr.size() != 4){
                error( "Sphere list file is not in the correct format." ,cout);
                in.close();
                assert(0);
                return true;
            }
        }
        in.close();
    }
    return true;
}

bool SimErrno::checkInitWalkerFile(Parameters &params)
{

    bool warning_ = false;

    ifstream in(params.ini_walkers_file);

    if(!in){
        error( "Walkers initial positions file cannot be open.",cout);
        in.close();
        assert(0);
        return true;
    }

    // check if all particles are inside the voxel.
    unsigned long count = 0;
    double temp;
    int limit_index=0;
    while(in >> temp){
        count++;
        if(params.voxels_list.size()>0){
            if(temp <= params.voxels_list[0].first[limit_index] || temp >= params.voxels_list[0].second[limit_index]){
                error( "At least one walker initial position is not inside the defined voxel. "+ to_string(temp) ,cout);
                assert(0);
            }
        }
        limit_index = (limit_index<2)?limit_index+1:0;
    }

    params.ini_walkers_file_count = uint(count/3.0);

    if( (count/3.0) < 1.0 ){
        error( "No initial positions found in the walkers positions file." ,cout);
        assert(0);
        return true;
    }

    if(count %3 != 0){
        error( "List of initial position should include x,y,z position for all walkers, check positions format: " ,cout);
        assert(0);
        return true;
    }

    if(count/3.0 <= params.num_walkers){
        warning("Positions file has less positions than initialized walkers. Positions will be repeated.",cout);
        warning_ = true;
    }


    in.close();

    return warning_;

}

bool SimErrno::checkVoxelLimits(Parameters &params)
{
    if(params.voxels_list.size()>1){
        error( "Only single voxel simulations are supported.",cout);
        assert(0);
        return true;
    }

    return false;

}

bool SimErrno::checkConfigurationFile(const char* configuration_file)
{
    string scheme_ = configuration_file;
    std::string scheme_filename = scheme_.substr(scheme_.find_last_of("/\\") + 1);
    cout << SH_FG_PURPLE<< '/' << setfill('*') << setw(49) << "  " +scheme_filename+ " " << setw(30) << "/\n" SH_DEFAULT ;

    info("Checking configuration file labels...",cout);

    int count_tag_obstacle=0,count_tag_voxels=0,count_tag_log=0,count_tag_delta=0;
    int count_tag_phase = 0, count_tag_positions = 0, count_hexa_obstacle_tag=0;
    int count_tag_sampling_area=0;


    ifstream in(configuration_file);

    if(!in){
        error( "Cannot open the configuration file: " + string(configuration_file) ,cout);
        assert(0);
        return true;
    }


    string tmp="";
    bool ended = false;

    bool voxel_defined = false;
    bool fixed_configuration = false;

    while(in >> tmp ){

        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        if(Parameters::str_dist(tmp,"<voxels>") <= 0){
            count_tag_voxels++;
            voxel_defined=true;
        }
        else if(Parameters::str_dist(tmp,"</voxels>") <= 0){
            count_tag_voxels--;
        }
        else if(Parameters::str_dist(tmp,"<delta>") <= 0){
            count_tag_delta++;
        }
        else if(Parameters::str_dist(tmp,"</delta>") <= 0){
            count_tag_delta--;
        }
        else if(Parameters::str_dist(tmp,"<obstacle>") <= 0){
            count_tag_obstacle++;
        }
        else if(Parameters::str_dist(tmp,"</obstacle>") <= 0){
            count_tag_obstacle--;
        }
        else if(Parameters::str_dist(tmp,"<log>") <= 0){
            count_tag_log++;
        }
        else if(Parameters::str_dist(tmp,"<cylinder_hex_packing>") <= 0){
            count_hexa_obstacle_tag++;
            fixed_configuration = true;
        }
        else if(Parameters::str_dist(tmp,"<sphere_hex_packing>") <= 0){
            count_hexa_obstacle_tag++;
            fixed_configuration = true;
        }
        else if(Parameters::str_dist(tmp,"</log>") <= 0){
            count_tag_log--;
        }
        else if(Parameters::str_dist(tmp,"<positions>") <= 0 ){
            count_tag_positions++;
        }
        else if(Parameters::str_dist(tmp,"</positions>") <= 0 ){
            count_tag_positions--;
        }
        else if(Parameters::str_dist(tmp,"<phase>") <= 0 ){
            count_tag_phase++;
        }
        else if(Parameters::str_dist(tmp,"</phase>") <= 0 ){
            count_tag_phase--;
        }
        else if(Parameters::str_dist(tmp,"</cylinder_hex_packing>") <= 0){
            count_hexa_obstacle_tag--;
        }
        else if(Parameters::str_dist(tmp,"</sphere_hex_packing>") <= 0){
            count_hexa_obstacle_tag--;
        }
        else if(Parameters::str_dist(tmp,"<spawning_area>") == 0){
            count_tag_sampling_area++;
        }
        else if(Parameters::str_dist(tmp,"</spawning_area>") == 0){
            count_tag_sampling_area--;
        }

        else if(Parameters::str_dist(tmp,"<end>") == 0){
            ended = true;
        }

    }

    if(count_tag_delta!= 0 ){
        error( "<delta> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }

    if(count_tag_voxels!= 0 ){
        error( "<voxels> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }

    if(count_tag_log!= 0 ){
        error( "<log> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }
    if(count_tag_obstacle!= 0 ){
        error( "<obstacle> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }
    if(count_tag_phase!= 0 ){
        error( "<phase> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }
    if(count_tag_positions!= 0 ){
        error( "<positions> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }
    if(count_hexa_obstacle_tag!= 0 ){
        error( "<obstacle_hex_packing> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }
    if(count_tag_sampling_area!= 0 ){
        error( "<spawning_area> tag is not properly set in: " + string(configuration_file),cout);
        assert(0);
        return true;
    }

    if(ended == false){
        warning("Configuration file did not end with the END tag, possible misbehaviour.",cout);
    }

    if(voxel_defined && fixed_configuration){
        warning("The defined voxel will be overridden by the substrate packing configuration",cout);
    }

    info("Done...",cout);

    return false;
}

void SimErrno::printSimulatinInfo(Parameters &params, ostream &out,bool color)
{
    string answer;
    if (color)
        out<< setfill('-')<< "\033[1;35m/********************   MC/DC Simulation Info:   *************************/" << SH_DEFAULT << "\n";
    else
        out<< setfill('-')<< "/***********************   MC/DC Simulation Info:   *************************/"  << "\n";

    infoMenu(" Software Version:      ------",  VERSION_ID,out, color,35);

    infoMenu(" Number of particles:   ------",  to_string(params.num_walkers),out, color,35);
    //out << SH_FG_GRAY << "[INFO]   " << SH_DEFAULT << " Number of particles:   ------" << setw(35) << params.num_walkers << endl;

    infoMenu(" Number of steps:       ------",  to_string(params.num_steps), out, color,35);

    infoMenu(" Number of cores:       ------",  to_string(params.num_proc ), out, color,35);

    if(params.scale_from_stu)
        infoMenu(" Diffusivity:           ------",  to_string(params.diffusivity*1e6)+"e-9 m^2/s",out, color,35);
    else
        infoMenu(" Diffusivity:           ------",  to_string(params.diffusivity*1e6)+"e-6 mm^2/ms",out, color,35);

    infoMenu(" Particle dynamics duration: -",  " " + to_string(params.sim_duration) +" ms" , out, color,35);

    answer = (params.PLY_files.size() > 0)?" true":" false";
    infoMenu(" PLY obstacles:         ------", answer, out, color,35);

    if(params.PLY_files.size() > 0)
        infoMenu(" Number of PLYs:        ------", to_string( params.PLY_files.size()),out, color,35);

    answer = (params.cylinders_files.size() > 0) || params.gamma_cyl_packing || params.hex_cyl_packing ?" true":" false";
    infoMenu(" Cylinder obstacles:    ------",  answer, out, color,35);

    if(params.hex_cyl_packing){
        infoMenu(" Hexagonal Configuration:  ---", "true", out, color,35);
        infoMenu(" Hex. radius:           ------",  " "+ to_string(params.hex_packing_radius*1e3)+" um",out, color,35);
        infoMenu(" Hex. ICVF:             ------",  " "+ to_string(params.hex_packing_icvf),out, color,35);
        //infoMenu(" Separation:            ------",  " "+ to_string(params.hex_packing_separation*1e3)+" um",out, color,35);
    }

    answer = (params.spheres_files.size() > 0) || params.gamma_sph_packing || params.hex_sphere_packing ?" true":" false";
    infoMenu(" Spherical obstacles:    ------",  answer, out, color,34);

    if(params.hex_sphere_packing){
        infoMenu(" Hex. sph. configuration:  ---", "true", out, color,35);
        infoMenu(" Hex. radius:           ------",  " "+ to_string(params.hex_packing_radius*1e3)+" um",out, color,35);
        infoMenu(" Hex. ICVF:             ------",  " "+ to_string(params.hex_packing_icvf),out, color,35);
        //infoMenu(" Separation:            ------",  " "+ to_string(params.hex_packing_separation*1e3)+" um",out, color,35);
    }

    if(params.gamma_cyl_packing){
        infoMenu(" Gamma Configuration:   ------", " true", out, color,35);
        infoMenu(" Gamma alpha:           ------",  " "+ to_string(params.gamma_packing_alpha)+" um",out, color,35);
        infoMenu(" Gamma scale:           ------",  " "+ to_string(params.gamma_packing_beta),out, color,35);
        infoMenu(" Target ICVF:           ------",  " "+ to_string(params.gamma_icvf),out, color,35);
        infoMenu(" Min. radius:           ------",  " "+ to_string(params.min_obstacle_radii)+" um",out, color,35);
    }
    if(params.gamma_sph_packing){
        infoMenu(" Gamma sph configuration: ----", " true", out, color,35);
        infoMenu(" Gamma alpha:           ------",  " "+ to_string(params.gamma_packing_alpha)+" um",out, color,35);
        infoMenu(" Gamma scale:           ------",  " "+ to_string(params.gamma_packing_beta),out, color,35);
        infoMenu(" Target ICVF:           ------",  " "+ to_string(params.gamma_icvf),out, color,35);
        infoMenu(" Min. radius:           ------",  " "+ to_string(params.min_obstacle_radii)+" um",out, color,35);
    }

    answer = (params.write_traj)?" true":" false";
    infoMenu(" Write trajfile:        ------",  answer, out, color,35);

    answer = (params.write_bin)?" true":" false";
    infoMenu(" Write to binary:       ------",  answer, out, color,35);

    answer = (params.write_txt)?" true":" false";
    infoMenu(" Write to txt:          ------",  answer, out, color,35);

    answer = (params.separate_signals)?" true":" false";
    infoMenu(" Separated signals      ------",  answer, out, color,35);

    answer = (params.scale_from_stu)?" true":" false";
    infoMenu(" Standard units:        ------",  answer, out, color,35);

    answer = (params.obstacle_permeability > 0)?" true":" false";
    infoMenu(" Permeability:          ------",  answer, out, color,35);

    if(params.obstacle_permeability > 0){
        infoMenu(" Permeability coeff:    ------", " " + std::to_string(params.obstacle_permeability), out, color,35);
    }

    answer = (params.seed != -1)?" true":" false";
    infoMenu(" Custom seed:           ------",  answer, out, color,35);

    if(params.seed != -1)
        infoMenu(" User's Seed:           -----", " " + to_string( params.seed), out, color,35);

    answer = (params.log_phase_shift)?" true":" false";
    infoMenu(" Write phase shift histogram: ",  answer, out, color,35);

    answer = (params.log_propagator)?" true":" false";
    infoMenu(" Write propagator file: ------",  answer, out, color,35);
    if(params.log_propagator){
        infoMenu(" Number of logged times: -----",  to_string(params.record_prop_times.size()), out, color,35);
    }
    answer = (params.ini_walkers_file.size() > 2)?" true":" false";
    infoMenu(" Walker initial position file: ", answer, out, color,34);

    answer = (params.record_pos_times.size() > 0)?" true":" false";
    infoMenu(" Save fixed walker positions: ",  answer, out, color,35);

    answer = (params.ini_delta_pos.size() > 0)?" true":" false";
    infoMenu(" Initial delta position: -----",  answer, out, color,35);

    if((params.ini_walker_flag.compare("intra")==0) || (params.ini_walker_flag.compare("extra")==0))
    infoMenu(" Walkers initial position: -----", " "+params.ini_walker_flag, out, color,33);


    infoMenu(" Number of voxels:      ------", " " + to_string( params.voxels_list.size()),out, color,35);

    if(params.custom_sampling_area)
        infoMenu(" Custom spawning area:  ------"," true" ,out, color,35);

    if(params.computeVolume == true)
        infoMenu(" Volume approximation:  ------"," true" ,out, color,35);

    if(params.subdivision_flag)
        infoMenu(" Number of subdivisions: -----", " " + to_string( params.subdivisions.size()),out, color,35);

    if(params.collision_sphere_distance > 0 ){
        infoMenu(" Collision sphere distance: --", " " + to_string( params.collision_sphere_distance),out, color,35);
    }

    answer = (params.discard_illegals == true)?" On":" Off";
    infoMenu(" Border Patrol          ------",  answer, out, color,35);

    answer = (params.discard_stucks == true)?" On":" Off";
    infoMenu(" Discard stuck spins    ------",  answer, out, color,35);


    if(params.max_simulation_time > 1){
        infoMenu(" Max simulation time:   --------", " " + to_string( params.max_simulation_time) +" secs",out, color,35);
    }

    if(params.scheme_file.length() > 1){
        std::string scheme_filename = params.scheme_file.substr(params.scheme_file.find_last_of("/\\") + 1);
        infoMenu(" Scheme file name:      ------", " "+scheme_filename,out,color,35);
    }

    infoMenu(" Date and Time:         ------", " " + currentDateTime() ,out,color,35);

    out << endl;

    if(params.write_txt){
        if(params.record_pos_times.size()){
            info("Text trajectory will be written for the following time steps: ",out,color);
            for (unsigned i = 0; i < params.record_pos_times.size(); i++){
                out << params.record_pos_times[i];

                if( i ==  params.record_pos_times.size()-1)
                    out << ".";
                else
                    out << ", ";

                if( !((i+1)%10) )
                    out << endl;

                if( i==10)
                    out << endl;
            }
            out << endl;
        }
    }

    if(params.write_traj){
        if(params.record_pos_times.size()){
            info("Binary trajectory will be written for the following time steps: ",out,color);
            for (unsigned i = 0; i < params.record_pos_times.size(); i++){
                out << params.record_pos_times[i];

                if( i< params.record_phase_times.size()-2)
                    out << ", ";
                else
                    out << ".";
            }
            out << endl;
        }
    }


    if(params.ini_delta_pos.size() > 0){
        info("Walkers to be initialize at: (" + to_string(params.ini_delta_pos[0]) + ", " + to_string(params.ini_delta_pos[1]) + ", " + to_string(params.ini_delta_pos[2]) + ")",out,color);
    }
}

void SimErrno::checkOuputPrefixAndWriteInfo(Parameters &params)
{
    info("Checking Ouput format...",cout);

    if(checkFileExist(params.output_base_name+"_simulation_info.txt")){
        appendRepetitionLabel(params);
    }

    ofstream out(params.output_base_name+"_simulation_info.txt");

    if(!out){
        error( "Cannot open write output files with the output <prefix> and location: "+ params.output_base_name,cout);
        assert(0);
        return;
    }

    printSimulatinInfo(params,out,false);

    out.close();
    info("Done...",cout);
}

bool SimErrno::checkGammaDistributionParamaters(Parameters &params)
{
    if( params.gamma_packing_alpha <= 0 || params.gamma_packing_beta <=0){
        error("Gamma distribution parameters are not on the supported range",cout,true);
        assert(0);
    }

    if( params.gamma_packing_alpha >= 20 ){
        warning("Shape parameter might be on a unsuitable range",cout,true);
    }

    if(params.gamma_packing_beta >=  10){
        warning("The inverse shape parameter might be on a unsuitable range",cout,true);
    }

    if(params.gamma_icvf <= 0 || params.gamma_icvf  > 1){
        error("ICVF should be in the range: (0,1] ",cout,true);
        assert(0);
    }

    if(params.voxels_list.size() > 0 ){
        warning("The voxel size will be overwritten by the gamma constructor",cout,true);
    }


    if(params.gamma_num_obstacles >= 1e6){
        warning("Number of cylinders to sample might be erroneous",cout,true);
    }

    return false;
}

void SimErrno::warning(string message, ostream &out, bool color)
{
    if(color)
        out <<  SH_BG_LIGHT_YELLOW <<  "[Warning]" << SH_DEFAULT << " " + message << endl;
    else
        out <<  "[Warning]" << " " + message << endl;

}

void SimErrno::info(string message, ostream &out, bool color)
{
        if(color)
            out <<  SH_FG_GREEN  <<  "[INFO]" << SH_DEFAULT << "    "+ message  << endl;
        else
            out <<  "[INFO]" << "     " + message  << endl;
}

void SimErrno::infoMenu(string message, string value,ostream &out, bool color,int space)
{
    if(color)
        out << SH_FG_GREEN<< "[INFO]   " << SH_DEFAULT << message << setw(space) <<  value << endl;
    else
        out << "[INFO]   " << message << setw(space) << value << endl;
}

void SimErrno::error(string message, ostream &out, bool color)
{
    if(color)
        out <<  SH_BG_RED <<  "[ERROR]" << SH_DEFAULT << "  " + message << SH_DEFAULT << endl;
    else
        out << " " + message << endl;
}

void SimErrno::expectedTime(string completed, string time, ostream & out, bool color, string steps_second ,string endl_str)
{
    if(color){
        std::string message = string(SH_FG_GREEN) + "[INFO]"+ SH_DEFAULT + "    [Completed: " + completed +"%]"  +  " [ETA: " + time + "]" +
                " ( "+ steps_second + " steps/second)" + endl_str;
        out << message;
        out.flush();
    }
    else{
        std::string message = "[INFO]    [Completed: " + completed +"%] "  +  " [ETA: " + time + "]"+" ( "+ steps_second + " steps/second)" + endl_str;
        out << message;
        out.flush();
    }

}

// Get current date/time, format is YYYY-MM-DD (HH:mm:ss)
std::string SimErrno::currentDateTime() {
    time_t     now = time(nullptr);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);

    strftime(buf, sizeof(buf), "%d-%m-%Y (%X)", &tstruct);

    return buf;
}

bool SimErrno::checkSubdivisionsFile(Parameters &params)
{

    if(!checkFileExist(params.subdivisions_file)){
        error( "Subdivision file cannot be open." ,cout);
        assert(0);
        return false;
    }

    ifstream in(params.subdivisions_file);

    unsigned count_lines = 0;
    bool flag_weird_scale = false;
    double min_pos[3];
    double max_pos[3];

    // checks if the file has full pairs of positions and format
    for( std::string line; getline( in, line ); )
    {
        count_lines++;
        std::vector<std::string> jkr = split_(line,' ');
        if(jkr.size() != 3){
            error("Incorrect format in subdivision file: 3 positions expected per line in: " + params.subdivisions_file,cout);
                        in.close();
            assert(0);
        }
    }
    in.close();

    if( (count_lines%2) != 0){
        error("Incorrect format in subdivision file: An even number of 3d positions are expected in: " + params.subdivisions_file,cout);
                    in.close();
        assert(0);
    }

    in.open(params.subdivisions_file);

    while (in >> min_pos[0] >> min_pos[1] >> min_pos[2] >> max_pos[0] >> max_pos[1] >> max_pos[2])
    {
        if ((min_pos[0] > max_pos[0]) || (min_pos[1] > max_pos[1]) || (min_pos[2] > max_pos[2])){
            error( "Incorrect format in subdivision file: Negative voxel size:",cout);
            in.close();
            assert(0);
            return false;
        }
        if( ((max_pos[0] - min_pos[0]) > 1) || ((max_pos[1] - min_pos[1]) > 1) || ((max_pos[2] - min_pos[2]) > 1)){
            flag_weird_scale = true;
        }
    }
    in.close();

    if(flag_weird_scale == 1)
        warning( "Huge sudivision voxel size: Subdivisions scale may be erroneus.",cout);

    return false;

}

void SimErrno::appendRepetitionLabel(Parameters &params)
{
    int rep_number = 0;
    string rep_label = "";
    //finds the next index
    while(true){
        rep_label = std::to_string(rep_number++);

        if(rep_label.size()==1){
            rep_label = "0"+rep_label;
        }

        if(!checkFileExist(params.output_base_name+"_rep_"+ rep_label +"_simulation_info.txt"))
            break;
    }

    params.output_base_name+=+"_rep_"+ rep_label ;

}
