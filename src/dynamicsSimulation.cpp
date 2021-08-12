/** Main class of the dynamic Simulation.
 *  Everything related with the diffusion and collision of the particles is either implemented
 *  in this or used in this class.
 *
 * startSimulation method carries most of the computational burden.
 *
 * Project MC/DC Simulator
 * @author Jonathan
 * @version 1.42 stable
 */


#include "dynamicsSimulation.h"
#include "constants.h"
#include <algorithm>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <assert.h>
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include <fstream>
#include "constants.h"
#include "collisionsphere.h"
#include "simerrno.h"
#include "simulablesequence.h"

using namespace Eigen;
using namespace std;
using namespace sentinels;

/**
 * DynamicsSimulation implementation
 */
DynamicsSimulation::DynamicsSimulation() {
    plyObstacles_list = nullptr;
    spheres_list      = nullptr;
    cylinders_list    = nullptr;

    params.num_walkers = 1; //N
    params.num_steps   = 1; //T
    params.traj_file   = "";
    ini_pos_file       = "";
    completed          = 0;
    max_simulation_time = 0;

    params.sim_duration = 1; //secs

    num_simulated_walkers = 0;

    params.diffusivity = DIFF_CONST;
    step_lenght = sqrt(6.0*params.diffusivity*params.sim_duration/params.num_steps);
    params.write_traj = trajectory.write_traj = false;
    params.write_txt = trajectory.write_txt   = false;

    if(params.seed > 0){
        mt.seed(ulong(params.seed));
    }
    else{
        //Random seed
        std::random_device rd;
        mt.seed(rd());
    }

    print_expected_time = 1;
    icvf=0;
    intra_tries=0;
    total_tries=0;
    aux_walker_index = 0;
}

/**
 * @param configuration file
 */
DynamicsSimulation::DynamicsSimulation(std::string conf_file) {
    plyObstacles_list = nullptr;
    spheres_list      = nullptr;
    cylinders_list    = nullptr;

    completed = 0;
    readConfigurationFile(conf_file);

    trajectory.initTrajectory(params);

    if(params.seed > 0){
        mt.seed(ulong(params.seed));
    }
    else{
        //Random seed
        std::random_device rd;
        mt.seed(rd());
    }
    print_expected_time = 1;

    this->step(1)=0;
    icvf=0;
    intra_tries=0;
    total_tries=0;
    aux_walker_index = 0;
}

/**
 * @param Parameter instance
 */
DynamicsSimulation::DynamicsSimulation(Parameters& params_) {

    plyObstacles_list = nullptr;
    spheres_list      = nullptr;
    cylinders_list    = nullptr;

    params = params_;
    completed = 0;
    trajectory.initTrajectory(params);

    if(params.seed > 0){
        mt.seed(ulong(params.seed));
    }
    else{
        //Random seed
        std::random_device rd;
        mt.seed(rd());
    }

    print_expected_time  = 1;
    this->step(1)=0;
    icvf=0;
    intra_tries=0;
    total_tries=0;
    aux_walker_index = 0;
}

void DynamicsSimulation::initObstacleInformation(){

    if(params.collision_sphere_distance<= 0){
        params.collision_sphere_distance = inner_col_dist_factor;
    }

    //Cylinders list of index initialization
    for(unsigned i= 0 ; i < (*cylinders_list).size();i++){
        cylinders_deque.push_back(i);

        if(params.obstacle_permeability > 0.0){
            (*cylinders_list)[i].percolation = params.obstacle_permeability;
        }
    }

    walker.cylinders_collision_sphere.collision_list        = &cylinders_deque;
    walker.cylinders_collision_sphere.list_size             = unsigned(cylinders_deque.size());
    walker.cylinders_collision_sphere.big_sphere_list_end   = walker.cylinders_collision_sphere.list_size;

    //Spheres list of index initialization
    for(unsigned i= 0 ; i < spheres_list->size();i++){
        spheres_deque.push_back(i);

        if(params.obstacle_permeability > 0.0){
            (*spheres_list)[i].percolation = params.obstacle_permeability;
        }
    }

    walker.spheres_collision_sphere .collision_list        = &spheres_deque;
    walker.spheres_collision_sphere.list_size             = unsigned(spheres_deque.size());
    walker.spheres_collision_sphere.big_sphere_list_end   = walker.spheres_collision_sphere.list_size;

    // PLY index list initialization
    for(unsigned i= 0 ; i < (*plyObstacles_list).size();i++){
        std::vector<unsigned> jkr;
        for(unsigned t =0; t < (*plyObstacles_list)[i].face_number;t++){
            jkr.push_back(t);
        }

        if(params.obstacle_permeability > 0.0){
            (*plyObstacles_list)[i].percolation = params.obstacle_permeability;
        }

        ply_deque.push_back(jkr);
        walker.ply_collision_sphere.small_sphere_list_end.push_back(0);
        walker.ply_collision_sphere.big_sphere_list_end.push_back(unsigned(jkr.size()));
        walker.ply_collision_sphere.list_size++;
    }
    walker.ply_collision_sphere.collision_list = &ply_deque;
}

void DynamicsSimulation::updatePropagator(Eigen::Matrix3Xd& log_pos_r)
{
    for ( uint t =0; t <  params.record_prop_times.size(); t++){
        auto time = params.record_prop_times[t];

        Eigen::Vector3f displacement;

        displacement[0] = float(log_pos_r(0,0) - log_pos_r(0,time));
        displacement[1] = float(log_pos_r(1,0) - log_pos_r(1,time));
        displacement[2] = float(log_pos_r(2,0) - log_pos_r(2,time));

        for (uint i=0; i < params.prop_dirs.size();i++){

            Eigen::Vector3f direction = params.prop_dirs[i];

            direction.normalize();

            Eigen::Vector3f d_projection = direction.dot(displacement)*direction;

            propagator.propagator_log[t][i] +=  d_projection.squaredNorm();
        }
    }
}

void DynamicsSimulation::normalizePropagator(float num_samples)
{

    if(num_samples<0)
        return;

    for (uint i =0 ; i < this->propagator.num_times; i++){
        for (uint j = 0 ; j < propagator.num_dirs; j++){
            propagator.propagator_log[i][j]/=num_samples;
        }
    }
}

void DynamicsSimulation::computeICVF()
{
    icvf = float(intra_tries)/float(total_tries);
}

bool DynamicsSimulation::finalPositionCheck()
{
    int cyl_id,ply_id,sph_id;

    if( ((*plyObstacles_list).size()>0) and sentinela.deport_illegals and params.obstacle_permeability <=0){

        bool isIntra = isInIntra(this->walker.pos_v,cyl_id,ply_id,sph_id,0);

        //cout << endl << endl << isIntra << " " << this->walker.location << "  " << walker.initial_location << endl;

        if((isIntra and this->walker.initial_location == Walker::extra) or ((!isIntra and this->walker.initial_location == Walker::intra))){
//            cout << "Im working" << endl;

//            cout << (this->walker.initial_location == Walker::intra) <<  "Intra"  << endl;
//            cout << isIntra << endl;
            return true;
        }
    }
    return false;
}

void DynamicsSimulation::writePropagator(std::string path)
{
    if(params.write_bin){
        ofstream bout;
        string path_bin = path+".bfloat" ;
        bout.open(path_bin.c_str(), std::ofstream::binary);

        if(!bout){
            std::cout << "Cannot open " << path << std::endl;
            return;
        }

        for (uint i =0 ; i < this->propagator.num_times; i++){
            for (uint j = 0 ; j < propagator.num_dirs; j++){
                float p = propagator.propagator_log[i][j];
                bout.write(reinterpret_cast<char *>(&p), sizeof(p));
            }
        }
        bout.close();
    }

    if(params.write_txt){

        ofstream tout;
        tout.open(path, std::ios::out);

        if(!tout){
            //TODO: Error handling
            std::cout << "Cannot open " << path.c_str()<< std::endl;
            return;
        }

        for (uint i =0 ; i < this->propagator.num_times; i++){
            for (uint j = 0 ; j < propagator.num_dirs; j++){
                float p = propagator.propagator_log[i][j];
                tout << p << " ";
            }
            tout << endl;
        }

        tout.close();
    }
}




void DynamicsSimulation::initSimulation()
{

    // Initial step length = sqrt(6*D*dt/T)
    step_lenght = sqrt(6.0*(params.diffusivity*params.sim_duration)/double(params.num_steps));

    // Writes the header file and opens .traj file (if required)
    trajectory.initTrajWriter();

    // Initialize the walker trajectory log
    walker.setNumberOfSteps(params.num_steps);

    //time step = dt/T
    time_step = params.sim_duration/double(params.num_steps);
    time_dt =    0;
    last_time_dt=0;

    //to predict time
    time(&start);

    //Initialize reading file for the walker position (if passed as parameter)
    std::ios_base::iostate exceptionMask = iniPos.exceptions() | std::ios::failbit;
    iniPos.exceptions(exceptionMask);

    // if a ini_walkers_file was passed
    if(params.ini_walkers_file.length() > 3){
        try {
            iniPos.open(params.ini_walkers_file);
            ini_pos_file_ini_index = (uint(id)*params.num_walkers)%params.ini_walkers_file_count;

            //We move the file start until we reach the initial position (MULTICORE SUPPPORT)
            for(unsigned i=0; i < ini_pos_file_ini_index; i++){
                double x,y,z;
                iniPos >> x; iniPos >> y; iniPos >> z;
            }
        }
        catch (std::ios_base::failure& e) {
            std::cerr << "Sim: " << id << " " << e.what() << '\n' << "Error loading walker initial positions. \n";
        }
    }

    initObstacleInformation();

    //Flags for the crossing and stuck particles. (numerical error sentinels)
    sentinela.deport_illegals = params.discard_illegals;
    sentinela.discard_stucks  = params.discard_stucks;

    if(params.log_propagator){
        propagator.num_dirs = uint(params.prop_dirs.size());
        propagator.num_times = uint(params.record_prop_times.size());
        propagator.initPropagator();
    }


    if(params.custom_sampling_area == false and voxels_list.size()>0){
        for(auto i = 0; i<3;i++){
           params.min_sampling_area[i]=voxels_list[0].min_limits[i];
           params.max_sampling_area[i]=voxels_list[0].max_limits[i];
        }
    }
}


bool DynamicsSimulation::expectedTimeAndMaxTimeCheck(unsigned w)
{
    /******  Some code to get a predicted time  **********/
    double completed_perc = int((w+1.0)/(params.num_walkers+1.0)*100.0);
    //unsigned tent_num_walkers = params.num_walkers/10 + 1;

    time(&now);

    second_passed = now - start;

    if(print_expected_time){
        if (completed_perc >= completed){

            if(this->completed <= 0.0){
                SimErrno::expectedTime(to_string(int(completed)), "Unknown",cout,true,"?","");
                cout.flush();
                completed+=5;
            }
            else if(completed_perc >= completed){
                cout << string(50,' ');
                cout << string(200,'\b');
                cout.flush();
                int steps_per_second = int(float(w*params.num_steps)/second_passed)*params.num_proc;

                if(steps_per_second > 0){
                string s_p =  std::to_string(steps_per_second);

                SimErrno::expectedTime(to_string(int(completed)), secondsToMinutes( (double(params.num_walkers - w+1))*second_passed/double(w+1)),
                                       cout,true,s_p,"");

                }
                else{
                    SimErrno::expectedTime(to_string(int(completed)), "Unknown",cout,true,"?","");
                    cout.flush();
                }
                completed = max(completed+5.0,completed_perc);
                cout.flush();
            }
        }
        else if( w == params.num_walkers-1){
            cout << string(50,' ');
            cout << string(200,'\b');
            cout.flush();
            SimErrno::expectedTime("100", "0 seconds",cout,true,"\n");
            cout.flush();
        }
    }
    /****** END Some code to get a predicted time  **********/


    if( (params.max_simulation_time>0) && (second_passed >= params.max_simulation_time) ){
        return true;
    }

    return false;
}

void DynamicsSimulation::writeDWSignal(SimulableSequence* dataSynth)
{
    //Writes the output data
    if(dataSynth){
        dataSynth->writeResultingData(params.output_base_name);
        if(params.log_phase_shift)
            dataSynth->writePhaseShiftDistribution(params.output_base_name);
    }
}

void DynamicsSimulation::iniWalkerPosition()
{
    walker.initial_location = Walker::unknown;
    walker.location         = Walker::unknown;
    walker.intra_extra_consensus = walker.intra_coll_count = walker.extra_coll_count=walker.rejection_count=0;

/*
    if(params.custom_ini_walker_pos.size()>0){
        double x,y,z;
        x= params.custom_ini_walker_pos[aux_walker_index][0];
        y= params.custom_ini_walker_pos[aux_walker_index][1];
        z= params.custom_ini_walker_pos[aux_walker_index][2];
        aux_walker_index++;
        //cout << x << ' ' << y << ' ' << z << endl;
        walker.setInitialPosition(x,y,z);

        bool intra_flag =isInIntra(walker.ini_pos, walker.in_obj_index,walker.in_ply_index, 0.0);
        walker.location = (intra_flag==1)?Walker::RelativeLocation::intra:Walker::RelativeLocation::extra;
        walker.initial_location = walker.location;

    }
    //If the number of positions is less than the walkers, it restarts.
    else */if(iniPos.is_open()){
        double x,y,z;

        iniPos >> x; iniPos >> y; iniPos >> z;
        walker.setInitialPosition(x,y,z);

        if(++ini_pos_file_ini_index >= params.ini_walkers_file_count ){
            iniPos.clear();
            iniPos.seekg(0);
            ini_pos_file_ini_index = 0;
        }
    }
    else if (params.ini_delta_pos.size() > 0){
        walker.setRandomInitialPosition(Vector3d(double(params.ini_delta_pos[0]),double(params.ini_delta_pos[1]),double(params.ini_delta_pos[2])),
                Vector3d(double(params.ini_delta_pos[0]),double(params.ini_delta_pos[1]),double(params.ini_delta_pos[2])));
    }
    else if(params.ini_walker_flag.compare("intra")== 0){
        Vector3d intra_pos;
        getAnIntraCellularPosition(intra_pos,walker.in_obj_index,walker.in_ply_index,walker.in_sph_index);
        walker.setInitialPosition(intra_pos);
        walker.intra_extra_consensus--;
        walker.initial_location = Walker::intra;
    }
    else if(params.ini_walker_flag.compare("extra")== 0){
        Vector3d extra_pos;
        getAnExtraCellularPosition(extra_pos);
        walker.setInitialPosition(extra_pos);
        walker.initial_location = Walker::extra;
        walker.intra_extra_consensus++;
    }
    //Todo: poner esto bien sin el caso de hexapacking
    else if(voxels_list.size() > 0 or params.custom_sampling_area){
        walker.setRandomInitialPosition(params.min_sampling_area,params.max_sampling_area);
        if(params.computeVolume){
            bool intra_flag =isInIntra(walker.ini_pos, walker.in_obj_index,walker.in_ply_index, walker.in_sph_index, 0.0);
            walker.location = (intra_flag==1)?Walker::RelativeLocation::intra:Walker::RelativeLocation::extra;
            walker.initial_location = walker.location;
        }
    }
    else{
        walker.setInitialPosition(Vector3d(0,0,0));
    }
}


void DynamicsSimulation::initWalkerObstacleIndexes()
{

     //* Cylinders Collision Sphere *//

    // The outer collision sphere has a radius r = l*T 
    float outer_col_dist_factor = float(params.num_steps*step_lenght);

    walker.initial_sphere_pos_v = walker.pos_v;
    walker.cylinders_collision_sphere.setBigSphereSize(outer_col_dist_factor);
    
    // The inner collision sphere has radius l*T*collision_sphere_distance
    float inner_col_dist_factor = step_lenght*sqrt(params.num_steps)*params.collision_sphere_distance;
    walker.cylinders_collision_sphere.setSmallSphereSize(inner_col_dist_factor);

    // New version Cylinders obstacle selection
    walker.cylinders_collision_sphere.small_sphere_list_end = 0;
    walker.cylinders_collision_sphere.big_sphere_list_end = unsigned(cylinders_deque.size());

    // We add and remove the cylinder indexes that are or not inside sphere.
    for(unsigned i = 0 ; i < walker.cylinders_collision_sphere.list_size; i++ ){
        unsigned index = walker.cylinders_collision_sphere.collision_list->at(i);
        float dist = float((*cylinders_list)[index].minDistance(walker));
        if (dist < walker.cylinders_collision_sphere.small_sphere_distance){
            walker.cylinders_collision_sphere.pushToSmallSphere(i);
        }
    }


    //* Spheres Collision Sphere *//

    // The outer collision sphere has a radius r = l*T
    walker.spheres_collision_sphere.setBigSphereSize(outer_col_dist_factor);
    // The inner collision sphere has radius l*T*collision_sphere_distance
    walker.spheres_collision_sphere.setSmallSphereSize(inner_col_dist_factor);

    // New version  obstacle selection
    walker.spheres_collision_sphere.small_sphere_list_end = 0;
    walker.spheres_collision_sphere.big_sphere_list_end = unsigned(spheres_deque.size());

    // We add and remove the sphere indexes that are or not inside sphere.
    for(unsigned i = 0 ; i < walker.spheres_collision_sphere.list_size; i++ ){
        unsigned index = walker.spheres_collision_sphere.collision_list->at(i);
        float dist = float((*spheres_list)[index].minDistance(walker));
        if (dist < walker.spheres_collision_sphere.small_sphere_distance){
            walker.spheres_collision_sphere.pushToSmallSphere(i);
        }
    }

    //* PLY Collision Sphere *//
    
    walker.ply_collision_sphere.setBigSphereSize(outer_col_dist_factor);
    walker.ply_collision_sphere.setSmallSphereSize(inner_col_dist_factor);

    //cout << outer_col_dist_factor << endl;
    //cout << inner_col_dist_factor << endl;

    for(unsigned i = 0 ; i < walker.ply_collision_sphere.list_size; i++ )
    {
        walker.ply_collision_sphere.small_sphere_list_end[i] = 0;
        walker.ply_collision_sphere.big_sphere_list_end[i] = (*plyObstacles_list)[i].face_number;
        for(unsigned t = 0 ; t < (*plyObstacles_list)[i].face_number; t++){

            unsigned index = walker.ply_collision_sphere.collision_list->at(i)[t];
            float dist = float((*plyObstacles_list)[i].minDistance(walker,index));

            if (dist > walker.ply_collision_sphere.big_sphere_distance)
            {
                walker.ply_collision_sphere.popFromBigSphere(i,t);
            }

            if (dist < walker.ply_collision_sphere.small_sphere_distance)
            {
                walker.ply_collision_sphere.pushToSmallSphere(i,t);
            }
        }
    }


}


void DynamicsSimulation::updateCollitionSphere(unsigned t)
{
    float inner_ball_size = walker.ply_collision_sphere.small_sphere_distance;
    float outher_ball_size = walker.ply_collision_sphere.big_sphere_distance;

    float sphere_sqrd_displacement = float((walker.initial_sphere_pos_v-walker.pos_v).norm());

    if(sphere_sqrd_displacement  + float(step_lenght)   > outher_ball_size){
        initWalkerObstacleIndexes();
    }
    else if(sphere_sqrd_displacement + float(step_lenght) > inner_ball_size    ){
        updateWalkerObstacleIndexes(t);
    }
}

void DynamicsSimulation::getAnIntraCellularPosition(Vector3d &intra_pos,int &cyl_ind, int& ply_ind, int& sph_ind)
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);


    if(cylinders_list->size() <=0 and plyObstacles_list->size() <= 0 and spheres_list->size() <=0){
        SimErrno::error("Cannot initialize intra-axonal walkers within the given substrate.",cout);
        SimErrno::error("There's no defined intra-axonal compartment (missing obstacles?)",cout);
        assert(0);
    }
    if(voxels_list.size()<=0){
        SimErrno::error("Cannot initialize intra-cellular walkers within the given substrate, no voxel.",cout);
        assert(0);
    }

    unsigned count = 0;
    while(true){

        if(count > 100000){
            SimErrno::error("Cannot initialize intra-axonal walkers within the given substrate",cout);
            SimErrno::error("Max. number of tries to find an intra-celular compartment reached",cout);
            assert(0);
        }

        double x = double(udist(gen));
        double y = double(udist(gen));
        double z = double(udist(gen));

        x = x*(params.min_sampling_area[0]) + ( 1.0-x)*params.max_sampling_area[0];
        y = y*(params.min_sampling_area[1]) + ( 1.0-y)*params.max_sampling_area[1];
        z = z*(params.min_sampling_area[2]) + ( 1.0-z)*params.max_sampling_area[2];


       // cout << initialization_gap[2] << endl;
        Vector3d pos_temp = {x,y,z};

        if(checkIfPosInsideVoxel(pos_temp) && (isInIntra(pos_temp,cyl_ind,ply_ind, sph_ind,-0.1))){
            intra_pos = pos_temp;
            return;
        }
        count++;
    }
}

void DynamicsSimulation::getAnExtraCellularPosition(Vector3d &extra_pos)
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);
    int dummy_a,dummy_b,dummy_c;
    if(voxels_list.size()<=0){
        SimErrno::error("Cannot initialize extra-cellular walkers within the given substrate, no voxel.",cout);
        assert(0);
    }

    unsigned count = 0;

    while(true){

        if(count > 10000){
            SimErrno::error("Cannot initialize extra-cellular walkers within the given substrate",cout);
            SimErrno::error("Max. number of tries to find an extra-celular compartment reached",cout);
            assert(0);
        }

        double x = double(udist(gen));
        double y = double(udist(gen));
        double z = double(udist(gen));

        x = x*(params.min_sampling_area[0]) + ( 1.0-x)*params.max_sampling_area[0];
        y = y*(params.min_sampling_area[1]) + ( 1.0-y)*params.max_sampling_area[1];
        z = z*(params.min_sampling_area[2]) + ( 1.0-z)*params.max_sampling_area[2];

        Vector3d pos_temp = {x,y,z};

        if(checkIfPosInsideVoxel(pos_temp) && (!isInIntra(pos_temp, dummy_a,dummy_b,dummy_c,barrier_tickness))){
            extra_pos = pos_temp;
            return;
        }
        count++;
    }
}

bool DynamicsSimulation::checkIfPosInsideVoxel(Vector3d &pos)
{

    for(unsigned v = 0; v < voxels_list.size();v++){
        if(     pos[0] - voxels_list[v].min_limits[0] > barrier_tickness &&
                pos[1] - voxels_list[v].min_limits[1] > barrier_tickness &&
                pos[2] - voxels_list[v].min_limits[2] > barrier_tickness &&
                pos[0] - voxels_list[v].max_limits[0] < barrier_tickness &&
                pos[1] - voxels_list[v].max_limits[1] < barrier_tickness &&
                pos[2] - voxels_list[v].max_limits[2] < barrier_tickness)
            return true;
    }

    return false;
}

//TODO: Use t to decrease the size of the sphere.
void DynamicsSimulation::updateWalkerObstacleIndexes(unsigned t_)
{

    float outher_col_dist_factor = float(params.num_steps-t_+1.0*step_lenght);
    walker.ply_collision_sphere.setBigSphereSize(outher_col_dist_factor);

    walker.initial_sphere_pos_v = walker.pos_v;

    //Cylinders obstacle update.
    walker.cylinders_collision_sphere.small_sphere_list_end = 0;

    for(unsigned i = 0 ; i < walker.cylinders_collision_sphere.big_sphere_list_end; i++ )
    {
        unsigned index = walker.cylinders_collision_sphere.collision_list->at(i);
        float dist    = float((*cylinders_list)[index].minDistance(walker));

        if (dist > walker.cylinders_collision_sphere.big_sphere_distance)
        {
            walker.cylinders_collision_sphere.popFromBigSphere(i);
        }
        if (dist < walker.cylinders_collision_sphere.small_sphere_distance)
        {
            walker.cylinders_collision_sphere.pushToSmallSphere(i);
        }
    }

    //Spheres obstacle update.
    walker.spheres_collision_sphere.small_sphere_list_end = 0;

    for(unsigned i = 0 ; i < walker.spheres_collision_sphere.big_sphere_list_end; i++ )
    {
        unsigned index = walker.spheres_collision_sphere.collision_list->at(i);
        float dist    = float((*spheres_list)[index].minDistance(walker));

        if (dist > walker.spheres_collision_sphere.big_sphere_distance)
        {
            walker.spheres_collision_sphere.popFromBigSphere(i);
        }
        if (dist < walker.spheres_collision_sphere.small_sphere_distance)
        {
            walker.spheres_collision_sphere.pushToSmallSphere(i);
        }
    }

    //PLY update obstacle
    for(unsigned i = 0 ; i < walker.ply_collision_sphere.list_size; i++ )
    {
        walker.ply_collision_sphere.small_sphere_list_end[i] = 0;


        for(unsigned t = 0 ; t < walker.ply_collision_sphere.big_sphere_list_end[i]; t++){
            float dist  = INFINITY_VALUE;
            if((walker.in_ply_index <=0) || walker.in_ply_index == int(i)){
                unsigned triangle_index = walker.ply_collision_sphere.collision_list->at(i)[t];
                dist = float((*plyObstacles_list)[i].minDistance(walker,triangle_index));
            }

            if (dist > walker.ply_collision_sphere.big_sphere_distance)
            {
                walker.ply_collision_sphere.popFromBigSphere(i,t);
            }

            if (dist < walker.ply_collision_sphere.small_sphere_distance)
            {
                walker.ply_collision_sphere.pushToSmallSphere(i,t);
            }
        }
    }

    // cout  << " " << walker.collision_sphere_ply.small_sphere_list_end[0] << endl;
}

string DynamicsSimulation::secondsToMinutes(double t)
{
    if(t < 60){
        return std::to_string(int(t)) + " seconds";
    }

    int mins    = int(t/60.0);

    int seconds = int(t - mins*60);

    string min_s =  (mins>1)?" minutes":" minute";
    return std::to_string(mins) +  min_s +  " and " + std::to_string(seconds) + " seconds";

}


bool DynamicsSimulation::isInsideSpheres(Vector3d &position, int& sph_id,double distance_to_be_inside)
{
    Walker tmp;
    tmp.setInitialPosition(position);

    //track the number of positions checks for intra/extra positions

    for(unsigned i = 0 ; i < spheres_list->size(); i++){

        double dis = (*spheres_list)[i].minDistance(tmp);

        if( dis <= distance_to_be_inside ){
            intra_tries++;
            sph_id = i;
            return true;
        }
    }
    sph_id = -1;

    return false;
}

bool DynamicsSimulation::isInsideCylinders(Vector3d &position, int& cyl_id,double distance_to_be_inside)
{
    Walker tmp;
    tmp.setInitialPosition(position);

    //track the number of positions checks for intra/extra positions

    for(unsigned i = 0 ; i < cylinders_list->size(); i++){

        double dis = (*cylinders_list)[i].minDistance(tmp);

        if( dis <= distance_to_be_inside ){
            intra_tries++;
            cyl_id = i;
            return true;
        }
    }
    cyl_id = -1;

    return false;
}

bool DynamicsSimulation::isInsidePLY(Vector3d &position, int &ply_id,double distance_to_be_inside)
{
    ply_id= -1;

    //1) We find the closest PLY and triangle based on the triangle's center
    Walker tmp;
    tmp.setInitialPosition(position);

    double t,min_t = 1e6;
    unsigned min_j_index = 0;
    int min_i_index = -1;
    for (unsigned i=0; i < (*plyObstacles_list).size(); i++){
        for (unsigned j=0; j < (*plyObstacles_list)[i].face_number; j++){
            t = (position - (*plyObstacles_list)[i].faces[j].center).squaredNorm();
            // cout << t<< endl;
            if(t< min_t){
                min_i_index = i;
                min_j_index = j;
                min_t = t;
            }
        }
    }

    //2) We corroborate by casting an infinite ray and checking collisions

    Eigen::Vector3d ray = (-position + (*plyObstacles_list)[min_i_index].faces[min_j_index].center).normalized();
    Collision colision_temp;

    double new_min_t = 1e6;
    for (unsigned i=0; i < (*plyObstacles_list).size(); i++){
        for (unsigned j=0; j < (*plyObstacles_list)[i].face_number; j++){
            (*plyObstacles_list)[i].faces[j].stepIntersects_MT(tmp,ray,1e8,colision_temp);

            if(colision_temp.type == Collision::hit and new_min_t > colision_temp.t){
                new_min_t = colision_temp.t;
                min_i_index = i;
                min_j_index = j;
            }
        }
    }

    //3) Finally we check the sign of the closest collision. The sign indicates either intra or extra.
    if(min_i_index >= 0){
        Eigen::Vector3d normal;
        (*plyObstacles_list)[min_i_index].faces[min_j_index].getNormal(normal);
        //Orientation respect the triangle
        double dot = ((position - (*plyObstacles_list)[min_i_index].faces[min_j_index].center).normalized()).dot(normal);
        if (dot < distance_to_be_inside){
            intra_tries++;
            ply_id = min_i_index;
            return true;
        }
    }

    return false;
}


bool DynamicsSimulation::isInIntra(Vector3d &position, int& cyl_id,  int& ply_id, int& sph_id, double distance_to_be_intra_ply)
{
    bool isIntra = false;
    total_tries++;
    if(cylinders_list->size()>0){
        isIntra|= this->isInsideCylinders(position,cyl_id,barrier_tickness);
    }

    if(plyObstacles_list->size()>0){
        isIntra|=isInsidePLY(position,ply_id,distance_to_be_intra_ply);
    }

    if(spheres_list->size()>0){
        isIntra|=isInsideSpheres(position,sph_id,barrier_tickness);
    }
    return isIntra;
}




void DynamicsSimulation::startSimulation(SimulableSequence *dataSynth) {

    //Initialize values, arrays and files.
    initSimulation();

    //Alias of the step length, may vary when the time step is dynamic.
    double l = step_lenght;
    bool back_tracking;

    /*********************   WARNING  **********************/
    /*                                                     */
    /*                 DYNAMIC SIMULATION CORE             */
    /*                                                     */
    /*********************   WARNING  **********************/
    unsigned w=0;
    for (w = 0 ; w < params.num_walkers; w++)
    {
        //flag in case there was any error with the particle.
        back_tracking = false;

        walker.setIndex(w);

        // Initialize the walker initial position
        iniWalkerPosition();

        // Selects only obstacles that are close enough to collide and the ones inside a collision sphere
        initWalkerObstacleIndexes();

        //Initial position;
        walker.setRealPosLog(walker.pos_r,0);
        walker.setVoxPosLog (walker.pos_v,0);

        //cout << "\n Iniatial postionl                                                                  ";
        //cout << walker.ini_pos[0] << " "  << walker.ini_pos[1] << " "  << walker.ini_pos[2] << endl;

        for(unsigned t = 1 ; t <= params.num_steps; t++) //T+1 steps in total (avoid errors)
        {
            //Get the time step in milliseconds
            getTimeDt(last_time_dt,time_dt,l,dataSynth,t,time_step);

            //Generates a random oriented step of size l
            generateStep(step,l);

            // Moves the particle. Checks collision and handles bouncing.
            try{
                updateWalkerPosition(step);
            }
            catch(Sentinel::ErrorCases error){

                // Possible errors, or numerical un-handed cases should end here.
                sentinela.deportationProcess(walker,w,t,back_tracking,params,id);

                if ( (error == Sentinel::ErrorCases::stuck) || (error == Sentinel::ErrorCases::crossed)){
                    //w--;
                    break;
                }

                if ( error == Sentinel::rejected  )
                    continue;
            }

            // Saves the final particle position after bouncing in the time t.
            walker.setRealPosLog(walker.pos_r,t);
            walker.setVoxPosLog (walker.pos_v,t);

            //updates the collision neighborhood (if any)
            updateCollitionSphere(t);

            walker.steps_count++;
            walker.rejection_count = 0;
        }// end for t

        if(!back_tracking)
            if(finalPositionCheck()){
                back_tracking=true;
                sentinela.illegal_count++;
                w--;
            }

        //If there was an error, we don't compute the signal or write anything.
        if(back_tracking){
            continue;
        }

        //updates the phase shift.
        if(dataSynth)
            dataSynth->update_phase_shift(this->time_step,walker.pos_r_log);

        //Update de DWI signal
        if(dataSynth)
            dataSynth->update_DWI_signal(walker);

        //Write the positions.
        trajectory.writePosition(walker.pos_r_log);

        if(params.log_propagator){
            //Update Propagator
            updatePropagator(walker.pos_r_log);
        }

        //Displays the remained expected time and check for the time limit.
        if(expectedTimeAndMaxTimeCheck(w)){
            cout << "\n" << SH_BG_LIGHT_YELLOW <<  "[Warning]" << SH_DEFAULT << "  Sim: " << id << " "
                 << "Max time limit reached: Simulation halted after "<< ++w << " spins" << endl;
            break;
        }

    }// for w


    /*********************   WARNING  **********************/
    /*                                                     */
    /*         END OF THE DYNAMIC SIMULATION CORE          */
    /*                                                     */
    /*********************   WARNING  **********************/

    num_simulated_walkers = w;

    if(num_simulated_walkers<= params.num_walkers){

        trajectory.reWriteHeaderFile(num_simulated_walkers);

        this->params.num_walkers = num_simulated_walkers;
    }

    if(params.log_propagator){
        normalizePropagator(num_simulated_walkers);
    }

    //computes the ICVF from the initialization.
    computeICVF();

    // Info display.
    time(&now);
    second_passed = difftime(now,start);
    if(params.verbatim)
        SimErrno::info(" Sim: " + to_string(id) + " Simulation ended after: "
                       + secondsToMinutes(second_passed) + " seconds" ,cout);


    // Writes the final DWI signal, and the phase shift.
    if(params.log_opp)
        writeDWSignal(dataSynth);

    return;
}

DynamicsSimulation::~DynamicsSimulation() {
    if(iniPos.is_open())
        iniPos.close();
}

/**
 * @param conf_file_path
 * @return void
 */
void DynamicsSimulation::readConfigurationFile(std::string conf_file_path) {
    params.readSchemeFile(conf_file_path);
}

/**
 * @return void
 */
void DynamicsSimulation::generateStep(Vector3d & step, double l) {

    if(walker.status == Walker::on_object){
        step = walker.next_direction.normalized();
        return;
    }

    std::uniform_real_distribution<double> dist(0,1);

    /* Unbiased random direction*/
    double theta  = 2.0*M_PI*dist(mt);
    double cosPhi = 2.0*dist(mt)-1.0;
    double cosTh  = cos(theta);
    double sinTh  = sin(theta);
    double sinPhi = sqrt(1.0-cosPhi*cosPhi);

    step(0) = l*cosTh*sinPhi;
    step(1) = l*sinTh*sinPhi;
    step(2) = l*cosPhi;

    step.normalize();

}

void DynamicsSimulation::generateDirectedStep(Vector3d &new_step, Vector3d &direction){

    std::uniform_real_distribution<double> dist(0,1);

    /* Unbiased random direction*/
    double theta  = 2.0*M_PI*dist(mt);
    double cosPhi = 2.0*dist(mt)-1.0;
    double cosTh  = cos(theta);
    double sinTh  = sin(theta);
    double sinPhi = sqrt(1.0-cosPhi*cosPhi);

    new_step(0) = cosTh*sinPhi;
    new_step(1) = sinTh*sinPhi;
    new_step(2) = cosPhi;

    new_step.normalize();

    double rn = direction.dot(new_step);

    if(rn<0.0)
        new_step*=-1.0;
}

/**
 * @return True if a bouncing in needed.
 */
bool DynamicsSimulation::updateWalkerPosition(Eigen::Vector3d& step) {


    //new step to take
    Vector3d bounced_step = step.normalized(),end_point;
    Vector3d previous_real_position, previous_voxel_position, real_pos, voxel_pos;

    //Save the walker initial position.
    walker.getVoxelPosition(previous_voxel_position);
    walker.getRealPosition(previous_real_position);

    // Collision instance to save manage the collision (in Spanish).
    Collision colision;

    // True when the particle needs to be bounced and updates.
    bool bounced = false;
    bool update_walker_status  = false;

    // Maximum displacement. Is updated after each bouncing (if any)
    double tmax = step_lenght;

    // Clears the status of the sentinel.
    sentinela.clear();

    unsigned bouncing_count = 0;
    do{
        bounced = false;
        bouncing_count++;
        
        // Checks the number of bouncing per step.
        walker.steps_count++;

        // True if there was a collision and the particle needs to be bounced.
        update_walker_status |= checkObstacleCollision(bounced_step, tmax, end_point, colision);

        // Updates the position and bouncing direction.
        if(update_walker_status){
            bounced = updateWalkerPositionAndHandleBouncing(bounced_step,tmax,colision);
            // restarts the variables.
            update_walker_status = false;
            colision.type = Collision::null;
        }
        else{
            if (colision.type == Collision::null){
                // clear from previous status
                walker.status = Walker::free;
                walker.next_direction = {0,0,0};
            }
        }
        sentinela.checkErrors(walker,params,((*plyObstacles_list).size() == 0),bouncing_count);

    }while(bounced);


    if(tmax >= 0.0){

        // Update the walker position after the bouncing (or not)
        walker.getRealPosition(real_pos);
        walker.setRealPosition(real_pos  + tmax*bounced_step);

        walker.getVoxelPosition(voxel_pos);
        walker.setVoxelPosition(voxel_pos+ tmax*bounced_step);
    }

    return false;
}

bool DynamicsSimulation::checkObstacleCollision(Vector3d &bounced_step,double &tmax, Eigen::Vector3d& end_point,Collision& colision)
{

    Collision colision_tmp;
    colision_tmp.type = Collision::null;
    colision_tmp.t = INFINITY_VALUE;

    //Origin O
    Eigen::Vector3d ray_origin;
    walker.getVoxelPosition(ray_origin);

    //To keep track of the closest collision
    double max_collision_distance = tmax;


    // The collision checks the three possible obstacles in this order: Voxel, Cylinders, PLY.
    // The closest collision is kept at the end.
    
    //Check Voxel limits
    for(unsigned int i = 0 ; i < voxels_list.size(); i++ )
    {
        voxels_list[i].CheckCollision(walker,bounced_step,tmax,colision_tmp);
        handleCollisions(colision,colision_tmp,max_collision_distance,i);
    }

    //For each Cylinder Obstacles
    for(unsigned int i = 0 ; i < walker.cylinders_collision_sphere.small_sphere_list_end; i++ )
    {
        unsigned index = walker.cylinders_collision_sphere.collision_list->at(i);

        (*cylinders_list)[index].checkCollision(walker,bounced_step,tmax,colision_tmp);
        handleCollisions(colision,colision_tmp,max_collision_distance,index);
    }

    //For each Spehere Obstacle
    for(unsigned int i = 0 ; i < walker.spheres_collision_sphere.small_sphere_list_end; i++ )
    {
        unsigned index = walker.spheres_collision_sphere.collision_list->at(i);

        (*spheres_list)[index].checkCollision(walker,bounced_step,tmax,colision_tmp);
        handleCollisions(colision,colision_tmp,max_collision_distance,index);
    }

    //For each PLY Obstacles
    for(unsigned int i = 0 ; i < walker.ply_collision_sphere.collision_list->size(); i++ )
    {

        if((walker.in_ply_index >=0) && walker.in_ply_index != int(i)){
            continue;
        }

        (*plyObstacles_list)[i].checkCollision(walker,bounced_step,tmax,colision_tmp, walker.ply_collision_sphere.collision_list->at(i),
                                            walker.ply_collision_sphere.small_sphere_list_end[i]);

        handleCollisions(colision,colision_tmp,max_collision_distance,i);
    }


    return colision.type != Collision::null;
}


void DynamicsSimulation::handleCollisions(Collision &colision, Collision &colision_2, double &max_collision_distance, unsigned indx)
{
    // nothing to do;
    if (colision_2.type == Collision::null)
        return;

    colision_2.obstacle_ind = int(indx);

    if (colision.type == Collision::hit || colision.type == Collision::boundary){
        if(colision_2.doIHaveMorePiorityThan(colision)){
            colision = colision_2;
            max_collision_distance = colision_2.t;
            colision.obstacle_ind = int(indx);
        }
        return;
    }

    if(colision.type == Collision::near ){
        if (colision_2.type == Collision::hit || colision_2.type == Collision::boundary){
            colision = colision_2;
            max_collision_distance = colision_2.t;
            colision.obstacle_ind = int(indx);
        }
        return;
    }

    // if we get here means that colision.type = 'null'
    if(colision_2.type == Collision::near){

        colision = colision_2;
        colision.obstacle_ind = int(indx);

        return;
    }

    colision = colision_2;
}


void DynamicsSimulation::mapWalkerIntoVoxel(Eigen::Vector3d& bounced_step, Collision &colision,double barrier_thicknes)
{

    walker.setRealPosition(walker.pos_r + colision.t*bounced_step);

    Eigen::Vector3d voxel_pos = walker.pos_v + (colision.t)*bounced_step;

    bool mapped = false;
    for(int i = 0 ; i < 3; i++)
    {
        if ( fabs(voxel_pos[i] -  voxels_list[0].min_limits[i]) <= EPS_VAL){
            voxel_pos[i] = voxels_list[0].max_limits[i];
            mapped = true;
        }
        else if ( fabs(voxel_pos[i] - voxels_list[0].max_limits[i]) <= EPS_VAL){
            voxel_pos[i] = voxels_list[0].min_limits[i];
            mapped = true;
        }
    }

    walker.setVoxelPosition(voxel_pos);

    if (mapped){
        initWalkerObstacleIndexes();
    }
}

void DynamicsSimulation::getTimeDt(double &last_time_dt, double &time_dt, double &l, SimulableSequence* dataSynth, unsigned t, double time_step)
{
    last_time_dt = time_step*(t-1);
    time_dt = time_step*(t);

    if(dataSynth){
        if(dataSynth->dynamic){
            last_time_dt = dataSynth->time_steps[t-1];
            time_dt = dataSynth->time_steps[t];
            l = sqrt(6.0*(params.diffusivity*(time_dt - last_time_dt)));
        }
    }
}

bool DynamicsSimulation::updateWalkerPositionAndHandleBouncing(Vector3d &bounced_step, double &tmax, Collision &colision)
{

    // To avoid numerical errors.
    double min_step_length = barrier_tickness;

    Eigen::Vector3d real_pos, voxel_pos;
    walker.getRealPosition(real_pos);
    walker.getVoxelPosition(voxel_pos);

    bool bounced =false;

    //Sets the status of the walker
    walker.status = Walker::free;


    if (tmax < min_step_length)
    {
        tmax = 0;
        return false;
    }

    if(colision.type == Collision::hit && colision.col_location != Collision::voxel)
    {
        bounced = true;
        double displ = colision.t;

        //If the collision was really close we can't trust the normal direction;
        if(colision.t < 1e-10 && walker.status != walker.bouncing){
            sentinela.rejected_step = true;
            Eigen::Vector3d direction = -step;
            generateDirectedStep(walker.next_direction,direction);
            walker.status = Walker::on_object;
            return false;
        }

        walker.status = Walker::bouncing;

        tmax -= displ;

        // Labels the walker w/r it's orientation.
        if(colision.col_location == Collision::inside){
            walker.intra_extra_consensus--;
            walker.location = Walker::intra;
        }
        if(colision.col_location == Collision::outside){
            walker.intra_extra_consensus++;
            walker.location = Walker::extra;
        }
        if(walker.initial_location == Walker::unknown){
            walker.initial_location = walker.location;
        }

        //We update the position.
        walker.setRealPosition (real_pos   + displ*bounced_step);
        walker.setVoxelPosition(voxel_pos  + displ*bounced_step);

        bounced_step = colision.bounced_direction;
    }
    else if(colision.type == Collision::hit && colision.col_location == Collision::voxel)
    {

        bounced = true;

        walker.status = Walker::on_voxel;

        mapWalkerIntoVoxel(bounced_step,colision,barrier_tickness);
        bounced_step = colision.bounced_direction;
        tmax-=colision.t;
    }
    else if(colision.type == Collision::near){
        //sentinela.rejected_step   = true;
        Eigen::Vector3d direction = -bounced_step; //WARNING: deberiamos usar el bounced step.
        generateDirectedStep(walker.next_direction,direction);
        walker.status = Walker::on_object;
        return false;
    }
    else if(colision.type == Collision::degenerate){
        sentinela.rejected_step = true;
        Eigen::Vector3d direction = -step;
        generateDirectedStep(walker.next_direction,direction);
        walker.status = Walker::on_object;
        return false;
    }

    return bounced;
}



void DynamicsSimulation::setDuration(const double &duration)
{
    params.sim_duration = duration;
    trajectory.dyn_duration = duration;
}

void DynamicsSimulation::setWalkersNum(const unsigned &N)
{
    params.num_walkers  = N;
    trajectory.N = N;
}

void DynamicsSimulation::setStepsNum(const unsigned &T)
{
    params.num_steps = T;
    trajectory.T = T;
}


