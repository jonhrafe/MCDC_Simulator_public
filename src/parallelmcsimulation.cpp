#include "parallelmcsimulation.h"
#include <iomanip>
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include "cylindergammadistribution.h"
#include "spheregammadistribution.h"

//* Auxiliare method to split words in a line using the spaces*//
template<typename Out>
void split_(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split_(s, delim, std::back_inserter(elems));
    return elems;
}


using namespace std;

ParallelMCSimulation::ParallelMCSimulation(std::string config_file)
{

    mean_second_passed = 0;

    total_sim_particles = 0;

    SimErrno::checkConfigurationFile(config_file.c_str());

    params.readSchemeFile(config_file);

    SimErrno::checkSimulationParameters(params);
    //printSimulationInfo();
    initializeUnitSimulations();

    SimErrno::printSimulatinInfo(params,std::cout);

}

ParallelMCSimulation::ParallelMCSimulation(Parameters &params)
{
    this->params = params;
    mean_second_passed = 0;
    total_sim_particles = 0;
    SimErrno::checkSimulationParameters(params);
    this->params = params;
    initializeUnitSimulations();
    SimErrno::printSimulatinInfo(params,std::cout);
    icvf=0;
}

ParallelMCSimulation::~ParallelMCSimulation()
{
    for(unsigned i = 0; i < simulations.size(); i++){
        delete simulations[i];
    }
}

void ParallelMCSimulation::startSimulation()
{
    cout<<setfill('-');
    cout << SH_FG_PURPLE << "/********************   MC/DC Simulation START:  *************************/" << SH_DEFAULT << "\n";

    for(unsigned int i =0; i < simulations.size(); i++){
        sim_threads.push_back(std::thread(&MCSimulation::startSimulation,(simulations[i])));
    }

    for(unsigned i =0; i < simulations.size(); i++){
        sim_threads[i].join();
    }

    for(unsigned i =0; i < simulations.size(); i++){
        mean_second_passed  += simulations[i]->dynamicsEngine->second_passed/double(simulations.size());
        total_sim_particles += simulations[i]->dynamicsEngine->num_simulated_walkers;
    }

    cout<<setfill('-');

    SimErrno::info( " All " + to_string(params.num_proc) +  " simulations ended after: " + DynamicsSimulation::secondsToMinutes(mean_second_passed)
                    + " in average",cout  );

    SimErrno::info( " Joining resulting data... ",cout  );

    jointResults();

    SimErrno::info( " Done.",cout  );

    ofstream out(params.output_base_name+"_simulation_info.txt", std::ofstream::app);
    SimErrno::info("All " + to_string(params.num_proc) +  " simulations ended after: " + DynamicsSimulation::secondsToMinutes(mean_second_passed)
                   + " in average",out,false);
    SimErrno::info("Number of particles labeled as stuck: "        + to_string(stuck_count)  ,out,false);
    SimErrno::info("Number of particles eliminated due crossings: "+ to_string(illegal_count),out,false);

    if(params.max_simulation_time > 0){
        SimErrno::info("Number of simulated particles: "+ to_string(total_sim_particles),out,false);
        SimErrno::info("Mean simulation speed: "+ to_string( unsigned(params.num_steps*total_sim_particles/mean_second_passed)) + " steps/second",out,false);
    }
    else{
        SimErrno::info("Mean simulation speed: "+ to_string( unsigned(params.num_steps*params.num_walkers/mean_second_passed)) + " steps/second",out,false);
    }
    if(params.voxels_list.size()>0){
        string message = "Voxel limits:";
        SimErrno::info(message,out,false);

        out << std::setprecision(10) << "( " <<  params.voxels_list[0].first[0] <<
                " " << params.voxels_list[0].first[1]<<
                " " << params.voxels_list[0].first[2]<< " )\n";

        out << std::setprecision(10) << "( " <<  params.voxels_list[0].second[0]<<
                " " << params.voxels_list[0].second[1] <<
                " " << params.voxels_list[0].second[2] << " )\n";
    }

    if(params.custom_sampling_area){
        out << std::fixed;
        string message="Custom sampled area (mm):";
        SimErrno::info(message,out,false);
        out << std::setprecision(10) << "( " <<  params.min_sampling_area[0] <<
                " " << params.min_sampling_area[1]<<
                " " << params.min_sampling_area[2]<< " )\n";

        out << std::setprecision(10) << "( " <<  params.max_sampling_area[0]<<
                " " << params.max_sampling_area[1] <<
                " " << params.max_sampling_area[2] << " )\n";
    }


    if((params.custom_sampling_area or params.voxels_list.size()>0) and (params.computeVolume)){
        string message="Estimated Intra-axonal volume from sampling (mm^3):";
        out << std::scientific;
        SimErrno::info(message,out,false);
        out << std::scientific;
        out << aprox_volumen << endl;
    }


    if(params.hex_cyl_packing == true){
        SimErrno::info("Hex. packing radius: "+ to_string( params.hex_packing_radius)    ,out,false);
        SimErrno::info("Hex. packing separation: " +  to_string( params.hex_packing_separation) ,out,false);
        SimErrno::info("Hex. packing ICVF: " +  to_string( params.hex_packing_icvf) ,out,false);
    }
    out.close();
}


void ParallelMCSimulation::initializeUnitSimulations()
{

    //Build anything that needs to be syn between simulations.
    specialInitializations();

    // The number of walker N devided between all the processes


    unsigned N_per_sim = params.num_walkers/params.num_proc;

    for(unsigned i = 0; i < params.num_proc-1; i++){

        //Parameters for each simulation simulation
        Parameters params_temp       = params;
        params_temp.num_walkers      = N_per_sim;
        params_temp.output_base_name+= "_"+std::to_string(i);


        MCSimulation* simulation_ = new MCSimulation(params_temp);
        simulation_->plyObstacles_list = &this->plyObstacles_list;
        simulation_->cylinders_list = &this->cylinders_list;
        simulation_->sphere_list = &this->spheres_list;
        simulation_->dynamicsEngine->print_expected_time = 0;
        simulations.push_back(simulation_);

        if(params.verbatim)
            SimErrno::info( " Sim: " + to_string(simulation_->dynamicsEngine->id) + " Initialized",cout);
    }

    //We set thet remaining walkers so it sums up the desired total

    //Parameters for each simulation simulation
    Parameters params_temp = params;
    params_temp.num_walkers = params.num_walkers - N_per_sim *(params.num_proc-1);
    params_temp.output_base_name+= "_"+std::to_string(params.num_proc-1);


    MCSimulation* simulation_ = new MCSimulation(params_temp);
    //simulation_->dynamicsEngine->print_expected_time = 0;
    simulation_->plyObstacles_list = &this->plyObstacles_list;
    simulation_->cylinders_list = &this->cylinders_list;
    simulation_->sphere_list    = &this->spheres_list;
    simulations.push_back(simulation_);

    if(params.verbatim)
        SimErrno::info( " Sim: " + to_string(simulation_->dynamicsEngine->id) + " Initialized",cout);
}

void ParallelMCSimulation::jointResults()
{

    //====================================================================================================
    //Complete DWI data.
    //====================================================================================================

    //Writes the ouput data
    if(simulations[0]->dataSynth){


        std::string outDWI    = params.output_base_name  + "_DWI.txt";
        std::string outDWI_intra    = params.output_base_name  + "_DWI_intra.txt";
        std::string outDWI_extra    = params.output_base_name  + "_DWI_extra.txt";
        std::string outDWIi   = params.output_base_name  + "_DWI_img.txt";
        std::string outPhase  = params.output_base_name  + "_phase_shift.txt";

        std::string boutDWI    = params.output_base_name  + "_DWI.bfloat";
        std::string boutDWI_intra    = params.output_base_name  + "_DWI_intra.bfloat";
        std::string boutDWI_extra    = params.output_base_name  + "_DWI_extra.bfloat";
        std::string boutDWIi   = params.output_base_name  + "_DWI_img.bfloat";
        std::string boutPhase  = params.output_base_name  + "_phase_shift.bfloat";

        std::ofstream phase_out,phase_bout, dwi_out, dwii_out, dwi_bout, dwii_bout,
                      dwi_intra_out,dwi_extra_out,dwi_intra_bout,dwi_extra_bout;

        if(params.log_phase_shift){
            if(params.write_bin)
                phase_bout.open(boutPhase, std::ofstream::binary);   //phase shift (binary)

            if(params.write_txt)
                phase_out.open(outPhase, std::ofstream::out);     //phase shift  (txt)
        }

        if(params.write_txt){
            dwi_out.open(outDWI  ,std::ofstream::out);       //real part


            if(params.img_signal == true)
             dwii_out.open(outDWIi,std::ofstream::out);       //img part

            if(params.separate_signals == true){
                dwi_intra_out.open(outDWI_intra ,std::ofstream::out);       //intra part
                dwi_extra_out.open(outDWI_extra ,std::ofstream::out);       //extra part
            }
        }

        if(params.write_bin){
            dwi_bout.open(boutDWI  ,std::ofstream::binary);       //real part (binary)

            if(params.img_signal == true)
                {dwii_bout.open(boutDWIi,std::ofstream::binary);}       //img part  (binary)

            if(params.separate_signals == true){
                dwi_intra_bout.open(boutDWI_intra ,std::ofstream::binary);       //intra part
                dwi_extra_bout.open(boutDWI_extra ,std::ofstream::binary);       //extra part
            }
        }


        for (unsigned i = 0 ; i < simulations[0]->dataSynth->DWI.size(); i++ )
        {
            double DWI   = 0;
            double DWIi  = 0;
            double DWI_intra   = 0;
            double DWI_extra   = 0;

            std::vector<int> phase(3600, 0);
            for(unsigned p = 0; p < params.num_proc; p++)
            {
                DWI+=  simulations[p]->dataSynth->DWI[i];       //real part

                if(params.img_signal == true)
                    DWIi+= simulations[p]->dataSynth->DWIi[i];      //img  part

                if(params.separate_signals){
                    DWI_intra+=  simulations[p]->dataSynth->DWI_intra[i];      //intra part
                    DWI_extra+=  simulations[p]->dataSynth->DWI_extra[i];      //extra part
                }

                // for each histogram bin, 36000 fixed size
                for(unsigned h = 0; h < 3600; h++)
                {
                    phase[h]+=simulations[p]->dataSynth->phase_shift_distribution(i,h);
                }

            }


            if(params.write_txt){
                dwi_out  << DWI << std::endl;

                if(params.img_signal == true)
                    dwii_out << DWIi << std::endl;

                if(params.separate_signals){
                    dwi_intra_out << DWI_intra << std::endl;
                    dwi_extra_out << DWI_extra << std::endl;
                }
            }

            if(params.write_bin){
                float holder = float(DWI);
                dwi_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));

                if(params.img_signal == true)
                {
                    holder = float(DWIi);
                    dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                }

                if(params.separate_signals){
                    float holder = float(DWI_intra);
                    dwi_intra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                    holder = float(DWI_extra);
                    dwi_extra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                }
            }

            if(params.log_phase_shift){

                if(params.write_txt){
                    //write the full histogram for the gradient i
                    for(unsigned h = 0; h < 3600; h++)
                    {
                        phase_out << phase[h];
                        if(h < 3600-1)
                            phase_out << " ";
                        else
                            phase_out << endl;
                    }
                }

                if(params.write_bin){
                    //write the full histogram for the gradient i
                    for(unsigned h = 0; h < 3600; h++)
                    {
                        float holder = phase[h];
                        phase_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                    }
                }
            }
        }

        if(params.write_txt){
            dwi_out.close();

            if(params.img_signal == true)
                dwii_out.close();

            if(params.separate_signals){
                dwi_intra_out.close();
                dwi_extra_out.close();
            }
        }
        if(params.write_bin){
            dwi_bout.close();

            if(params.img_signal == true)
                dwii_bout.close();

            if(params.separate_signals){
                dwi_intra_bout.close();
                dwi_extra_bout.close();
            }
        }

        if(params.log_phase_shift)
            phase_out.close();
    }

    stuck_count = illegal_count = 0;


     // ICVF for subdivisions

    for(unsigned i = 0 ; i < simulations.size();i++){
        stuck_count   += simulations[i]->dynamicsEngine->sentinela.stuck_count;
        illegal_count += simulations[i]->dynamicsEngine->sentinela.illegal_count;
        icvf+=           simulations[i]->dynamicsEngine->num_simulated_walkers*simulations[i]->dynamicsEngine->icvf;
    }

    icvf/= float(this->total_sim_particles);

    if(params.custom_sampling_area == false and params.voxels_list.size()>0){
        for(auto i = 0; i<3;i++){
           params.min_sampling_area[i]= params.voxels_list[0].first[i];
           params.max_sampling_area[i]= params.voxels_list[0].second[i];
        }
    }

    float sampling_volume =1;

    for (auto i =0; i < 3;i++ ){
        sampling_volume*= (params.max_sampling_area[i]-params.min_sampling_area[i]);
    }


    aprox_volumen = (icvf)*sampling_volume;



    //====================================================================================================
    //Subdivision results join.
    //====================================================================================================


    if(simulations[0]->dataSynth && this->params.subdivision_flag)
    {
        std::string out_vox_DWI    = params.output_base_name + "_voxels_DWI.txt";
        std::string out_vox_DWIi   = params.output_base_name + "_voxels_DWI_img.txt";
        std::string out_vox_DWI_intra   = params.output_base_name + "_voxels_DWI_intra.txt";
        std::string out_vox_DWI_extra   = params.output_base_name + "_voxels_DWI_extra.txt";

        std::string bout_vox_DWI    = params.output_base_name  + "_voxels_DWI.bfloat";
        std::string bout_vox_DWIi   = params.output_base_name  + "_voxels_DWI_img.bfloat";
        std::string bout_vox_DWI_intra   = params.output_base_name + "_voxels_DWI_intra.bfloat";
        std::string bout_vox_DWI_extra   = params.output_base_name + "_voxels_DWI_extra.bfloat";

        std::ofstream dwi_out, dwii_out, dwi_bout, dwii_bout,
                      dwi_intra_out, dwi_extra_out, dwi_intra_bout,dwi_extra_bout ;

        if(params.write_txt){
            dwi_out.open(out_vox_DWI  ,std::ofstream::out);       //real part

            if(params.img_signal == true)
                dwii_out.open(out_vox_DWIi,std::ofstream::out);       //img part

            if(params.separate_signals){
                dwi_intra_out.open(out_vox_DWI_intra  ,std::ofstream::out);      //intra part
                dwi_extra_out.open(out_vox_DWI_extra ,std::ofstream::out);       //extra part
            }
        }

        if(params.write_bin){
            dwi_bout.open(bout_vox_DWI  ,std::ofstream::binary);       //real part (binary)
            if(params.img_signal == true)
                dwii_bout.open(bout_vox_DWIi,std::ofstream::binary);       //img part  (binary)

            if(params.separate_signals){
                dwi_intra_bout.open(bout_vox_DWI_intra ,std::ofstream::binary);       //intra part (binary)
                dwi_extra_bout.open(bout_vox_DWI_extra ,std::ofstream::binary);       //extra part (binary)
            }

        }


        std::string outPDen   = params.output_base_name  + "_volFractions.txt";
        std::ofstream PDen_out(outPDen,std::ofstream::out);   //Particle Density map

        //## HEADER
        if(params.write_txt){

            // #Num_voxels
            dwi_out  << simulations[0]->dataSynth->subdivisions.size() << endl;

            if(params.img_signal == true)
            {
                dwii_out << simulations[0]->dataSynth->subdivisions.size() << endl;
                dwii_out << simulations[0]->dataSynth->DWI.size() << endl;
            }

            // #Size_DWI
            dwi_out  << simulations[0]->dataSynth->DWI.size() << endl;

            if(params.separate_signals){
                // #Num_voxels
                dwi_intra_out << simulations[0]->dataSynth->subdivisions.size() << endl;
                dwi_extra_out << simulations[0]->dataSynth->subdivisions.size() << endl;
                // #Size_DWI
                dwi_intra_out  << simulations[0]->dataSynth->DWI.size() << endl;
                dwi_extra_out  << simulations[0]->dataSynth->DWI.size() << endl;
            }
        }
        //## HEADER (Binary)
        if(params.write_bin){
            // #Num_voxels
            float holder = float(simulations[0]->dataSynth->subdivisions.size());
            dwi_bout.write (reinterpret_cast<char *>(&holder), sizeof(float));

            if(params.img_signal == true)
                dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));

            // #Size_DWI
            holder = float(simulations[0]->dataSynth->DWI.size());
            dwi_bout.write (reinterpret_cast<char *>(&holder), sizeof(float));

            if(params.img_signal == true)
                dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));

            if(params.separate_signals){
                // #Num_voxels
                float holder = float(simulations[0]->dataSynth->subdivisions.size());
                dwi_intra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                dwi_extra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                // #Size_DWI
                holder = float(simulations[0]->dataSynth->DWI.size());
                dwi_intra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                dwi_extra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
            }
        }

        for(uint s = 0; s < simulations[0]->dataSynth->subdivisions.size(); s++)
        {
            double density_intra = 0;
            double density_extra = 0;
            double density = 0;
            for(unsigned p = 0; p < params.num_proc; p++){

                density_intra +=  simulations[p]->dataSynth->subdivisions[s].density_intra;
                density_extra +=  simulations[p]->dataSynth->subdivisions[s].density_extra;
                density       +=  simulations[p]->dataSynth->subdivisions[s].density;
            }

            PDen_out << density << endl;
            PDen_out << density_intra << endl;
            PDen_out << density_extra << endl;


            for (unsigned i = 0 ; i < simulations[0]->dataSynth->DWI.size(); i++ )
            {
                double sub_DWI   = 0;
                double sub_DWIi  = 0;
                double sub_DWI_intra  = 0;
                double sub_DWI_extra  = 0;

                for(unsigned p = 0; p < params.num_proc; p++)
                {
                    sub_DWI +=  simulations[p]->dataSynth->sub_DWI[s][i];       //real part

                    if(params.img_signal == true)
                        sub_DWIi+=  simulations[p]->dataSynth->sub_DWIi[s][i];      //img  part

                    if(params.separate_signals){
                        sub_DWI_intra +=  simulations[p]->dataSynth->sub_DWI_intra[s][i];      //intra part
                        sub_DWI_extra +=  simulations[p]->dataSynth->sub_DWI_extra[s][i];      //extra  part
                    }
                }


                if(params.write_txt){
                    dwi_out  << sub_DWI  << std::endl;

                    if(params.img_signal == true)
                        dwii_out << sub_DWIi << std::endl;

                    if(params.separate_signals){
                        dwi_intra_out << sub_DWI_intra << std::endl;
                        dwi_extra_out << sub_DWI_extra << std::endl;
                    }
                }
                if(params.write_bin){
                    float holder = float(sub_DWI);
                    dwi_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));

                    if(params.img_signal == true)
                    {
                        holder = float(sub_DWIi);
                        dwii_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                    }

                    if(params.separate_signals){
                        float holder = float(sub_DWI_intra);
                        dwi_intra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                        holder = float(sub_DWI_extra);
                        dwi_extra_bout.write(reinterpret_cast<char *>(&holder), sizeof(float));
                    }
                }
            }
        }

        if(params.write_txt){
            dwi_out.close();

            if(params.img_signal == true)
                dwii_out.close();

            if(params.separate_signals){
                dwi_intra_out.close();
                dwi_extra_out.close();
            }
        }
        if(params.write_bin){
            dwi_bout.close();
            if(params.img_signal == true)
                dwii_bout.close();

            if(params.separate_signals){
                dwi_intra_bout.close();
                dwi_extra_bout.close();
            }
        }

        PDen_out.close();
    }


    //====================================================================================================
    // NEW Join propagator join solution
    //====================================================================================================

    if (this->params.log_propagator){
        uint num_sims = uint(simulations.size());
        for (uint i =0; i< num_sims; i++){
            for (uint j = 0; j < this->simulations[i]->dynamicsEngine->propagator.num_times; j++){
                for (uint k = 0; k < this->simulations[i]->dynamicsEngine->propagator.num_dirs; k++){
                    if(i == 0){
                        this->simulations[0]->dynamicsEngine->propagator.propagator_log[j][k]/=float(num_sims);
                        }
                    else{
                        this->simulations[0]->dynamicsEngine->propagator.propagator_log[j][k]+= this->simulations[i]->dynamicsEngine->propagator.propagator_log[j][k]/float(num_sims);
                    }
                }
            }
        }
    }


    std::string prop_file_name = params.output_base_name + "_propagator.txt";

    if (this->params.log_propagator)
        this->simulations[0]->dynamicsEngine->writePropagator(prop_file_name);
}

void ParallelMCSimulation::specialInitializations()
{
    addObstacleConfigurations();
    addObstaclesFromFiles();

    if(params.number_subdivisions>1){
        params.addSubdivisions();
    }

   //std::cout << params.img_signal << std::endl;
    for(unsigned i = 0; i < params.PLY_files.size(); i++){
        //std::cout << i << std::endl;
        //plyObstacles_list.push_back(PLYObstacle(params.PLY_files[i],centers,max_distance,params.PLY_scales[i]));
        plyObstacles_list.push_back(PLYObstacle(params.PLY_files[i],params.PLY_scales[i]));
        plyObstacles_list.back().id=i;
        plyObstacles_list.back().percolation = params.PLY_percolation[i];
    }
}


void ParallelMCSimulation::addObstaclesFromFiles()
{
    for(unsigned i = 0; i < params.cylinders_files.size(); i++){

        bool z_flag = false;
        std::ifstream in(params.cylinders_files[i]);

        if(!in){
            //std::cout << "\033[1;37m[INFO]\033[0m Sim: " << count << " " << "[ERROR] Unable to open:" << params.cylinders_files[i] << std::endl;
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

        in.open(params.cylinders_files[i]);

        if(z_flag){
            double x,y,z,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> r)
            {
                cylinders_list.push_back(Cylinder(Eigen::Vector3d(x,y,z),Eigen::Vector3d(x,y,z+1.0),r,scale));
            }
            in.close();
        }
        else{
            double x,y,z,ox,oy,oz,r;
            double scale;
            in >> scale;
            while (in >> x >> y >> z >> ox >> oy >> oz >> r)
            {
                cylinders_list.push_back(Cylinder(Eigen::Vector3d(x,y,z),Eigen::Vector3d(ox,oy,oz),r,scale));
            }
            in.close();
        }
    }

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

        in.open(params.spheres_files[i]);
        double x,y,z,r;
        double scale;
        in >> scale;

        while (in >> x >> y >> z >> r)
        {
            spheres_list.push_back(Sphere(Eigen::Vector3d(x,y,z),r,scale));
        }
        in.close();
    }
}

void ParallelMCSimulation::addObstacleConfigurations()
{

    if(params.hex_cyl_packing){

        double rad = params.hex_packing_radius,sep = params.hex_packing_separation;

        // h = sqrt(3)/2 * sep
        double h =0.8660254037844386*sep;

        cylinders_list.push_back(Cylinder(Eigen::Vector3d(0,0,0),Eigen::Vector3d(0,0,1.0),rad));
        cylinders_list.push_back(Cylinder(Eigen::Vector3d(sep,0,0),Eigen::Vector3d(sep,0,1.0),rad));

        cylinders_list.push_back(Cylinder(Eigen::Vector3d(0,2.0*h,0),Eigen::Vector3d(0,2.0*h,1.0),rad));
        cylinders_list.push_back(Cylinder(Eigen::Vector3d(sep,2.0*h,0),Eigen::Vector3d(sep,2.0*h,1.0),rad));

        cylinders_list.push_back(Cylinder(Eigen::Vector3d(0.5*sep,h,0),Eigen::Vector3d(0.5*sep,h,1.0),rad));

        // To avoid problems with the boundaries
        cylinders_list.push_back(Cylinder(Eigen::Vector3d(-0.5*sep,h,0),Eigen::Vector3d(-0.5*sep,h,1.0),rad));
        cylinders_list.push_back(Cylinder(Eigen::Vector3d(1.5*sep,h,0),Eigen::Vector3d(1.5*sep,h,1.0),rad));

        if(params.voxels_list.size()>0)
            params.voxels_list.clear();

        pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(Eigen::Vector3d(0,0,0),Eigen::Vector3d(sep,2.0*h,2.0*h));
        params.voxels_list.push_back(voxel_);
    }

    if(params.hex_sphere_packing){

         double rad = params.hex_packing_radius,sep = params.hex_packing_separation;
         double h = sqrt(3.)/2.0*sep;
         //double h = 0.866025404*sep;

        // // [0.0, 0.0, 0.0, 1.0],[s, 0.0, 0.0, 1.0],
         spheres_list.push_back(Sphere(Eigen::Vector3d(0.0, 0.0, 0.0),rad));
         spheres_list.push_back(Sphere(Eigen::Vector3d(sep, 0.0, 0.0),rad));

        // // [s/2.0, h, 0.0, 1.0],
         spheres_list.push_back(Sphere(Eigen::Vector3d(sep/2.0, h, 0.0),rad));

        // // [0.0, 2*h, 0.0, 1.0],[s, 2*h, 0.0, 1.0],
         spheres_list.push_back(Sphere(Eigen::Vector3d(0.0, 2.*h, 0.0),rad));
         spheres_list.push_back(Sphere(Eigen::Vector3d(sep, 2.*h, 0.0),rad));

        // // [0.0, h, h, 1.0],[s, h, h, 1.0],
        spheres_list.push_back(Sphere(Eigen::Vector3d(0.0, h, h),rad));
        spheres_list.push_back(Sphere(Eigen::Vector3d(sep, h, h),rad));

        // // [s/2, 0, h, 1.0],[s/2, 2*h, h, 1.0],
        spheres_list.push_back(Sphere(Eigen::Vector3d(sep/2., 0.0, h),rad));
        spheres_list.push_back(Sphere(Eigen::Vector3d(sep/2., 2.0*h, h),rad));

        // // [0.0, h, -h, 1.0],[s, h, -h, 1.0],
        spheres_list.push_back(Sphere(Eigen::Vector3d(0.0, h, -h),rad));
        spheres_list.push_back(Sphere(Eigen::Vector3d(sep,h,-h),rad));

        // // [s/2, 0, -h, 1.0],[s/2, 2*h, -h, 1.0],
        spheres_list.push_back(Sphere(Eigen::Vector3d(sep/2.0 ,0.0 ,-h),rad));
        spheres_list.push_back(Sphere(Eigen::Vector3d(sep/2.0,2.0*h,-h),rad));

        pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(Eigen::Vector3d(0,0,-h),Eigen::Vector3d(sep,2.0*h,h));
        params.voxels_list.push_back(voxel_);

    }


    if(params.gamma_cyl_packing == true){

        string message = "Initialializing Gamma distribution (" + std::to_string(params.gamma_packing_alpha) + ","
                + std::to_string(params.gamma_packing_beta) + ").\n";
        SimErrno::info(message,cout);

        CylinderGammaDistribution gamma_dist(params.gamma_num_obstacles,params.gamma_packing_alpha, params.gamma_packing_beta,params.gamma_icvf
                                             ,params.min_limits, params.max_limits,params.min_obstacle_radii);

        gamma_dist.displayGammaDistribution();

        gamma_dist.createGammaSubstrate();

        params.max_limits = gamma_dist.max_limits;
        params.min_limits = gamma_dist.min_limits;

        if(params.voxels_list.size()<=0){
            pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params.min_limits,params.max_limits);
            params.voxels_list.push_back(voxel_);
        }
        else{
            params.voxels_list[0].first =  params.min_limits;
            params.voxels_list[0].second = params.max_limits;
        }

        string file = params.output_base_name + "_gamma_distributed_cylinder_list.txt";

        ofstream out(file);

        gamma_dist.printSubstrate(out);

        this->cylinders_list = gamma_dist.cylinders;

        //params.cylinders_files.push_back(file);

        out.close();

        SimErrno::info("Done.\n",cout);
    }


    if(params.gamma_sph_packing == true){

        string message = "Initialializing Gamma distribution (" + std::to_string(params.gamma_packing_alpha) + ","
                + std::to_string(params.gamma_packing_beta) + ").\n";
        SimErrno::info(message,cout);

        SphereGammaDistribution gamma_dist(params.gamma_num_obstacles,params.gamma_packing_alpha, params.gamma_packing_beta,params.gamma_icvf
                                             ,params.min_limits, params.max_limits,params.min_obstacle_radii);

        gamma_dist.displayGammaDistribution();

        gamma_dist.createGammaSubstrate();

        params.max_limits = gamma_dist.max_limits;
        params.min_limits = gamma_dist.min_limits;

        if(params.voxels_list.size()<=0){
            pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params.min_limits,params.max_limits);
            params.voxels_list.push_back(voxel_);
        }
        else{
            params.voxels_list[0].first =  params.min_limits;
            params.voxels_list[0].second = params.max_limits;
        }

        string file = params.output_base_name + "_gamma_distributed_sphere_list.txt";

        ofstream out(file);

        gamma_dist.printSubstrate(out);

        this->spheres_list = gamma_dist.spheres;

        //params.cylinders_files.push_back(file);
        out.close();
        SimErrno::info("Done.\n",cout);
    }


}
