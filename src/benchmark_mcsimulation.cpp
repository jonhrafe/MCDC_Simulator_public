#include <benchmark_mcsimulation.h>
#include <iostream>

using namespace std;


Benchmark_mcsimulation::Benchmark_mcsimulation()
{
  benchmark_id = 0; 

}


Benchmark_mcsimulation::Benchmark_mcsimulation(unsigned int id_)
{
  benchmark_id = id_; 
}

void Benchmark_mcsimulation::startBenchmark()
{

  /*
    1. Load bencmark parameters
    2. launch bencmark
    3. Save results
    */


  // 1. Load benchmark parameters

  selectBenchmark(); 


  //2. Launch benchmark

  parallelsimulation =  new ParallelMCSimulation(benchmark_params);    

  cout << parallelsimulation->spheres_list.size() << parallelsimulation->plyObstacles_list.size() << endl;

  parallelsimulation->startSimulation();


  /*
    simulation = new MCSimulation(benchmark_params);    

    simulation->sphere_list         = &this->spheres_list;
    simulation->plyObstacles_list   =  &this->plyObstacles_list;
    simulation->cylinders_list      = &this->cylinders_list;

    simulation->startSimulation();

    // 3. Save results
    writeBenchmark();
    */
}

void Benchmark_mcsimulation::selectBenchmark()
{

  if(benchmark_id==0){
    // No benchmark
    cout << "No benchmark";
    return;
  }

  // Common parameters 
  Parameters params_tmp;

  params_tmp.scale_from_stu       = true;
  params_tmp.write_txt            = true;
  params_tmp.write_bin            = false;
  params_tmp.write_traj           = false;
  params_tmp.separate_signals     = false;


  params_tmp.ini_walker_flag      = "delta";
  params_tmp.scheme_file          = "benchmark/scheme_files/PGSE_sample_scheme.scheme";

  params_tmp.num_proc            = 1; 
  params_tmp.num_walkers         = 10;
  params_tmp.num_steps           = 10; 
  params_tmp.sim_duration        = 50;
  params_tmp.diffusivity         = 2e-9;



  // Benchmark-specific parameters
  if(benchmark_id==1){

    // Sphere case
    cout << "Sphere benchmark" << endl;

    params_tmp.output_base_name     = "benchmark/output/sphere/sphere_gamma_packing";

    // Sphere list
    string path_sphere_list = params_tmp.output_base_name + "_gamma_distributed_sphere_list.txt";
    params_tmp.spheres_files.push_back(path_sphere_list); 


    //loadSpheres(path_sphere_list, params_tmp);

  }

  else if(benchmark_id==2){

      // PLY case
      cout << "PLY benchmark" << endl;

      params_tmp.output_base_name = "benchmark/output/ply/hexagonal_packed_spheres";
      params_tmp.PLY_files.push_back("benchmark/output/ply/hexagonal_packed_spheres.ply");

    }
  else if (benchmark_id==3){
      cout << "PLY cilinders" << endl;

      params_tmp.output_base_name = "benchmark/output/cylinderPly_";
      params_tmp.PLY_files.push_back("benchmark/output/cylinderPly/50L_70icvf.ply");

    }
  else if (benchmark_id==4){
      cout << "List cylinders" << endl;

      params_tmp.output_base_name = "benchmark/output/cylinderList";
      params_tmp.PLY_files.push_back("benchmark/output/ply/cylinderList/50L_70icvf_endpoints.txt");
      params_tmp.cylinders_files.push_back(" ")       /*!< file paths with a list of cilinders obstacles                              */

  }
  params_tmp.traj_file            = params_tmp.output_base_name;                                                

  // Voxel
  string path_voxel = params_tmp.output_base_name + "_voxel.txt";
  loadVoxel(path_voxel, params_tmp);

  benchmark_params = params_tmp;

  return;
}

void Benchmark_mcsimulation::loadSpheres(string path_sphere_list, Parameters &params_tmp){

  params_tmp.spheres_files.push_back(path_sphere_list); 

  std::ifstream in(params_tmp.spheres_files[0]);

  if(!in){
    return;
  }
else{
    double x,y,z,r;
    double scale;
    in >> scale;

    while (in >> x >> y >> z >> r)
      {
        spheres_list.push_back(Sphere(Eigen::Vector3d(x,y,z),r,scale));
      }
    in.close();
    return;
  }
}

void Benchmark_mcsimulation::loadVoxel(string path_voxel, Parameters &params_tmp){

  std::ifstream in(path_voxel);

  if(!in){
    return;
  }
else{
    string tmp="";

    while(in >> tmp){
      if ((params_tmp.str_dist(tmp,"<voxels>") == 0)){

        double x,y,z;

        in >> x >> y >> z;
        params_tmp.min_limits   = {x,y,z};  
        in >> x >> y >> z;
        params_tmp.max_limits   = {x,y,z};

        pair<Eigen::Vector3d,Eigen::Vector3d> voxel_(params_tmp.min_limits, params_tmp.max_limits);
        params_tmp.voxels_list.push_back(voxel_);

      }
      if ((params_tmp.str_dist(tmp,"</voxels>") == 0)){
        in.close();
        return;
      }
    }
    in.close();
    return;
  }
}

void Benchmark_mcsimulation::writeBenchmark()
{
  std::string outDWI          = benchmark_params.output_base_name  + "_DWI.txt";
  std::string outDWI_intra    = benchmark_params.output_base_name  + "_DWI_intra.txt";
  std::string outDWI_extra    = benchmark_params.output_base_name  + "_DWI_extra.txt";
  std::string outDWIi         = benchmark_params.output_base_name  + "_DWI_img.txt";

  std::ofstream dwi_out, dwii_out, dwi_intra_out,dwi_extra_out;

  dwi_out.open(outDWI  ,std::ofstream::out);       //real part
  dwii_out.open(outDWIi, std::ofstream::out);       //img part
  dwi_intra_out.open(outDWI_intra ,std::ofstream::out);       //intra part
  dwi_extra_out.open(outDWI_extra ,std::ofstream::out);       //extra part

  for (unsigned i = 0 ; i < simulation->dataSynth->DWI.size(); i++ ){

    double DWI   = 0;
    double DWIi  = 0;
    double DWI_intra   = 0;
    double DWI_extra   = 0;

    DWI+=  simulation->dataSynth->DWI[i];       //real part

    DWIi+= simulation->dataSynth->DWIi[i];      //img  part

    DWI_intra+=  simulation->dataSynth->DWI_intra[i];      //intra part
    DWI_extra+=  simulation->dataSynth->DWI_extra[i];      //extra part

    dwi_out  << DWI << std::endl;

    dwii_out << DWIi << std::endl;

    dwi_intra_out << DWI_intra << std::endl;
    dwi_extra_out << DWI_extra << std::endl;
  }

  dwi_out.close();
  dwii_out.close();
  dwi_intra_out.close();
  dwi_extra_out.close();

  return;   
}


Benchmark_mcsimulation::~Benchmark_mcsimulation()
{
}
