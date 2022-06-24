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
    1. Load benchmark parameters
    2. launch bencmark
    */

  // 1. Load benchmark parameters
  selectBenchmark(); 


  //2. Launch benchmark
  parallelsimulation =  new ParallelMCSimulation(benchmark_params);    

  parallelsimulation->startSimulation();

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
  params_tmp.ini_walkers_file     = "";
  
  params_tmp.cpu_computation      = true;
  params_tmp.gpu_computation      = true;

  params_tmp.ini_walker_flag      = "delta";
  params_tmp.scheme_file          = "benchmark/scheme_files/10shell_fixed_DdTe.scheme";

  params_tmp.num_proc            = 10; 
  
  /*
  params_tmp.num_walkers         = pow(50, 3);
  params_tmp.num_steps           = 11400; 
  params_tmp.sim_duration        = 57;
  params_tmp.diffusivity         = 2e-6;
  */ 

  params_tmp.num_walkers         = 256*params_tmp.num_proc;
  params_tmp.num_steps           = 1000; 
  params_tmp.sim_duration        = 57;
  params_tmp.diffusivity         = 2e-6;

  // Benchmark-specific parameters

  if(benchmark_id==1){

    // Sphere case
    cout << "Sphere benchmark" << endl;

    params_tmp.output_base_name     = "benchmark/output/sphere/R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing";

    // Sphere list
    string path_sphere_list = params_tmp.output_base_name + "_sphere_list.txt";
    params_tmp.spheres_files.push_back(path_sphere_list); 
        
    }

  else if (benchmark_id==2){
      cout << "Sphere PLY" << endl;

      params_tmp.output_base_name =        "benchmark/output/spherePly/R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing_sphere_list";
      params_tmp.PLY_files.push_back("benchmark/output/spherePly/R_2_R_4_v_50_ICVF_0.57_gaussian_sphere_packing_sphere_list.ply");     /*!< file paths with a list of spheres obstacles                              */;;
      params_tmp.PLY_scales.push_back(1e-3);
      params_tmp.PLY_percolation.push_back(0.0);
  }
  
  else if (benchmark_id==3){
      cout << "PLY cilinders" << endl;

      params_tmp.output_base_name =  "benchmark/output/cylinderPly/50L_70icvf";
      params_tmp.PLY_files.push_back("benchmark/output/cylinderPly/50L_70icvf.ply");
      params_tmp.PLY_scales.push_back(1e-3);
      params_tmp.PLY_percolation.push_back(0.0);

    }
  else if (benchmark_id==4){
      cout << "List cylinders" << endl;

      params_tmp.output_base_name =        "benchmark/output/cylinderList/50L_70icvf";
      params_tmp.cylinders_files.push_back("benchmark/output/cylinderList/50L_70icvf.txt")  ;     /*!< file paths with a list of cilinders obstacles                              */;;

  }
  else if(benchmark_id==5){
    printf("Free diffusion\n");
    params_tmp.output_base_name =        "benchmark/output/free/free";

  }


  params_tmp.traj_file            = params_tmp.output_base_name;                                                

  // Voxel
  string path_voxel = params_tmp.output_base_name + "_voxel.txt";
  loadVoxel(path_voxel, params_tmp);

  benchmark_params = params_tmp;

  return;
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

Benchmark_mcsimulation::~Benchmark_mcsimulation()
{
}
