#include "benchmark.h"
#include "constants.h"
#include <filesystem> // Add this include for filesystem
#include <chrono>

using namespace std;


std::string getCurrentDateTime() {
    // Get the current time
    std::time_t now_time = std::time(nullptr);

    // Format the time into a string
    char buffer[20];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d_%H-%M-%S", std::localtime(&now_time));

    return std::string(buffer);
}

Benchmark::Benchmark(string output)
{
    output_dir = output;
}

void Benchmark::start()
{

    // Benchmark for 2 different packing densities and multi-processors.
    // Hexagonal packing only, extra cellular and a PLY mesh.
    // 1 processor and N-1 processors, where N is the total number of processors in your system.
    /*
     * 2.1 separation and 1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************   First experiment:  *************************/" << SH_DEFAULT << "\n";

    SimErrno::info("Hexagonal Packing: 1 um radius and 2.1 um separation. Single processor",cout);

    Parameters params_h1;

    //create this folder: params_h1.output_base_name = output_dir + "/test1/
    
    params_h1.scheme_file = output_dir + "/benchmark.scheme";
    std::string timestamp = getCurrentDateTime();
    params_h1.output_base_name = output_dir + "/benchmark_test1_single_Hexagonal_" + timestamp;
    params_h1.num_walkers = 10000;
    params_h1.num_steps   = 1000;
    params_h1.sim_duration= 36;
    params_h1.diffusivity = 2.0e-6;
    params_h1.write_txt      = true;
    params_h1.scale_from_stu = true;
    params_h1.hex_cyl_packing = true;
    params_h1.hex_packing_radius = 1e-3;
    params_h1.ini_walker_flag = "extra";
    params_h1.hex_packing_icvf = 0.6;
    params_h1.num_proc = 1;
    params_h1.seed = 42;

    ParallelMCSimulation Sim1(params_h1);

    Sim1.startSimulation();

    /*
     * 2.01 separation and n-1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************  2nd experiment:  *************************/" << SH_DEFAULT << "\n";

    SimErrno::info("Hexagonal Packing: 1 um radius and 2.1 um separation. N-1 processors",cout);

    Parameters params_h3 = params_h1;
    params_h3.num_walkers *= 10;
    params_h3.output_base_name = output_dir + "/benchmark_test2_N-1_"+ timestamp;
    params_h3.hex_packing_separation = 2.1e-3;
    params_h3.num_proc = std::thread::hardware_concurrency()-1;
    params_h3.write_txt      = true;
    ParallelMCSimulation Sim3(params_h3);

    Sim3.startSimulation();

    /*
     * 2.01 separation and n-1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************  3th experiment:  *************************/" << SH_DEFAULT << "\n";

    SimErrno::info("PLY Mesh File",cout);

    Parameters params_h5;

    params_h5.scheme_file = output_dir + "/benchmark.scheme";
    params_h5.output_base_name = output_dir + "/benchmark_single_PLY" + timestamp;
    params_h5.num_walkers = 10000;
    params_h5.num_steps   = 1000;
    params_h5.sim_duration= 36;
    params_h5.diffusivity = 0.8e-6;
    params_h5.write_txt      = true;
    params_h5.scale_from_stu = true;
    params_h5.PLY_files.push_back("instructions/meshes/hexagonal_packed_spheres.ply");
    params_h5.PLY_scales.push_back(1e-3);
    params_h5.PLY_percolation.push_back(0.0);
    params_h5.num_proc = std::thread::hardware_concurrency()-1;
    params_h5.ini_delta_pos = {0,0,0};
    params_h5.ini_walker_flag = "delta";

    std::pair<Eigen::Vector3d,Eigen::Vector3d> voxel(Eigen::Vector3d(-1e-3,-1e-3,-2e-3),Eigen::Vector3d(1e-3,1e-3,2e-3));

    params_h5.voxels_list.push_back(voxel);
    ParallelMCSimulation Sim5(params_h5);

    Sim5.startSimulation();


}


