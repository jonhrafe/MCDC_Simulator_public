#include "benchmark.h"
#include "constants.h"
using namespace std;

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

    params_h1.scheme_file = output_dir + "/benchmark.scheme";
    params_h1.output_base_name = output_dir + "/test1/benchmark_test1_single_2.1";
    params_h1.num_walkers = 10000;
    params_h1.num_steps   = 1000;
    params_h1.sim_duration= 36;
    params_h1.diffusivity = 2.0e-6;
    params_h1.write_txt      = true;
    params_h1.scale_from_stu = true;
    params_h1.hex_cyl_packing = true;
    params_h1.hex_packing_radius = 1e-3;
    params_h1.ini_walker_flag = "extra";

    params_h1.hex_packing_separation = 2.1e-3;
    params_h1.num_proc = 1;

    ParallelMCSimulation Sim1(params_h1);

    Sim1.startSimulation();


    /*
     * 2.1 separation and 2 processors:
    */
    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************   Second experiment:  *************************/" << SH_DEFAULT << "\n";


    SimErrno::info("Hexagonal Packing: 1 um radius and 2.01 um separation. Single processor",cout);

    Parameters params_h2 = params_h1;
    params_h2.output_base_name = output_dir + "/test2/benchmark_test2_double_2.01";
    params_h2.num_proc = 1;
    params_h2.hex_packing_separation = 2.01e-3;

    ParallelMCSimulation Sim2(params_h2);

    Sim2.startSimulation();


    /*
     * 2.01 separation and n-1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************  3rd experiment:  *************************/" << SH_DEFAULT << "\n";



    SimErrno::info("Hexagonal Packing: 1 um radius and 2.1 um separation. N-1 processors",cout);

    Parameters params_h3 = params_h1;
    params_h3.num_walkers *= 10;
    params_h3.output_base_name = output_dir + "/test3/benchmark_test3_N-1_2.1";
    params_h3.hex_packing_separation = 2.1e-3;
    params_h3.num_proc = std::thread::hardware_concurrency()-1;
    params_h3.write_txt      = false;
    ParallelMCSimulation Sim3(params_h3);

    Sim3.startSimulation();



    /*
     * 2.01 separation and n-1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************  4th experiment:  *************************/" << SH_DEFAULT << "\n";



    SimErrno::info("Hexagonal Packing: 1 um radius and 2.01 um separation. N-1 processors",cout);

    Parameters params_h4 = params_h1;
    params_h4.num_walkers *= 10;
    params_h4.output_base_name = output_dir + "/test4/benchmark_test4_N-1_2.01";
    params_h4.hex_packing_separation = 2.01e-3;
    params_h4.num_proc = std::thread::hardware_concurrency()-1;
    params_h4.write_txt      = false;
    ParallelMCSimulation Sim4(params_h4);

    Sim4.startSimulation();


    /*
     * 2.01 separation and n-1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************  5th experiment:  *************************/" << SH_DEFAULT << "\n";



    SimErrno::info("PLY Mesh File",cout);



    Parameters params_h5;

    params_h5.scheme_file = output_dir + "/benchmark.scheme";
    params_h5.output_base_name = output_dir + "/test5/benchmark_test5_single_PLY";
    params_h5.num_walkers = 10000;
    params_h5.num_steps   = 1000;
    params_h5.sim_duration= 36;
    params_h5.diffusivity = 0.8e-6;
    params_h5.write_txt      = true;
    params_h5.scale_from_stu = true;
    params_h5.PLY_files.push_back(output_dir+"/Cylinder_0.5radius_15_units_long.ply");
    params_h5.num_proc = std::thread::hardware_concurrency()-1;
    params_h5.ini_delta_pos = {0,0,0};
    params_h5.ini_walker_flag = "delta";

    std::pair<Eigen::Vector3d,Eigen::Vector3d> voxel(Eigen::Vector3d(-1e-3,-1e-3,-2e-3),Eigen::Vector3d(1e-3,1e-3,2e-3));

    params_h5.voxels_list.push_back(voxel);
    ParallelMCSimulation Sim5(params_h5);

    Sim5.startSimulation();


}


