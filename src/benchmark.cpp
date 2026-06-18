#include "benchmark.h"
#include "constants.h"
#include <filesystem> // Add this include for filesystem
#include <chrono>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

using namespace std;

// Helper function to read doubles from a file
std::vector<double> readDoublesFromFile(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }
    std::vector<double> values;
    double value;
    while (file >> value) {
        values.push_back(value);
    }
    return values;
}

// Helper function to calculate Mean Absolute Difference
double calculateMeanAbsoluteDifference(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw std::runtime_error("Vectors have different sizes");
    }
    double sum_diff = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        sum_diff += std::abs(v1[i] - v2[i]);
    }
    return sum_diff / v1.size();
}


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

    params_h1.scheme_file =  "dev/Scheme_3shells3_40ms.txt";
    std::string timestamp = getCurrentDateTime();
    params_h1.output_base_name = output_dir + "/benchmark_test1_single_Hexagonal_" + timestamp;
    params_h1.num_walkers = 2000;
    params_h1.num_steps   = 1000;
    params_h1.sim_duration= 41;
    params_h1.diffusivity = 2.0e-6;
    params_h1.write_bin      = false;
    params_h1.write_txt      = true;
    params_h1.use_mm_ms = false;   // values are programmatic; SI scheme scaling preserved
    params_h1.hex_cyl_packing = true;
    params_h1.hex_packing_radius = 1e-3;
    //params_h1.ini_walker_flag = "extra";
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
    params_h3.num_walkers *= 5;
    params_h3.output_base_name = output_dir + "/benchmark_test2_N-1_"+ timestamp;
    params_h3.num_proc = std::thread::hardware_concurrency()-1;
    params_h3.write_txt      = true;
    ParallelMCSimulation Sim3(params_h3);

    Sim3.startSimulation();

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************   3rd experiment:  *************************/" << SH_DEFAULT << "\n";

    SimErrno::info("Sphere distribution. Single processor",cout);

    Parameters params_h4;

    params_h4.scheme_file =  "dev/Scheme_3shells3_40ms.txt";
    timestamp = getCurrentDateTime();
    params_h4.output_base_name = output_dir + "/benchmark_test3_spheres_" + timestamp;
    params_h4.num_walkers = 10000;
    params_h4.num_steps   = 1000;
    params_h4.sim_duration= 41;
    params_h4.diffusivity = 2.0e-6;
    params_h4.write_txt      = true;
    params_h4.write_bin      = false;
    params_h4.use_mm_ms = false;   // values are programmatic; SI scheme scaling preserved
    params_h4.write_traj    = false;
    params_h4.gamma_sph_packing = true;
    params_h4.gamma_packing_alpha = 5.0;
    params_h4.gamma_packing_beta = 0.5;
    params_h4.gamma_num_obstacles = 100;
    //params_h4.t2_extra = 1000000;
    //params_h4.t2_intra = 1000000;
    params_h4.gamma_icvf = 0.50;
    //params_h4.ini_walker_flag = "intra";
    params_h4.num_proc = 1;
    params_h4.seed = 42;

    ParallelMCSimulation Sim4(params_h4);

    Sim4.startSimulation();
    /*
     * 2.01 separation and n-1 processor:
    */

    cout << endl << endl;
    cout << SH_FG_PURPLE << "/********************  4th experiment:  *************************/" << SH_DEFAULT << "\n";

    SimErrno::info("PLY Mesh File",cout);

    Parameters params_h5;

    params_h5.scheme_file =  "dev/Scheme_3shells3_40ms.txt";
    params_h5.output_base_name = output_dir + "/benchmark_single_PLY_" + timestamp;
    params_h5.num_walkers  = 1000;
    params_h5.num_steps    = 1000;
    params_h5.sim_duration = 41;
    params_h5.diffusivity  = 0.800e-6;
    params_h5.diff_intra   = 0.800e-6;
    params_h5.diff_extra   = 0.600e-6;
    params_h5.write_bin      = false;
    params_h5.write_txt      = true;
    params_h5.use_mm_ms = false;   // values are programmatic; SI scheme scaling preserved
    params_h5.write_traj     = true;
    params_h5.PLY_files.push_back("instructions/meshes/hexagonal_packed_spheres_INT.ply");
    params_h5.PLY_scales.push_back(1e-3);
    params_h5.PLY_percolation.push_back(0.0);
    params_h5.num_proc = 1;std::thread::hardware_concurrency()-1;
    //params_h5.ini_delta_pos = {0,0,0};
    params_h5.ini_walker_flag = "intra";
    //params_h5.obstacle_permeability = 0.0001;
    //params_h5.discard_illegals = false;
    std::pair<Eigen::Vector3d,Eigen::Vector3d> voxel(Eigen::Vector3d(-6e-3,-6e-3,-6e-3),Eigen::Vector3d(6e-3,6e-3,6e-3));
    params_h5.voxels_list.push_back(voxel);
    ParallelMCSimulation Sim5(params_h5);

    Sim5.startSimulation();
    
    std::vector<std::string> expected_files = {
        "dev/benchmark_test1_single_Hexagonal_DWI.txt",
        "dev/benchmark_test2_N-1_DWI.txt",
        "dev/benchmark_test3_spheres_DWI.txt",
        "dev/benchmark_test4_single_PLY_DWI.txt"
    };

    // Corresponding output files
    std::vector<std::string> output_files = {
        params_h1.output_base_name + "_DWI.txt",
        params_h3.output_base_name + "_DWI.txt",
        params_h4.output_base_name + "_DWI.txt",
        params_h5.output_base_name + "_DWI.txt"
    };

    for (size_t i = 0; i < expected_files.size(); ++i) {
        try {
            std::vector<double> expected = readDoublesFromFile(expected_files[i]);
            std::vector<double> output = readDoublesFromFile(output_files[i]);

            // normalize the vectors by it's maximum value
            double max_expected = *std::max_element(expected.begin(), expected.end());
            double max_output = *std::max_element(output.begin(), output.end());
            for (size_t j = 0; j < expected.size(); ++j) {
                expected[j] /= max_expected;
                output[j] /= max_output;
            }
            double mad = calculateMeanAbsoluteDifference(expected, output);

            std::cout << "Experiment " << i + 1 << ": Mean Absolute Difference = " << mad << "\n";
            std::cout << "Passed: " << (mad < 0.1 ? "\033[32myes\033[0m" : "\033[31mno\033[0m") << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Error in experiment " << i + 1 << ": " << e.what() << "\n";
        }
    }

}
