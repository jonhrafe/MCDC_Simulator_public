// o   o   o-o     oo-o     o-o      o-o               o      o
// |\ /|  /       / |  \   /        |     o            |      |
// | O | O       o  |   O O          o-o    o-O-o o  o |  oo -o- o-o o-o
// |   |  \     /   |  /   \            | | | | | |  | | | |  |  | | |
// o   o   o-o o    o-o     o-o     o--o  | o o o o--o o o-o- o  o-o o

#include <iostream>
#include <thread>
#include <fstream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "dynamicsSimulation.h"
#include "pgsesequence.h"
#include "parallelmcsimulation.h"
#include "voxel.h"
#include "cylinder.h"
#include "simerrno.h"
#include "benchmark.h"

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

void printUsage();

int main(int argn, char* argv[])
{

    string conf = "";
    string output_benchmark = "";

    if(argn == 2){
        conf = argv[1];
        ParallelMCSimulation simulation(conf);
        simulation.startSimulation();
    }
    else if(argn == 3){
        string secondParam = argv[1];
        if (secondParam == "--conf") {
            conf = argv[2];

            ParallelMCSimulation simulation(conf);

            simulation.startSimulation();

        } else if (secondParam == "--benchmark") {
            output_benchmark = argv[2];

            Benchmark bench(output_benchmark);
            bench.start();
        } else {
            printUsage();
            return -1;
        }
    } else {
        printUsage();
        return -1;
    }

    #ifdef __linux__

        string command = "notify-send -i emblem-default \"MC/DC\" \"Simulation Finished\"";
        system (command.c_str());
    #endif

    return 0;

}
void printUsage() {
    cout << R"(
   *                    (                           
 (  `      (            )\ )    (        )       )  
 )\))(     )\          (()/(    )\    ( /(    ( /(  
((_)()\  (((_)         /(_)) (((_)   )(_))   )\()) 
(_()((_) )\___     __  (_))_  )\___  ((_)    ((_)\  
|  \/  |((/ __|   / /  |   \((/ __| |_  )   / (_) 
| |\/| | | (__   / /   | |) || (__   / /  _| () |  
|_|  |_|  \___| /_/    |___/  \___| /___|(_)\__/

 MC-DC Simulator - Advanced Diffusion MRI Simulation
 Version: )" << VERSION_ID << R"(

 Usage:
   MC-DC_Simulator <configuration_file.conf>
   MC-DC_Simulator --conf <configuration_file.conf>
   MC-DC_Simulator --benchmark <output_file>

 Parameters (configuration file):
   Simulation Parameters:
   
     N <int>                 Number of particles in the simulation.
     T <int>                 Number of time steps.
     duration <float>        Diffusion duration (seconds).

   File Outputs:
     out_file_index <string> Output path and prefix for simulation results.
     write_txt <int>         1 for .txt output, 0 to disable.
     write_bin <int>         1 for .bin output, 0 to disable.
     write_traj_file <int>   1 to output trajectory files, 0 to disable.

   Obstacle and Protocol:
     scheme_file <string>    Path to the simulation protocol file.
     scale_from_stu <int>    Set to 1 if protocol uses scaled units, 0 otherwise.
     <obstacle>              Define obstacle configurations (e.g., sphere, cylinder).
     <cylinder_gamma_packing> Specify gamma cylinder obstacle settings.
     <ply_obstacle>          Define PLY mesh model obstacles.
     ini_walkers_pos <string> Custom initial positions for particles (e.g., intra, extra).

   Performance and Process Management:
     num_process <int>       Number of processors to use (for parallelism).

   Other:
     <END>                   Marks the end of the configuration file (mandatory).

 Examples:
   MC-DC_Simulator simulation.conf
   MC-DC_Simulator --conf simulation.conf
   MC-DC_Simulator --benchmark benchmark_results.txt

 For detailed instructions and examples, see:
   https://github.com/jonhrafe/MCDC_Simulator_public

)";
}