// o   o   o-o     oo-o     o-o      o-o               o      o
// |\ /|  /       / |  \   /        |     o            |      |
// | O | O       o  |   O O          o-o    o-O-o o  o |  oo -o- o-o o-o
// |   |  \     /   |  /   \            | | | | | |  | | | |  |  | | |
// o   o   o-o o    o-o     o-o     o--o  | o o o o--o o o-o- o  o-o o

#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "dynamicsSimulation.h"
#include <fstream>
#include "pgsesequence.h"
#include "parallelmcsimulation.h"
#include "voxel.h"
#include "cylinder.h"
#include <thread>
#include "simerrno.h"

typedef unsigned int uint;

using namespace std;
using namespace Eigen;

void printUsage();

int main(int argn, char* argv[])
{

    string conf = "";

    if(argn == 2){
        conf = argv[1];
    }
    else{
        printUsage();
        return -1;
    }


    ParallelMCSimulation simulation(conf);

    simulation.startSimulation();

    #ifdef __linux__

        string command = "notify-send -i emblem-default \"MC/DC\" \"Simulation Finished\"";
        system (command.c_str());
    #endif

    return 0;

}

void printUsage(){

    cout << " ███╗   ███╗ ██████╗    ██╗██████╗  ██████╗" << endl;
    cout << " ████╗ ████║██╔════╝   ██╔╝██╔══██╗██╔════╝" << endl;
    cout << " ██╔████╔██║██║       ██╔╝ ██║  ██║██║     " << endl;
    cout << " ██║╚██╔╝██║██║      ██╔╝  ██║  ██║██║     " << endl;
    cout << " ██║ ╚═╝ ██║╚██████╗██╔╝   ██████╔╝╚██████╗" << endl;
    cout << " ╚═╝     ╚═╝ ╚═════╝╚═╝    ╚═════╝  ╚═════╝" << endl;

    cout << endl; cout << endl;
    cout << " Version: " << VERSION_ID << endl;
    cout << " Usage: MC-DC_Simulator <configuration_file.conf>\n\n";
    cout << " <configuration_file.conf>  Plain .txt file with the simulation parameters (see https://github.com/jonhrafe/MCDC_Simulator_public):\n\n";

    cout << "   N <int>                      Number of particles.\n";
    cout << "   T <int>                      Number of time steps.\n";
    cout << "   duration <float>             Diffusion duration in seconds.\n";
    cout << "   out_file_index <string>      Simulation ouput path and prefix.\n";
    cout << "   scheme_file <string>         Simulation protocol.\n";
    cout << "   scale_from_stu <int>         not 0 if the protocol is in SU.\n";

    cout << "   write_txt <int>              not 0 for .txt ouput.\n";
    cout << "   write_bin <int>              not 0 for .bin ouput.\n";
    cout << "   write_traj_file              not 0 to write the trajfile.\n";

    cout << "   <obstacle>                   Obstacle definition tag.\n";
    cout << "   <cylinder_gamma_packing>     Gamma cylinders obstacles tag.\n";
    cout << "   <ply_obstacle>               ply-mesh-model obstacle tag.\n";
    cout << "   ini_walkers_pos <string>     Custom initial particles position (intra, extra).\n";
    cout << "   num_process <int>            Number of processors to use.\n";

    cout << "   <END>                        END of the conf-file parameters (needed).\n";
}
