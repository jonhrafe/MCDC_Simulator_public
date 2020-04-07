//███╗   ███╗ ██████╗    ██╗██████╗  ██████╗    ███████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗ ██████╗ ██████╗
//████╗ ████║██╔════╝   ██╔╝██╔══██╗██╔════╝    ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██╔═══██╗██╔══██╗
//██╔████╔██║██║       ██╔╝ ██║  ██║██║         ███████╗██║██╔████╔██║██║   ██║██║     ███████║   ██║   ██║   ██║██████╔╝
//██║╚██╔╝██║██║      ██╔╝  ██║  ██║██║         ╚════██║██║██║╚██╔╝██║██║   ██║██║     ██╔══██║   ██║   ██║   ██║██╔══██╗
//██║ ╚═╝ ██║╚██████╗██╔╝   ██████╔╝╚██████╗    ███████║██║██║ ╚═╝ ██║╚██████╔╝███████╗██║  ██║   ██║   ╚██████╔╝██║  ██║
//╚═╝     ╚═╝ ╚═════╝╚═╝    ╚═════╝  ╚═════╝    ╚══════╝╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝

#include <iostream>
#include "Eigen/Core"
#include "Eigen/Dense"
#include <dynamicsSimulation.h>
#include <fstream>
#include <pgsesequence.h>
#include <parallelmcsimulation.h>
#include "voxel.h"
#include "cylinder.h"
#include <thread>
#include <simerrno.h>
#include "win_types.h"  // quick hack for compilitation on Windows

using namespace std;
using namespace Eigen;

void printUsage();

int main(int argn, char* argv[])
{

    string conf = "/home/jonathan/Desktop/output/Percolacion/delta_center/conf_first_delta_perc.conf";


    if(argn == 2){
        conf = argv[1];
    }
    if(argn != 2){
        printUsage();
        return -1;
    }

    ParallelMCSimulation simulation_short_scheme(conf);

    simulation_short_scheme.startSimulation();

    #ifdef __linux__

        string command = "notify-send -i emblem-default \"MC/DC\" \"Simulation Finished\"";
        system (command.c_str());
    #endif

    return 0;

}

void printUsage(){

    cout << "// ███╗   ███╗ ██████╗    ██╗██████╗  ██████╗    ███████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗ ██████╗ ██████╗  " << endl;
    cout << "// ████╗ ████║██╔════╝   ██╔╝██╔══██╗██╔════╝    ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██╔═══██╗██╔══██╗ " << endl;
    cout << "// ██╔████╔██║██║       ██╔╝ ██║  ██║██║         ███████╗██║██╔████╔██║██║   ██║██║     ███████║   ██║   ██║   ██║██████╔╝ " << endl;
    cout << "// ██║╚██╔╝██║██║      ██╔╝  ██║  ██║██║         ╚════██║██║██║╚██╔╝██║██║   ██║██║     ██╔══██║   ██║   ██║   ██║██╔══██╗ " << endl;
    cout << "// ██║ ╚═╝ ██║╚██████╗██╔╝   ██████╔╝╚██████╗    ███████║██║██║ ╚═╝ ██║╚██████╔╝███████╗██║  ██║   ██║   ╚██████╔╝██║  ██║ " << endl;
    cout << "// ╚═╝     ╚═╝ ╚═════╝╚═╝    ╚═════╝  ╚═════╝    ╚══════╝╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝ " << endl;


    cout << endl; cout << endl;

    cout << " Usage: MC-DC_Simulator <configuration_file.conf>\n\n";

    cout << " <configuration_file.conf>   Plain .txt file, may contain the following parameters (see doc for examples):\n";

    cout << "  N <int>                      Number of particles.\n";
    cout << "  T <int>                      Number of time steps.\n";
    cout << "  duration <float>             Diffusion duration in seconds.\n";
    cout << "  out_file_index <string>      Simulation ouput path and prefix.\n";
    cout << "  scheme_file <string>         Simulation protocol.\n";
    cout << "  scale_from_stu <int>         not 0 if the protocol is in SU.\n";

    cout << "  write_txt <int>              not 0 for .txt ouput.\n";
    cout << "  write_bin <int>              not 0 for .bin ouput.\n";
    cout << "  write_traj_file              not 0 to write the trajfile.\n";

    cout << "  <obstacle>                   Obstacle definition tag.\n";
    cout << "  <cylinder_gamma_packing>     Gamma cylinders obstacles tag.\n";
    cout << "  <ply_obstacle>               ply-mesh-model obstacle tag.\n";
    cout << "  ini_walkers_pos <string>     Custom initial particles position (intra, extra).\n";
    cout << "  num_process <int>            Number of processors to use.\n";

    cout << "  <END>                        END of the conf-file parameters (needed).\n";




//    ifstream in("../misc/Params-ERROR.txt");
//    if(!in){
//        cout << "Error in the local Paths "<< endl;
//        return;
//     }
//    string temp;
//    while( std::getline(in, temp)){
//        cout << temp << endl;
//    }
//    in.close();
}
