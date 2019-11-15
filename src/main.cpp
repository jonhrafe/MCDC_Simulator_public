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

using namespace std;
using namespace Eigen;

void printUsage();

int main(int argn, char* argv[])
{

    string conf = "/home/jonathan/Desktop/output/Percolacion/delta_center/conf_first_delta_perc.conf";


    if(argn == 2 || argn == 3){
        conf = argv[1];
    }

    if(conf.compare("benchmark") == 0){
        if(argn!=3){
            printUsage();
            return -1;
        }
    }

    if(argn != 2){
        printUsage();
        return -1;
    }

    ParallelMCSimulation simulation_short_scheme(conf);

    simulation_short_scheme.startSimulation();

    #ifdef __linux__

        string command = "notify-send -i emblem-default \"MC/DC\" \"Simulation Finished\"";
        int out = system (command.c_str());
    #endif

    return 0;

}

void printUsage(){
    ifstream in("../misc/Params-ERROR.txt");
    if(!in){
        cout << "Error in the local Paths "<< endl;
        return;
     }
    string temp;
    while( std::getline(in, temp)){
        cout << temp << endl;
    }
    in.close();
}
