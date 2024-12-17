//███╗   ███╗ ██████╗    ██╗██████╗  ██████╗    ███████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗ ██████╗ ██████╗
//████╗ ████║██╔════╝   ██╔╝██╔══██╗██╔════╝    ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██╔═══██╗██╔══██╗
//██╔████╔██║██║       ██╔╝ ██║  ██║██║         ███████╗██║██╔████╔██║██║   ██║██║     ███████║   ██║   ██║   ██║██████╔╝
//██║╚██╔╝██║██║      ██╔╝  ██║  ██║██║         ╚════██║██║██║╚██╔╝██║██║   ██║██║     ██╔══██║   ██║   ██║   ██║██╔══██╗
//██║ ╚═╝ ██║╚██████╗██╔╝   ██████╔╝╚██████╗    ███████║██║██║ ╚═╝ ██║╚██████╔╝███████╗██║  ██║   ██║   ╚██████╔╝██║  ██║
//╚═╝     ╚═╝ ╚═════╝╚═╝    ╚═════╝  ╚═════╝    ╚══════╝╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝  ╚═╝

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
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
void printUsagedataSynth();
string secondsToMinutes(double t);
int main(int argn, char* argv[])
{

    string scheme_file = "";
    string trajfile = "";

    if(argn != 3){
        printUsagedataSynth();
        return -1;
    }

    if(argn == 2 || argn == 3){
        scheme_file = argv[2];
        trajfile = argv[1];
    }


    PGSESequence sequence(scheme_file.c_str(),trajfile.c_str());

    for (int i = 0 ; i < sequence.num_rep; i++){

        sequence.scheme[i][3]*=(1.0/1000.0); // Teslas/mimilimeters
        sequence.scheme[i][4]*=(1000.0);
        sequence.scheme[i][5]*=(1000.0);
        sequence.scheme[i][6]*=(1000.0);
    }



    time_t start,now;                                /*!< Auxiliar Variable for time recording and estimation for time*/

    time(&start);
    sequence.getDWISignal();
    time(&now);
    double seconds_passed = difftime(now,start);
    SimErrno::info(" Signal computed in: " + secondsToMinutes(seconds_passed),cout,false);

    for(int i = 0 ; i < sequence.num_rep; i++)
        cout << sequence.DWI[i] << endl;

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

void printUsagedataSynth(){
    cout << "dataSynth <trajfile> <scheme_file>" << endl;
}

string secondsToMinutes(double t)
{
    if(t < 60){
        return std::to_string(int(t)) + " seconds";
    }

    int mins    = int(t/60.0);

    int seconds = int(t - mins*60);

    string min_s =  (mins>1)?" minutes":" minute";
    return std::to_string(mins) +  min_s +  " and " + std::to_string(seconds) + " seconds";

}
