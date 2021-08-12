#include "simulablesequence.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "constants.h"
using namespace std;



void SimulableSequence::writeResultingData(std::string output_base_name)
{
    string outDWI   = output_base_name + "_DWI.txt";
    string outDWIi  = output_base_name + "_DWI_img.txt";

    ofstream dwi_out(outDWI,std::ofstream::out);
    ofstream dwii_out(outDWIi,std::ofstream::out);

    for (unsigned i = 0 ; i < this->DWI.size(); i++ ){
        dwi_out  << this->DWI[i]  << endl;
        dwii_out << this->DWIi[i] << endl;
    }

    dwi_out.close();
    dwii_out.close();

    if(this->subdivision_flag){
        string out_sub_DWI   = output_base_name + "_voxels_DWI.txt";
        string out_sub_DWIi  = output_base_name + "_voxels_DWI_img.txt";

        ofstream sub_dwi_out (out_sub_DWI,std::ofstream::out);
        ofstream sub_dwii_out(out_sub_DWIi,std::ofstream::out);

        //## HEADER
        // #Num_voxels
        sub_dwi_out  << subdivisions.size() << endl;
        sub_dwii_out << subdivisions.size() << endl;
        // #Size_DWI
        sub_dwi_out  << this->DWI.size() << endl;
        sub_dwii_out << this->DWI.size() << endl;

        for(uint s = 0; s < this->subdivisions.size(); s++ ){
            for (uint i = 0 ; i < this->DWI.size(); i++ ){
                sub_dwi_out  << this->sub_DWI[s][i]  << endl;
                sub_dwii_out << this->sub_DWIi[s][i] << endl;
            }
        }
        sub_dwi_out.close();
        sub_dwii_out.close();
    }//END subdivision_flag

}

void SimulableSequence::writePhaseShiftDistribution(string output_base_name)
{

    string outPhase = output_base_name + "_phase_shift.txt";

    ofstream phase_out(outPhase,std::ofstream::out);   //phase shift


    for (unsigned i = 0 ; i < this->DWI.size(); i++ ){
        //write the full histogram for the gradient i
        for(int h = 0; h < 3600; h++){
            phase_out << this->phase_shift_distribution(i,h);
            if(h < 3600-1)
                phase_out << " ";
            else
                phase_out << endl;
        }
    }

    phase_out.close();
}

void SimulableSequence::cleanPhaseShift()
{
    for(uint i = 0 ; i < DWI.size();i++){
        this->phase_shift[i] = 0;
    }
}

void SimulableSequence::cleanDWISignal(){
    for(uint s=0; s< num_rep; s++){
        DWI[s] = 0; // Real part
        DWIi[s]= 0; // Img part
        phase_shift[s] = 0;
    }
}

void SimulableSequence::initializeSubdivisionSignals(){

    for(uint s = 0; s < subdivisions.size(); s++){
        vector<double> tmp_DWI(num_rep,0);
        vector<double> tmp_DWIi(num_rep,0);
        this->sub_DWI.push_back(tmp_DWI);
        this->sub_DWIi.push_back(tmp_DWIi);

        if(separate_signal){
            vector<double> tmp_DWI_intra(num_rep,0);
            vector<double> tmp_DWI_extra(num_rep,0);
            this->sub_DWI_intra.push_back(tmp_DWI_intra);
            this->sub_DWI_extra.push_back(tmp_DWI_extra);
        }
    }
}

void SimulableSequence::initializeIntraExtraSignals()
{
    for(auto i =0; i < num_rep; i++){
        this->DWI_intra.push_back(0);
        this->DWI_extra.push_back(0);
    }
}
