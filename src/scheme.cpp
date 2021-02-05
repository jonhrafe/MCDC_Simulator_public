#include "scheme.h"
#include <fstream>
#include <iostream>
#include <assert.h>
#include "constants.h"
using namespace std;

Scheme::Scheme()
{
}

Scheme::Scheme(std::string scheme_file_)
{
    scheme_file = scheme_file_;
    readSchemeFile(scheme_file);
}

void Scheme::readSchemeFile(string scheme_file_,bool scale_from_stu)
{
    ifstream in(scheme_file_.c_str());

    if(!in.is_open()){
        std::cout << "Error loading the file" << std::endl;
        in.close();
        return;
    }

    this->scheme_file = scheme_file_;

    string header_;
    in >> header;
    in >> header_;
    header.append(header_);

    std::size_t found = header.find("STEJSKALTANNER");
    std::size_t found_APGSE = header.find("APGSE");
    std::size_t found_waveForm = header.find("WAVEFORM");

    if (found!=std::string::npos){
        type = "PGSE";
    }
    else if (found_APGSE!=std::string::npos){
        type = "APGSE";
    }
    else if (found_waveForm!=std::string::npos){
        type = "WAVEFORM";
    }
    else{
        cout << "Sequence type not valid\n\nAborting" << endl;
        type = "Unknown";
        assert(true);
    }

    if(type == "PGSE"){
        readPGSE(in,scale_from_stu);
    }
    else if (type == "APGSE"){
        readAPGSE(in,scale_from_stu);
    }
    else if (type == "WAVEFORM"){
        readWaveForm(in,scale_from_stu);
    }

    // cout << "[Satus] Scheme file successfully loaded"<< endl;
}

void Scheme::readPGSE(ifstream& in, bool scale_from_stu)
{
    vector<double> scheme_line;
    double tmp;

    num_rep = 0;
    while( in >> tmp){
        scheme_line.push_back(tmp);
        num_rep++;
        for(int i = 0 ; i < 6; i++){
            in >> tmp;
            scheme_line.push_back(tmp);
        }
        scheme.push_back(scheme_line);
        scheme_line.clear();
    }

    if(scale_from_stu)
        for(unsigned i = 0; i < scheme.size(); i++){
            scheme[i][3]/=m_to_mm; // G (T/m) to (T/mm)
            scheme[i][4]*=s_to_ms; // Delta (s) to (ms)
            scheme[i][5]*=s_to_ms; // delta (s) to (ms)
            scheme[i][6]*=s_to_ms; // TE    (s) to (ms)
        }

    in.close();
}

void Scheme::readAPGSE(ifstream &in, bool scale_from_stu)
{
    vector<double> scheme_line;
    double tmp;

    num_rep = 0;
    while( in >> tmp){
        scheme_line.push_back(tmp);
        num_rep++;
        for(int i = 0 ; i < 8; i++){
            in >> tmp;
            scheme_line.push_back(tmp);
        }
        scheme.push_back(scheme_line);
        scheme_line.clear();
    }

    if(scale_from_stu)
        for(unsigned i = 0; i < scheme.size(); i++){
            scheme[i][3]/=m_to_mm; // G1 (T/m) to (T/mm)
            scheme[i][4]/=m_to_mm; // G2 (T/m) to (T/mm)
            scheme[i][5]*=s_to_ms; // Delta  (s) to (ms)
            scheme[i][6]*=s_to_ms; // delta1 (s) to (ms)
            scheme[i][7]*=s_to_ms; // delta2 (s) to (ms)
            scheme[i][8]*=s_to_ms; // TE     (s) to (ms)
        }
    in.close();
}

void Scheme::readWaveForm(ifstream &in, bool scale_from_stu)
{
    this->scale_from_stu = scale_from_stu;
    in.close();
}
