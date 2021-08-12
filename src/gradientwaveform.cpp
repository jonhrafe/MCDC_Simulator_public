//!  GradientWaveform Sequence Class  =============================================================/
/*!
  Derived Class.
  Implementation of the the General Wavefroms

  \date   March 2018
  \author Jonathan Rafael
  \version 1.42.0
*=====================================================================================*/

#include "gradientwaveform.h"
#include "constants.h"
#include "Eigen/Dense"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <assert.h>

using namespace std;



/*! \class  GradientWaveform
 *  \brief  Implementation of the the General Wavefroms
 */

GradientWaveform::GradientWaveform()
{
    wave_duration  = 0;
    dt = 0;
    wave_bins = 0;
    separate_signal=false;
    num_rep=0;
}

GradientWaveform::GradientWaveform(Scheme &scheme_)
{
    this->num_rep=0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    readSchemeParameters(scheme_);
    //phase_shift_distribution.resize(num_rep,3600);
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
    separate_signal=false;

}

//TODO correjir esto
GradientWaveform::GradientWaveform(Scheme& scheme_, const char *traj_file_name)
{
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    readSchemeParameters(scheme_);
    trajectory.setTrajFile(traj_file_name);
    T = uint(trajectory.T);
    dyn_duration = trajectory.dyn_duration;
    //phase_shift_distribution.resize(scheme_.num_rep,3600); // Dynamic duration no es igual que el waveform duration
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
    separate_signal=false;
}

double GradientWaveform::getNumericalbValue(unsigned)
{
    return 0;
}


void GradientWaveform::readSchemeParameters(Scheme &scheme_)
{
    scheme_file = scheme_.scheme_file;
    this->scale_from_stu = scheme_.scale_from_stu;
    readSchemeFile();

    for(unsigned i = 0 ; i < uint(num_rep); i++){
        DWI.push_back(0);
        if(this->img_signal == true)
            DWIi.push_back(0);
        phase_shift.push_back(0);
    }
}

void GradientWaveform::readSchemeFile()
{
    ifstream in(this->scheme_file.c_str());

    //Reads the 2 string header
    string header_;
    in >> header_;
    in >> header_;

    in >> this->wave_duration;
    in >> this->wave_bins;
    float holder;
    in >> holder;
    this->num_rep = int(holder);

    if(scale_from_stu){
        wave_duration*= s_to_ms;
    }
    this->dt = wave_duration/(wave_bins-1);

    for (uint i = 0 ; i < uint(num_rep); i++)
    {
        for(uint t = 0; t< uint(this->wave_bins); t++){

            vector<float> wave_vector={0,0,0};
            in >> wave_vector[0];
            in >> wave_vector[1];
            in >> wave_vector[2];

            if(scale_from_stu){
                wave_vector[0]/=m_to_mm; // Gx (T/m) to (T/mm);
                wave_vector[1]/=m_to_mm; // Gy (T/m) to (T/mm);
                wave_vector[2]/=m_to_mm; // Gz (T/m) to (T/mm));
            }

            this->waveform.push_back(wave_vector);
        }
    }
    in.close();
}

void GradientWaveform::getGradImpulse(int i, double t, double tLast, Eigen::Vector3d& Gdt){
    return;
}

void GradientWaveform::getInterpolatedGradImpulse(uint rep_index, double t_sim, double t_sim_last, Eigen::Vector3d& Gdt)
{

    // If the simulation time is bigger than the waveform we fill with zeros.
    if(t_sim > wave_duration){
        Gdt = {0.,0.,0.};
        return;
    }

    //The index in the Waveform time resolution.
    // index of the waveform duration time (initial one) (used for interpolation)
    uint wt_1 = uint(t_sim/this->dt);

    // index of the wave direction and maginitud
    uint index = rep_index*uint(wave_bins) + wt_1;

    uint index_next = (wt_1+1>= uint(this->wave_bins))?index:index+1;

    Eigen::Vector3d iniG(waveform[index][0],waveform[index][1],waveform[index][2]);
    Eigen::Vector3d nextG(waveform[index_next][0],waveform[index_next ][1],waveform[index_next ][2]);


    // percentaje of iniG  = (1 - (dt_sim - t_1*dt)/dt) = ini_perc;
    // percentaje of nextG = (1.0-ini_perc) = sec_perc;

    double ini_perc = (1.0) - (t_sim - wt_1*dt)/(this->dt);

    //Linear interpolation scaled by the applied time
    Gdt = (ini_perc*iniG + (1.0-ini_perc)*nextG)*(t_sim-t_sim_last);

//    if(iniG[0]!= 0){
//        cout << iniG     << endl;
//        cout << ini_perc << endl;
//        cout << nextG    << endl;
//        cout << Gdt    << endl;
//    }
}

void GradientWaveform::update_phase_shift(double time_step, Eigen::Matrix3Xd trajectory)
{
    Eigen::Vector3d xt;
    Eigen::Vector3d Gdt;
    double dt,dt_last;

    for (uint t=1; t <= uint(this->T) ;t++){ //TODO: checar si deberia ser <= T
        //Displacement
        xt[0] = trajectory(0,t) - trajectory(0,0);
        xt[1] = trajectory(1,t) - trajectory(1,0);
        xt[2] = trajectory(2,t) - trajectory(2,0);

        double dos_pi = 2.0*M_PI;

        dt      = time_step*(t);
        dt_last = time_step*(t-1);
        for(uint s=0; s < uint(num_rep) ;s++){

            getInterpolatedGradImpulse(s,dt,dt_last,Gdt);
            double val = giro*(Gdt[0]*xt[0]+Gdt[1]*xt[1]+Gdt[2]*xt[2]);
            val = fmod(val,dos_pi);
            phase_shift[s] = fmod(phase_shift[s] + val,dos_pi);
        }
    }
}

void GradientWaveform::update_phase_shift(double dt, double dt_last, Walker walker)
{
    // Deprecated for General Forms.
    return;
}


void GradientWaveform::getDWISignal()
{
    trajectory.initTrajReaderFile();

    trajectory.readTrajectoryHeader();

    double N        = trajectory.N;
    double T        = trajectory.T;
    double duration = trajectory.dyn_duration;
    double rt       = duration/T;
    double dos_pi   = 2.0*M_PI;
    double dt,xt[3],dt_last;


    Eigen::Matrix3Xd steps_log; // complete trajectory of one walker

    Eigen::Vector3d Gdt;
    Eigen::VectorXd phase_shift;

    steps_log.resize(3,unsigned(T+1));
    phase_shift.resize(num_rep);

    for (int w = 0; w < N; w++)
    {
        trajectory.readCurrentWalkersTrajectory(steps_log);
        for (uint t = 1; t <= uint(trajectory.T); t++)
        {
            dt = rt*(t);
            dt_last = rt*(t-1.0);

            xt[0] = steps_log(0,t) - steps_log(0,0);
            xt[1] = steps_log(1,t) - steps_log(1,0);
            xt[2] = steps_log(2,t) - steps_log(2,0);

            for(uint s=0; s < num_rep ;s++)
            {
                getInterpolatedGradImpulse(s,dt,dt_last,Gdt);
                double val = giro*(Gdt[0]*xt[0]+Gdt[1]*xt[1]+Gdt[2]*xt[2]);

                val = fmod(val,2*M_PI);
                //printf("%d - %1.25f \n",w,val );
                phase_shift[s] = fmod(phase_shift[s]+ val,dos_pi);
            }
        }

        for(uint s=0; s < num_rep; s++){
            DWI[s] += cos(phase_shift[s]); // Real part
            if(this->img_signal == true)
                DWIi[s]+= sin(phase_shift[s]); // Img part

            phase_shift[s] = 0;
        }
    }
}



void GradientWaveform::update_DWI_signal(Walker& walker)
{
    for(uint s=0; s< uint(num_rep); s++){

        double cos_phase_shift = cos(phase_shift[s]);
        double sin_phase_shift = sin(phase_shift[s]);

        DWI[s] += cos_phase_shift; // Real part
        if(this->img_signal == true)
            DWIi[s]+= sin_phase_shift; // Img part

        if(this->separate_signal){

            if(walker.location == Walker::RelativeLocation::intra){
                DWI_intra[s]+=cos_phase_shift;
            }
            else if(walker.location == Walker::RelativeLocation::extra){
                DWI_extra[s]+=cos_phase_shift;
            }
            else{
                cout << walker.location << endl;

            }
        }

        if(save_phase_shift){
            //Index between 0 and 3600, this give us a histogram with 3600 bins
            unsigned index = (phase_shift[s])>0?uint(phase_shift[s]*1800.0/M_PI):uint(-phase_shift[s]*1800.0/M_PI);
            phase_shift_distribution(s,index)+=1;
        }

        if(subdivision_flag){
            for(uint i = 0 ; i < subdivisions.size(); i++){

                if( subdivisions[i].isInside(walker.pos_v)){
                    sub_DWI[i][s] += cos_phase_shift; // Real part
                    if(this->img_signal == true)
                        sub_DWIi[i][s]+= sin_phase_shift; // Img part

                    if(separate_signal){
                        if(walker.location == Walker::RelativeLocation::intra){
                            sub_DWI_intra[i][s]+=cos_phase_shift;
                        }
                        else if(walker.location == Walker::RelativeLocation::extra){
                            sub_DWI_extra[i][s]+=cos_phase_shift;
                        }
                    }

                    break;  //WARNING this break means that the subdivision are mutally exclusive
                }
            }
        }
        phase_shift[s] = 0;
    } //s


    //The for bellow is outside so it's not computed for each adquisition.
    if(subdivision_flag){
        for(uint i = 0 ; i < subdivisions.size(); i++){
            if( subdivisions[i].isInside(walker.pos_v)){

                subdivisions[i].density++;

                if(walker.intra_extra_consensus<0){
                        subdivisions[i].density_intra++;
                }
                else if(walker.intra_extra_consensus>0){
                    subdivisions[i].density_extra++;
                }
                break;  //WARNING this break means that the subdivision are mutally exclusive
            }
        }
    }
}

void GradientWaveform::setNumberOfSteps(unsigned T)
{
    this->T = T;
}
