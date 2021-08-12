#include "pgsesequence.h"
#include "constants.h"
#include "Eigen/Dense"
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <assert.h>

using namespace std;

PGSESequence::PGSESequence()
{
    num_rep = 0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    T = 0;
    separate_signal=false;
}


PGSESequence::PGSESequence(Scheme scheme_)
{
    num_rep=0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    readSchemeParameters(scheme_);
    phase_shift_distribution.resize(num_rep,3600);
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
}

PGSESequence::PGSESequence(Scheme scheme_, const char *traj_file_name)
{
    num_rep = 0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    readSchemeParameters(scheme_);
    trajectory.setTrajFile(traj_file_name);
    T = int(trajectory.T);
    dyn_duration = trajectory.dyn_duration;
    phase_shift_distribution.resize(scheme_.num_rep,3600);
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
}

PGSESequence::PGSESequence(const char *scheme_file_name)
{
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    scheme_file = scheme_file_name;
    readSchemeFile();
    dyn_duration = scheme[0][6];
    phase_shift_distribution.resize(num_rep,3600);
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
    T = 10; //dummy number
}

PGSESequence::PGSESequence(const char *scheme_file_name, const char *traj_file_name)
{
    num_rep = 0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    scheme_file = scheme_file_name;
    readSchemeFile();
    trajectory.setTrajFile(traj_file_name);
    T = int(trajectory.T);
    dyn_duration = trajectory.dyn_duration;
    phase_shift_distribution.resize(num_rep,3600);
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
}

PGSESequence::~PGSESequence()
{
}

void PGSESequence::getGradImpulse(int grad_index, double t, double tLast, Eigen::Vector3d& Gdt){

    for(int i = 0; i < 3; i++)
        Gdt[i] = 0;

    double g[3]  = {scheme[grad_index][0],scheme[grad_index][1],scheme[grad_index][2]};
    double G     =  scheme[grad_index][3];
    double Delta =  scheme[grad_index][4];
    double delta =  scheme[grad_index][5];
    double te    =  scheme[grad_index][6];

    //printf("%.25f - %.25f - %.25f - %.25f - %.25f - %.25f - %.25f - \n",g[0],g[1],g[2],G,Delta,delta,te );
    //cout << " " << g[0] << " " << g[1] << " " << g[2] << " " << G << " " << Delta << " " << delta << " " << te << endl;

    if (!(t >= 0.0 && t<=te)){
        return;
    }

    //    printf("%d - %.25f - %.25f \n",grad_index,t,tLast);

    double pad = (te - Delta - delta)/2.0;
    if ( (t < pad) || (t > te-pad)){
        return;
    }

    double firstBlockStart = pad;
    double firstBlockEnd = pad+delta;
    //between pad and first block
    double sgn = 1;
    if( t >=firstBlockStart && t < firstBlockEnd){
        if(tLast < firstBlockStart){
            double dt = t - firstBlockStart;
            for (int j=0; j < 3; j++){
                Gdt[j] = sgn*G*dt*g[j];
            }
            return;
        }
    }

    double secondBlockStart = pad+Delta;
    double secondBlockEnd = pad+Delta+delta;

    //between the 2 blocks
    sgn = 1;
    if( t >=firstBlockEnd && t < secondBlockStart){
        if(tLast < firstBlockEnd){
            double dt= firstBlockEnd-tLast;
            for (int j=0; j < 3; j++){
                Gdt[j] = sgn*G*dt*g[j];
            }
            return;
        }
        return;
    }

    //segundo bloque
    if (t >= secondBlockStart){
        sgn=-1;
    }

    //if after second block
    if (t >= secondBlockEnd){
        if (tLast<secondBlockEnd){
            // the block ended between this call and the last one
            // so need to calculate the partial contribution
            double dt = secondBlockEnd-tLast;
            for (int j=0; j < 3; j++){
                Gdt[j] = sgn*G*dt*g[j];
            }
            return;
        }
        return;
    }

    if((t >= secondBlockStart)&&( tLast < secondBlockStart)){
        // the block ended between this call and the last one
        // so need to calculate the partial contribution
        double dt= t-secondBlockStart;

        for (int j=0; j < 3; j++){
            Gdt[j] = sgn*G*dt*g[j];
        }
        return;
    }
    for (int j=0 ; j < 3; j++){
        Gdt[j] = sgn*G*(t-tLast)*g[j];
    }
}


void PGSESequence::readSchemeParameters(Scheme scheme_){

    scheme_file = scheme_.scheme_file;
    dyn_duration = scheme_.scheme[0][6];

    num_rep = scheme_.scheme.size();

    for(unsigned i = 0 ; i < num_rep; i++){
        DWI.push_back(0);
        DWIi.push_back(0);

        phase_shift.push_back(0);
        scheme.push_back(scheme_.scheme[i]);
    }
}

void PGSESequence::readSchemeFile()
{
    ifstream in(scheme_file.c_str());

    //TODO: Error handling
    if(!in.is_open()){
        cout << "[ERROR] Can't open the scheme file " << endl;
        in.close();
        return;
    }

    vector<double> scheme_line;
    double tmp;
    string header;
    in >> header;
    in >> header;
    num_rep = 0;
    while( in >> tmp){
        scheme_line.push_back(tmp);
        DWI.push_back(0);

        if(this->img_signal == true)
            DWIi.push_back(0);

        if(separate_signal){
            DWI_extra.push_back(0);
            DWI_intra.push_back(0);
        }

        num_rep++;
        for(int i = 0 ; i < 6; i++){
            in >> tmp;
            scheme_line.push_back(tmp);
        }
        scheme.push_back(scheme_line);
        scheme_line.clear();
    }

    in.close();
}

void PGSESequence::update_phase_shift(double dt, double dt_last, Walker walker)
{
    Eigen::Vector3d xt;
    Eigen::Vector3d Gdt;
    //Displacement
    xt[0] = walker.pos_r[0] - walker.ini_pos[0];
    xt[1] = walker.pos_r[1] - walker.ini_pos[1];
    xt[2] = walker.pos_r[2] - walker.ini_pos[2];

    double dos_pi = 2.0*M_PI;

    for(int s=0; s < num_rep ;s++){
        getGradImpulse(s,dt,dt_last,Gdt);
        double val = giro*(Gdt[0]*xt[0]+Gdt[1]*xt[1]+Gdt[2]*xt[2]);
        val = fmod(val,dos_pi);
        phase_shift[s] = fmod(phase_shift[s] + val,dos_pi);
    }
}

void PGSESequence::update_phase_shift(double time_step, Eigen::Matrix3Xd trajectory)
{
    Eigen::Vector3d xt;
    Eigen::Vector3d Gdt;
    double dt,dt_last;

    for (uint t=1; t < this->T ;t++){ //TODO: checar si deberia ser <= T
        //Displacement
        xt[0] = trajectory(0,t) - trajectory(0,0);
        xt[1] = trajectory(1,t) - trajectory(1,0);
        xt[2] = trajectory(2,t) - trajectory(2,0);

        double dos_pi = 2.0*M_PI;

        if(this->dynamic){
            dt_last = this->time_steps[t-1];
            dt      = this->time_steps[t];
        }
        else{
            dt_last = time_step*(t-1);
            dt      = time_step*(t);
        }

        for(int s=0; s < num_rep ;s++){

            getGradImpulse(s,dt,dt_last,Gdt);
            double val = giro*(Gdt[0]*xt[0]+Gdt[1]*xt[1]+Gdt[2]*xt[2]);
            val = fmod(val,dos_pi);
            phase_shift[s] = fmod(phase_shift[s] + val,dos_pi);
        }
    }
}

void PGSESequence::update_DWI_signal(Walker& walker)
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

                    if(this->img_signal)
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


void PGSESequence::setNumberOfSteps(unsigned T)
{
    this->T = T;
}

void PGSESequence::computeDynamicTimeSteps()
{
    double Delta =  scheme[0][4];
    double delta =  scheme[0][5];
    double TE    =  scheme[0][6];
    double pad   = (TE - Delta - delta)/2.0;

    unsigned steps_in = percent_steps_in*T;

    //we want them to be even
    if(steps_in%2)
        steps_in++;

    delta = delta + delta*delta/20;

    int steps_pad = (2.0*pad)/(2.0*pad + Delta - delta) * (T-steps_in);

    //we want them to be even
    if(steps_pad%2)
        steps_pad++;

    int steps_out = T - steps_in - steps_pad;

    if( steps_in <= 0 || steps_out <= 0 || T<=0 || steps_pad <=0 || percent_steps_in <= 0){
        cout << "[Error] Incoherent number of steps inside the gradient pulse!" << endl;
        assert(0);
    }

    time_steps.resize(T+1,1);

    double dt_pad = (2*pad) / double(steps_pad);

    double dt_out = (Delta-delta)/double(steps_out);

    double dt_in  = (2.0*delta)/double(steps_in);

    ulong count    = 0.0;
    double time    = 0.0;


    for(int i=0;i < steps_pad/2.0; i++){
        time_steps[count++] = time;
        time += dt_pad;
    }


    for(int i=0;i < steps_in/2.0; i++){
        time_steps[count++] = time;
        time += dt_in;
    }


    for(int i=0;i < steps_out; i++){
        time_steps[count++] = time;
        time += dt_out;
    }

    for(int i=0;i < steps_in/2.0; i++){
        time_steps[count++] = time;
        time += dt_in;
    }

    for(int i=0;i <= steps_pad/2.0; i++){
        time_steps[count++] = time;
        time += dt_pad;
    }


    if(count != T+1){
        cout << "WARNING! T was not fullilled correctly in the dynamic setting!" <<endl;
    }

    //    for(int i = 0 ; i < T+1; i++)
    //        cout << time_steps[i] << endl;

}


double PGSESequence::getbValue(unsigned i)
{
    double G     =  scheme[i][3];
    double Delta =  scheme[i][4];
    double delta =  scheme[i][5];

    return (G*delta*giro)*(G*delta*giro)*(Delta - delta/3);
}

double PGSESequence::getFreeDecay(unsigned i,double D){
    double b = getbValue(i);

    return exp(-b*D);
}


double PGSESequence::getNumericalbValue(unsigned i)
{
    return -i;
}

void PGSESequence::getDWISignal()
{

    trajectory.initTrajReaderFile();

    trajectory.readTrajectoryHeader();

    double N        = trajectory.N;
    double T        = trajectory.T;
    double duration = trajectory.dyn_duration;
    double rt       = duration/T;
    double dos_pi   = 2.0*M_PI;
    double dt,dt_last,xt[3];

    Eigen::Matrix3Xd steps_log; // complete trajectory of one walker

    Eigen::Vector3d Gdt;
    Eigen::VectorXd phase_shift;

    steps_log.resize(3,unsigned(T+1));
    phase_shift.resize(num_rep);

    for (int w = 0; w < N; w++)
    {
        trajectory.readCurrentWalkersTrajectory(steps_log);
        for (uint t = 1; t <= uint(T); t++)
        {
            dt      = rt*(t);
            dt_last = rt*(t-1.0);

            xt[0] = steps_log(0,t) - steps_log(0,0);
            xt[1] = steps_log(1,t) - steps_log(1,0);
            xt[2] = steps_log(2,t) - steps_log(2,0);

            for(int s=0; s < num_rep ;s++)
            {
                getGradImpulse(s,dt,dt_last,Gdt);
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

}// END getDWISignal

double PGSESequence::get_adt(int grad_index, double t, double tLast){

    double Delta =  scheme[grad_index][4];
    double delta =  scheme[grad_index][5];
    double te    =  scheme[grad_index][6];

    //printf("%.25f - %.25f - %.25f - %.25f - %.25f - %.25f - %.25f - \n",g[0],g[1],g[2],G,Delta,delta,te );
    //cout << " " << g[0] << " " << g[1] << " " << g[2] << " " << G << " " << Delta << " " << delta << " " << te << endl;

    if (!(t >= 0.0 && t<=te)){
        return -INFINITY_VALUE;
    }

    //    printf("%d - %.25f - %.25f \n",grad_index,t,tLast);

    double pad = (te - Delta - delta)/2.0;
    if ( (t < pad) || (t > te-pad)){
        return 0;
    }

    double firstBlockStart = pad;
    double firstBlockEnd = pad+delta;
    //between pad and first block
    double sgn = 1;
    if( t >=firstBlockStart && t < firstBlockEnd){
        if(tLast < firstBlockStart){
            double dt = t - firstBlockStart;

            return sgn*dt;
        }
    }

    double secondBlockStart = pad+Delta;
    double secondBlockEnd = pad+Delta+delta;

    //between the 2 blocks
    sgn = 1;
    if( t >=firstBlockEnd && t < secondBlockStart){
        if(tLast < firstBlockEnd){
            double dt= firstBlockEnd-tLast;
            return sgn*dt;
        }
        return 0;
    }

    //segundo bloque
    if (t >= secondBlockStart){
        sgn=-1;
    }

    //if after second block
    if (t >= secondBlockEnd){
        if (tLast<secondBlockEnd){
            // the black ended between this call and the last one
            // so need to calculate the partial contribution
            double dt = secondBlockEnd-tLast;

            return sgn*dt;
        }
        return 0 ;
    }

    if((t >= secondBlockStart)&&( tLast < secondBlockStart)){
        // the block ended between this call and the last one
        // so need to calculate the partial contribution
        double dt= t-secondBlockStart;

        return sgn*dt;;
    }

    return sgn*(t-tLast);
}
