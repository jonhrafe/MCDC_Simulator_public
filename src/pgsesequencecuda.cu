#include "pgsesequencecuda.cuh"
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


#include <cuda.h>
#include "cuda_runtime.h"

using namespace std;

__global__ void update_phase_shift_cuda(int walker_batch_size, double *traj, double *grad_sequence, int num_rep, int T, double *phase_shift_cuda){

    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int stride = blockDim.x * gridDim.x;

    /*
     - Template parameter for walker_batch_size 
     - parallel unroll
    */

    for (int idx_walker=i; idx_walker < walker_batch_size; idx_walker+=stride){

        if (idx_walker<walker_batch_size){

            int idx_time_s = idx_walker*3*(T+1); 
            double *curr_traj = &(traj[idx_time_s]);
            double pos_init[3]={curr_traj[0], curr_traj[1], curr_traj[2]};


            int idx_phase_shift_s = idx_walker*num_rep;
            double *curr_phase_shift = &(phase_shift_cuda[idx_phase_shift_s]);

            double dos_pi = 2.0*M_PI;

            double xt[3];

            double val = 0.0;
            
            int grad_sequence_idx = 0;

            for (int tt=1; tt <= T; tt++){ 
                //Displacement
                xt[0] = curr_traj[tt*3] - pos_init[0];
                xt[1] = curr_traj[tt*3+1] - pos_init[1];
                xt[2] = curr_traj[tt*3+2] - pos_init[2];

                
                for(int ss=0; ss < num_rep ; ss++){
                    grad_sequence_idx = 3*(ss*T + (tt-1)); 
                    val = (giro*(grad_sequence[grad_sequence_idx]*xt[0]+grad_sequence[grad_sequence_idx+1]*xt[1]+grad_sequence[grad_sequence_idx+2]*xt[2]));                
                    val = fmod(val, dos_pi);                    
                    curr_phase_shift[ss] = fmod(curr_phase_shift[ss]+val,dos_pi);
                }
            }
        }
    }
    return;
}


__global__ void update_DWI_signal_cuda(int num_rep, double *phase_shift_cuda, int walker_batch_size, double* DWI_cuda){

    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int stride = blockDim.x * gridDim.x;

    for (int idx_phase_shift=i; idx_phase_shift<num_rep;idx_phase_shift+=stride){

        //cudaMemCheck --> Check no outer-range index

        if (idx_phase_shift<num_rep){

            for (int j=idx_phase_shift; j<(idx_phase_shift+(num_rep*walker_batch_size)); j+=num_rep){

                DWI_cuda[idx_phase_shift] += cos(phase_shift_cuda[j]);
                

            }    
        }
    }
    return; // What about last elements ? idx_phase_shift - (num_rep*walker_batch_size)
    
}
  

__global__ void cudaCreateGradSequence( double time_step, int T, int num_rep, double *scheme_vec, double *grad_sequence){

	double tcurr, tlast;
    
    int acq_index = blockIdx.x*blockDim.x + threadIdx.x;
    int stride = blockDim.x*gridDim.x;
    int scheme_index=0;

    for (int i=acq_index;i<num_rep;i+=stride){

        if (i<num_rep){

            scheme_index = i*7;

            double g[3]  = { scheme_vec[scheme_index],scheme_vec[scheme_index+1], scheme_vec[scheme_index+2]};
            double G     =  scheme_vec[scheme_index+3];
            double Delta =  scheme_vec[scheme_index+4];
            double delta =  scheme_vec[scheme_index+5];
            double te    =  scheme_vec[scheme_index+6];
            double pad = (te - Delta - delta)/2.0;

            double firstBlockStart = pad;
            double firstBlockEnd = pad+delta;

            double secondBlockStart = pad+Delta;
            double secondBlockEnd = pad+Delta+delta;

            double sgn = 0.0;
            double dt = 0.0;

            int grad_idx=3*i*T;


            for(int t = 1 ; t <= T; t++){

                tlast = time_step*(double)(t-1.0);
                tcurr = time_step*(double)(t);

                // impulse condition
                if (!(tcurr >= 0.0 && tcurr<=te)){
                    sgn = 0.0;
                    dt = 0.0;
                } 

                else if ( (tcurr < pad) || (tcurr > te-pad)){
                    sgn = 0.0;
                    dt = 0.0;
                }

                else if ((tcurr >=firstBlockStart) && (tcurr < firstBlockEnd)){                
                    //between pad and first block
                    sgn = 1.0;
                    if(tlast < firstBlockStart){
                        dt = tcurr - firstBlockStart;
                    }
                    else{
                        dt=time_step;
                    }
                }

                else if ((tcurr >=firstBlockEnd) && (tcurr < secondBlockStart)){            
                    //between the 2 blocks
                    if(tlast < firstBlockEnd){
                        sgn = 1.0;
                        dt= firstBlockEnd-tlast;
                    }
                    else{
                        sgn = 0.0;
                        dt  = time_step;
                    }
                }

                else if((tcurr >= secondBlockStart)&&(tcurr <= secondBlockEnd)){
                    sgn = -1.0;
                    if (tlast < secondBlockStart){
                    // the block ended between this call and the last one
                    // so need to calculate the partial contribution
                        dt  = tcurr-secondBlockStart;
                    }
                    else{
                        dt=time_step;
                    }
                }
                else if (tcurr >= secondBlockEnd){
                    if (tlast<secondBlockEnd){
                        // the block ended between this call and the last one
                        // so need to calculate the partial contribution
                        sgn = -1.0;
                        dt = secondBlockEnd-tlast;
                    }
                    else{
                        sgn = 0.0;
                        dt=time_step;
                    }
                }
                else{
                    sgn= 0.0;
                    dt=time_step;
                }
                if (sgn!=0.0){
                    for (int j=0 ; j < 3; j++){
                        grad_sequence[grad_idx+j] = sgn*G*dt*g[j];
                    }
                }
                grad_idx+=3;
                
            }// end t
        }//end grad_index 
    }  
    return;
}

__global__ void cudaInitGrad(int T, int num_rep, double* grad_sequence){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int stride = blockDim.x * gridDim.x;

    for (int idx_grad=i; idx_grad<(3*num_rep*T);idx_grad+=stride){
        if (idx_grad<(3*num_rep*T)){
            grad_sequence[idx_grad]=0.0;
        }
    }
}

__global__ void cudaInitDWI(int walker_batch_size, int num_rep, double* DWI_cuda){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int stride = blockDim.x * gridDim.x;

    for (int idx_phase_shift=i; idx_phase_shift<num_rep;idx_phase_shift+=stride){
        DWI_cuda[idx_phase_shift]=0.0;
    }
}

__global__ void cudaInitPhaseShift(int walker_batch_size, int num_rep, double *phase_shift_cuda){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x; 
    int stride = blockDim.x * gridDim.x;

    for (int idx_walker=i; idx_walker < (walker_batch_size*num_rep); idx_walker+=stride){
        phase_shift_cuda[idx_walker]=0.0;
    }    
}

PGSESequenceCuda::PGSESequenceCuda()
{
    num_rep = 0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    T = 0;
    separate_signal=false;
}


PGSESequenceCuda::PGSESequenceCuda(Scheme scheme_)
{
    num_rep=0;
    dynamic = false;
    save_phase_shift = true;
    percent_steps_in = -1;
    readSchemeParameters(scheme_);
    phase_shift_distribution.resize(num_rep,3600);
    phase_shift_distribution = Eigen::ArrayXXf::Zero(num_rep,3600);
}




PGSESequenceCuda::~PGSESequenceCuda()
{
}

void PGSESequenceCuda::getGradImpulse(int grad_index, double t, double tLast, Eigen::Vector3d& Gdt){

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


void PGSESequenceCuda::readSchemeParameters(Scheme scheme_){

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

void PGSESequenceCuda::readSchemeFile()
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


void PGSESequenceCuda::createGradSequence(double *grad_sequence){

	double time_step	= this->dyn_duration/this->T;
	double tcurr(0.0), tlast(0.0);
    
    int grad_idx=0;

    double g[3]  = {0.0, .0, .0};
    double G     =  0;
    double Delta =  0;
    double delta =  0;
    double te    =  0;
    double pad =  0;

    double firstBlockStart =  0;
    double firstBlockEnd =  0;

    double secondBlockStart =  0;
    double secondBlockEnd =  0;

    double sgn = 0.0;
    double dt = 0.0;

    for (int grad_index=0 ; grad_index < this->num_rep ; grad_index++){

        g[0]  = this->scheme[grad_index][0];
        g[1]  = this->scheme[grad_index][1];
        g[2]  = this->scheme[grad_index][2];
        G     =  this->scheme[grad_index][3];
        Delta =  this->scheme[grad_index][4];
        delta =  this->scheme[grad_index][5];
        te    =  this->scheme[grad_index][6];
        pad = (te - Delta - delta)/2.0;

        firstBlockStart = pad;
        firstBlockEnd = pad+delta;

        secondBlockStart = pad+Delta;
        secondBlockEnd = pad+Delta+delta;

        sgn = 0.0;
        dt = 0.0;

        for(int t = 1 ; t <= T; t++){

            tlast = time_step*(double)(t-1.0);
            tcurr = time_step*(double)(t);

            // impulse condition
            if (!((tcurr >= 0.0) && (tcurr<=te))){
               sgn = 0;
               dt = 0.0;
            } 

            else if ((tcurr < pad) || (tcurr > te-pad)){
                sgn = 0;
                dt = 0.0;
            }

            else if ((tcurr >=firstBlockStart) && (tcurr < firstBlockEnd)){                
                sgn = 1;
                if(tlast < firstBlockStart){
                    //between pad and first block
                    dt = tcurr - firstBlockStart;
                }
                else{
                    dt=time_step;
                }
            }

            else if ((tcurr >=firstBlockEnd) && (tcurr < secondBlockStart)){            
                if(tlast < firstBlockEnd){
                    sgn = 1;
                    dt= firstBlockEnd-tlast;
                }
                else{
                    //between the 2 blocks
                    sgn = 0;
                    dt  = time_step;
                }
            }

            else if((tcurr >= secondBlockStart)&&(tcurr <= secondBlockEnd)){
                sgn = -1;
                if (tlast < secondBlockStart){
                // the block ended between this call and the last one
                // so need to calculate the partial contribution
                    dt  = tcurr-secondBlockStart;
                }
                else{
                    dt=time_step;
                }
            }
            else if (tcurr >= secondBlockEnd){
                if (tlast<secondBlockEnd){
                    // the block ended between this call and the last one
                    // so need to calculate the partial contribution
                    sgn = -1;
                    dt = secondBlockEnd-tlast;
                }
                else{
                    sgn = 0.0;
                    dt=0.0;
                }
            }
            else{
                sgn= 0.0;
                dt=time_step;
            }
            
            for (int j=0 ; j < 3; j++){
                grad_sequence[grad_idx+j] = (sgn*G*dt*g[j]);
            }
            grad_idx+=3;
            
        }// end t
    }//end grad_index   
}


void PGSESequenceCuda::initCudaVariables(int walker_batch_size){

    printf("Cuda init...");

    int blockSize = 256;
    int numBlocks = 0;

    cudaError_t err_malloc;

    // Grad sequence 

    err_malloc = cudaMallocManaged(&this->grad_sequence,  (3*(this->T)*(this->num_rep))*sizeof(double));
    
    if (err_malloc!=cudaSuccess){
        printf("var %d %d ", this->T, this->num_rep);
        printf("Memory Error grad sequence");
        return;
    }
    


    /*
    numBlocks = ((3*this->num_rep*this->T) + blockSize - 1) / blockSize;

    cudaInitGrad<<<numBlocks, blockSize>>>(this->T, this->num_rep, this->grad_sequence);

    numBlocks = (this->num_rep + blockSize - 1) / blockSize;

    double time_step	= this->dyn_duration/this->T;
    int nb_elem = (this->num_rep*this->scheme[0].size());

    double scheme_vec[nb_elem];
        
    int curr_idx = 0;
    for (int curr_acq=0;curr_acq<this->num_rep;curr_acq++){
        for (int jj=0;jj<this->scheme[curr_acq].size();jj++){
            curr_idx = curr_acq*(this->scheme[curr_acq].size())+jj;
            scheme_vec[curr_idx] = (double)(this->scheme[curr_acq][jj]);
        }
    }

    double *scheme_vec_ptr = scheme_vec;

    cudaCreateGradSequence<<<numBlocks, blockSize>>>(time_step, this->T, this->num_rep, scheme_vec_ptr, this->grad_sequence);
    */

    this->createGradSequence(this->grad_sequence);

   
    // Phase shift
    err_malloc=cudaMallocManaged(&this->phase_shift_cuda,  (walker_batch_size*this->num_rep)*sizeof(double));

    if (err_malloc!=cudaSuccess){
        printf("Memory Error phase_shift_cuda");
        return;
    }
    numBlocks = (walker_batch_size + blockSize - 1) / blockSize;
    cudaInitPhaseShift<<<numBlocks, blockSize>>>(walker_batch_size, this->num_rep, this->phase_shift_cuda);

    //DWI signal
    err_malloc=cudaMallocManaged(&this->DWI_cuda,  (this->num_rep)*sizeof(double));
    if (err_malloc!=cudaSuccess){
        printf("Memory Error DWI_cuda");
        return;
    }

    numBlocks = (this->num_rep + blockSize - 1) / blockSize;
    cudaInitDWI<<<numBlocks, blockSize>>>(walker_batch_size, this->num_rep, this->DWI_cuda);
    
    // Wait synchronization
    cudaDeviceSynchronize();

    printf("Done\n");

    return;

}

void PGSESequenceCuda::freeCudaVariables(){
    // Free memory
    cudaFree(this->phase_shift_cuda);
    cudaFree(this->grad_sequence);
    cudaFree(this->DWI_cuda);
    return;
    
}

void PGSESequenceCuda::resetCudaVariables(int walker_batch_size){

    int blockSize = 1;
    int numBlocks = 1;

    cudaInitPhaseShift<<<numBlocks, blockSize>>>(walker_batch_size, this->num_rep, this->phase_shift_cuda);
    cudaDeviceSynchronize();

    return;

}

void PGSESequenceCuda::update_phase_shift_DWI_signal(double* traj_mat, int walker_batch_size){   
    
    int blockSize = 256;
    int numBlocks = (walker_batch_size + blockSize - 1) / blockSize;


    printf("Phase shift update <%d, %d>...\n ", blockSize, numBlocks);
    //GPU Phase shift update
    update_phase_shift_cuda<<<numBlocks, blockSize>>>(walker_batch_size, traj_mat,  this->grad_sequence, this->num_rep, this->T, this->phase_shift_cuda);
    
    // Wait synchronization
    cudaDeviceSynchronize();
    printf("Done\n");

    // GPU signal update
    numBlocks = (this->num_rep + blockSize - 1) / blockSize;

    printf("DWI update<%d, %d>...\n ", blockSize, numBlocks);
    update_DWI_signal_cuda<<<numBlocks, blockSize>>>(this->num_rep, this->phase_shift_cuda, walker_batch_size, this->DWI_cuda);

    // Wait synchronization
    cudaDeviceSynchronize();
    printf("Done\n");

    for (unsigned int ii=0;ii<(this->num_rep);ii++){
        this->DWI[ii] = this->DWI_cuda[ii];
    }    
    return;
    
}



void PGSESequenceCuda::setNumberOfSteps(unsigned T)
{
    this->T = T;
}

void PGSESequenceCuda::computeDynamicTimeSteps()
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


double PGSESequenceCuda::getbValue(unsigned i)
{
    double G     =  scheme[i][3];
    double Delta =  scheme[i][4];
    double delta =  scheme[i][5];

    return (G*delta*giro)*(G*delta*giro)*(Delta - delta/3);
}

double PGSESequenceCuda::getFreeDecay(unsigned i,double D){
    double b = getbValue(i);

    return exp(-b*D);
}


double PGSESequenceCuda::getNumericalbValue(unsigned i)
{
    return -i;
}

void PGSESequenceCuda::getDWISignal()
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

double PGSESequenceCuda::get_adt(int grad_index, double t, double tLast){

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
