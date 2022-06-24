//!  PGSE Sequence Class  =============================================================/
/*!
  Derived Class.
  Implementation of the PGSE protocol

  \date   May 2016
  \author Jonathan Rafael
  \version 0.1.0
*=====================================================================================*/

#ifndef PGSESEQUENCECUDA_H
#define PGSESEQUENCECUDA_H

#include <vector>
#include <iostream>
#include "trajectory.h"
#include <Eigen/Core>
#include "simulablesequencecuda.h"
#include "scheme.h"


#include <cuda.h>
#include "cuda_runtime.h"

typedef unsigned long ulong;


__global__ void update_phase_shift_cuda(int N, double *traj, double *grad_sequence, int num_rep, int T, double *phase_shift_cuda);

__global__ void update_DWI_signal_cuda(int num_rep, double *phase_shift_cuda, int walker_batch_size, double *DWI_cuda );

__global__ void cudaCreateGradSequence(double time_step, int T, int num_rep, double *scheme_vec, double *grad_sequence);

__global__ void cudaInitGrad(int T, int num_rep, double* grad_sequence);

__global__ void cudaInitDWI(int walker_batch_size, int num_rep, double* DWI_cuda);

__global__ void cudaInitPhaseShift(int walker_batch_size, int num_rep, double *phase_shift_cuda);


/*! \class  PGSESequenceCuda
 *  \brief  Implementation of the PGSE protocol on CUDA
 */



class PGSESequenceCuda: public SimulableSequenceCuda{
public:

    double TE;           /*!< Time Echo.                                        */

    int T;               /*!< num bins (time steps)                             */

    double dyn_duration; /*!< simulation duration (miliseconds)                 */

    std::vector< std::vector<double> > scheme;  /*!< Scheme file values         */

    Trajectory trajectory;    /*!< If the signal is computed from a .trajfile   */





    //constructors

    /**
     * @brief Default constructor, set default NULL values. Not to be used.
     */
    PGSESequenceCuda();
    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file.
     */
    PGSESequenceCuda(Scheme scheme_);

    /**
     * @brief Destuctor. Does nothing.
     */
    ~PGSESequenceCuda();

    /**
     * @brief For using w/o the adt array
     */
    void getGradImpulse(int i, double t, double tLast, Eigen::Vector3d &Gdt);


    /**
     * @brief For using with the adt array
     */
    void getGradImpuse(int i,  double t, Eigen::Vector3d Gdt);

    /**
     * @brief Analytical defined b-value
     */
    double getbValue(unsigned);

    /**
     * @brief Expected free Decay
     */
    double getFreeDecay(unsigned i,double D);


    //TODO: to implement
    /**
     * @brief \warning not implemented yet.
     */
    double getNumericalbValue(unsigned);

    /**
     * @brief Computes de DW signal from a trajfile
     */
    void getDWISignal();

    /**
     * @brief reads the scheme files
     */
    void readSchemeFile();

    /**
     * @brief Updates the phase shift using the full stored trajectory
     */
    virtual void update_phase_shift_DWI_signal(double *traj_mat, int walker_batch_size);

    /**
     * @brief computes de signal value and sign in a certain time step.
     */
    double get_adt(int grad_index, double t, double tLast);

    /**
     * @brief prints the array adt in the format ().
     */
    double print_adt_and_dt(int grad_index, double t, double tLast);

    virtual void setNumberOfSteps(unsigned T);

    virtual void computeDynamicTimeSteps();

    virtual void freeCudaVariables();

    virtual void initCudaVariables(int walker_batch_size);

    virtual void resetCudaVariables(int walker_batch_size);


private:
    virtual void readSchemeParameters(Scheme scheme_);

    void createGradSequence(double *grad_sequence);

    






};
#endif // PGSESEQUENCECUDA_H
