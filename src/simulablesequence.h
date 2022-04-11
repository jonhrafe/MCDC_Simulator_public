//!  MR Sequence Primary Class =============================================================/
/*!
  Elemental base clase.
  Abstract class to synthesise the MRI signal

  \author Jonathan Rafael
*==========================================================================================*/

#ifndef SIMULABLESEQUENCE_H
#define SIMULABLESEQUENCE_H
#include <string>
#include <vector>
#include "walker.h"
#include "subdivision.h"

class SimulableSequence{
public:

    std::string scheme_file;                    /*!< Scheme file path                                                           */

    std::vector<double> DWI;                    /*!< Real part of the  DWI signal                                               */

    std::vector<double> DWIi;                   /*!< imaginary partof the  DWI signal                                           */

    std::vector<double> phase_shift;            /*!< auxiliar phase shift for signal computations.                              */

    int num_rep;                                /*!< number of repetitions                         .                            */

    bool save_phase_shift;                      /*!< flag, if true, saves the pahse shift distribution.                         */

    bool dynamic;                               /*!< Flag to indicate if the time steps are non-uniform                         */

    double percent_steps_in;                    /*!< percentage of steps that should be inside the gradient times               */

    std::vector<double> time_steps;             /*!< Auxiliary array to save the time steps                                     */

    Eigen::ArrayXXf phase_shift_distribution;   /*!< Matrix to save the phase shif distribution                                 */

    std::vector<std::vector<double>> sub_DWI;   /*!< Real part of the DWI signal for each subDivision                           */

    std::vector<std::vector<double>> sub_DWIi;  /*!< Imaginary part of the  DWI signal for each subdivision                     */

    bool subdivision_flag           = false;    /*!< flag to check if we have several voxel subdivision to compute the signal   */

    std::vector<Subdivision> subdivisions;      /*!< saves the actual positions of the subdivision to compute the signal        */

    SimulableSequence(){}

    virtual ~SimulableSequence(){}

    /**
     * @param i: Walker index
     * @param t: current time step  (in milisenconds)
     * @param tLast: last time step (in milisenconds)
     * @param Gdt: vector to compute de G*dt impulse
     */
    virtual void getGradImpulse(int i, double t, double tLast, Eigen::Vector3d& Gdt) = 0;

    /**
     * @param i: index of the gradient in the scheme_file (0,N-1)
     * @return b-value
     */
    virtual double getbValue(unsigned i) {return i;} /*WARNING:Needs to be overloaded*/

    /**
     * @brief Expected free Decay
     */
    virtual double getFreeDecay(unsigned i,double D){return i*D;}  /*WARNING: Needs to be overloaded*/

    /**
     * @param i: updated walker
     */
    virtual void update_phase_shift(double dt,double dt_last,Walker walker) = 0;

    /**
     * @param i: updated the phase shift over a whole trajectory
     */
    virtual void update_phase_shift(double time_step,Eigen::Matrix3Xd trajectory) = 0;

    /**
     * @brief Updates the DWI signal using the cumulated phase shift
     */
    virtual void update_DWI_signal(Walker& walker) = 0;

    /**
     * @brief Set the number of time steps if they are known
     */
    virtual void setNumberOfSteps(unsigned T) = 0;

    /**
     * @brief Compute the time for all the steps when they are not constant.
     */
    virtual void computeDynamicTimeSteps() {}

    /**
     * @brief Initialize the DWI signals for each subdivision.
     */
    virtual void initializeSubdivisionSignals();

    virtual void writeResultingData(std::string output_base_name);

    virtual void writePhaseShiftDistribution(std::string output_base_name);

    virtual void cleanPhaseShift();

    virtual void cleanDWISignal();
};


#endif // SIMULABLESEQUENCE_H
