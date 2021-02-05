//!  Gradient Wavefroms =============================================================/
/*!
*   \details   Main implementation of the gradient waveforms protocol
*   \author    Jonathan Rafael
*   \date      November 2017
*   \version   1.42
*==============================================================================================*/

#ifndef GRADIENTWAVEFORM_H
#define GRADIENTWAVEFORM_H

#include "pgsesequence.h"
#include "Eigen/Dense"
#include "constants.h"

class GradientWaveform: public SimulableSequence
{
public:

    double TE;                  /*!< Time Echo.                                           */

    uint T;                      /*!< num bins (time steps)                                */

    double dyn_duration;        /*!< simulation duration (miliseconds)                    */

    int wave_bins;              /*!< Wave discretization                                  */

    double wave_duration;       /*!< Wave duration (should be less qeual than dyn_dur.)   */

    double dt;                 /*!< individual time steps (miliseconds) of the wave       */

    bool scale_from_stu;       /*!< True if the input is in standar units                 */

    std::vector< std::vector<float> >  waveform; /*!< Defined waveforms                   */

    Trajectory trajectory;     /*!< If the signal is computed from a .trajfile            */

    /**
     * @brief Default constructor, set default NULL values. Not to be used.
     */
    GradientWaveform();

    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file
     */
    GradientWaveform(Scheme& scheme);

    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file and a traj file name.
     *        if this argument is passed a traj file is should be written.
     */
    GradientWaveform(Scheme& scheme_,const char* traj_file_name);


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
      * @brief reads the waveform
      */
    void readSchemeFile();

    /**
      * @brief For using with the adt array
      */
    void getInterpolatedGradImpulse(uint rep_index, double dt_sim, double t_sim_last, Eigen::Vector3d &Gdt);

    /*
     * @brief Updates the phase shift using the full stored trajectory
     */
     void update_phase_shift(double time_step, Eigen::Matrix3Xd trajectory);


     void update_phase_shift(double dt,double dt_last,Walker walker);
    /**
      * @brief Updates the DWI signal using the cumulated phase shift
      */
     void update_DWI_signal(Walker &walker);

     void setNumberOfSteps(unsigned T);

     void getGradImpulse(int i, double t, double tLast, Eigen::Vector3d& Gdt);


private:
    /**
      * @brief reads the scheme file parameters
      */
    void readSchemeParameters(Scheme &scheme_);


};

#endif // GRADIENTWAVEFORM_H
