//!  PGSE Sequence Class  =============================================================/
/*!
  Derived Class.
  Implementation of the PGSE protocol

  \date   May 2016
  \author Jonathan Rafael
  \version 0.1.0
*=====================================================================================*/

#ifndef PGSESEQUENCE_H
#define PGSESEQUENCE_H

#include <vector>
#include <iostream>
#include "trajectory.h"
#include "Eigen/Core"
#include "simulablesequence.h"
#include "scheme.h"
#include "constants.h"   // uint / ulong


/*! \class  ParallelMCSimulation
 *  \brief  Implementation of the PGSE protocol
 */
class PGSESequence: public SimulableSequence{
public:

    double TE;           /*!< Time Echo.                                        */

    int T;               /*!< num bins (time steps)                             */

    double dyn_duration; /*!< simulation duration (miliseconds)                 */

    std::vector< std::vector<double> > scheme;  /*!< Scheme file values         */

    /*!< perf B4: list of timesteps where the gradient impulse is non-zero for some
     *   direction. The PGSE gradient is on only during the two short delta-pulses,
     *   so for most timesteps Gdt==0 and the phase update is a no-op (phase is
     *   already range-reduced). These steps can be skipped -> bit-exact. The list is
     *   walker-independent: built ONCE and shared read-only across threads. When
     *   null, every timestep is processed (original behaviour).                    */
    const std::vector<unsigned>* grad_active_t = nullptr;

    Trajectory trajectory;    /*!< If the signal is computed from a .trajfile   */

    //constructors

    /**
     * @brief Default constructor, set default NULL values. Not to be used.
     */
    PGSESequence();
    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file.
     */
    PGSESequence(Scheme scheme_);
    /**
     * @brief Main constructor. Takes a pre-loaded Scheme file and a traj file name.
     *        if this argument is passed a traj file is should be written.
     */
    PGSESequence(Scheme scheme_,const char* traj_file_name);
    /**
     * @brief Constructor. Takes a the scheme file name to be loaded.
     */
    PGSESequence(const char* scheme_file_name);
    /**
     * @brief Constructor. Takes a scheme file name to be loaded and atraj file name.
     *        if this argument is passed a traj file is should be written.
     */
    PGSESequence(const char* scheme_file_name,const char* traj_file_name);

    /**
     * @brief Destuctor. Does nothing.
     */
    ~PGSESequence();

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

    /**
     * @brief Get echo time..
     */
    double getTE(unsigned);


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
     * @param i: updated walker
     */
    virtual void update_phase_shift(double dt,double dt_last,Walker walker);

    /**
     * @brief Updates the phase shift using the full stored trajectory
     */
    virtual void update_phase_shift(double time_step, Eigen::Matrix3Xd trajectory);

    /**
     * @brief perf B4: precompute (walker-independent) the list of timesteps that
     *        carry a non-zero gradient impulse for at least one direction. Other
     *        timesteps contribute exactly zero to the phase and are skipped, which
     *        is bit-exact (the phase is already range-reduced).
     */
    void buildActiveTimesteps(double time_step, std::vector<unsigned>& out) const;

    /** @brief Point this sequence at a shared, read-only active-timestep list. */
    void setActiveTimesteps(const std::vector<unsigned>* t){ grad_active_t = t; }

    /**
     * @brief Updates the DWI signal using the cumulated phase shift
     */
    virtual void update_DWI_signal(Walker &walker,double dt);

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

private:
    virtual void readSchemeParameters(Scheme scheme_);



};
#endif // PGSESEQUENCE_H
