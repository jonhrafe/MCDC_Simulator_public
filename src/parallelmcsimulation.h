/*! ==============================================================================/
*   \details   Class to handle multiprocessor paralellisation
*   \author    Jonathan Rafael
*   \date      November 2016
*
*================================================================================*/
#ifndef PARALLELMCSIMULATION_H
#define PARALLELMCSIMULATION_H

#include "mcsimulation.h"
#include <thread>


/*! \class ParallelMCSimulation
 * \brief  Class to handle multiprocessor paralellisation. This class basicly controls and syncronize several
 *         initializations of MonteCarlo simulations and add up the results. It's a way of soft paralelization.
 */
class ParallelMCSimulation
{
public:

    Parameters params;                           /*!< Parameters instance \see :Parameters:                              */

    double mean_second_passed;                   /*!< Simualation total time in seconds                                  */
    unsigned total_sim_particles;                /*!< Total number of simulated particles                                */
    unsigned stuck_count;                        /*!< Counts the number of particles stuck in the simulations            */
    unsigned illegal_count;                      /*!< Counts the number of particles that attempt to cross               */
    double icvf;                                 /*!< Stores the ICVF based on the particles sampling                    */
    double aprox_volumen;                        /*!< Stores the volumen based on ICVF and the voxel size                */
    std::vector<MCSimulation*> simulations;      /*!< vector of pointers to MCSimulation instances                       */
    std::vector<std::thread> sim_threads;        /*!< Number of threads (instances and processors) to be used            */
    std::vector <PLYObstacle> plyObstacles_list; /*!< vector with all the instances of PLYObstacles                      */
    std::vector <Cylinder> cylinders_list;       /*!< vector with all the instances of cylinders                         */
    std::vector <Sphere> spheres_list;          /*!< vector with all the instances of cylinders                         */

    std::vector<Eigen::Vector3f> total_ini_walker_pos; /*!< Number of threads (instances and processors) to be used      */

    /*!
     *  \param config_file .conf file name with the full set of experiments parameters.
     *  \brief Main constructor.
     */
    ParallelMCSimulation(std::string config_file);

    /*!
     *  \param parameters of the simulation. Read form a conf file or given by the user.
     *  \brief Constructor.
     */
    ParallelMCSimulation(Parameters &params);


    ParallelMCSimulation(){}

    ~ParallelMCSimulation();

    /*! \fn  startSimulation
     *  \brief Warp function. Calls the MCSimulation's native function for all the instances. \see :MCSimulation:.
     */
    void startSimulation();

private:
    /*!
     *  \brief return the number of processors in your machine
     */
    void getNumberOfProcessors();

    /*!
     *  \brief Initialize every individual MCSimulation
     */
    void initializeUnitSimulations();

    /*!
     *  \brief Joints the ouputs in a single file.
     */
    void jointResults();

    /*!
     *  \brief Initialize anythin that needs to be sync between simulations.
    */
    void specialInitializations();

    /*!
     *  \brief Initialize all the Cylinders from a file.
    */
    void addObstaclesFromFiles();

    /*!
     *  \brief Used to initialize the hexagonal packing
    */
    void addObstacleConfigurations();

};

#endif // PARALLELMCSIMULATION_H
