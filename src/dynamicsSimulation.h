//!  Dynamic simulation main class =============================================================/
/*!
*   \details   Main implementation of the particles dynamics. Handles collisions and bouncing
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.42
*==============================================================================================*/

#ifndef DynamicsSimulation_H
#define DynamicsSimulation_H

#include "walker.h"
#include <string>
#include "Eigen/Core"
#include <iostream>
#include <random>
#include "trajectory.h"
#include "simulablesequence.h"
#include "parameters.h"
#include "plyobstacle.h"
#include "voxel.h"
#include "cylinder.h"
#include "sentinel.h"
#include "propagator.h"
#include "sphere.h"


/*! \class DynamicsSimulation
 *  \brief  Main class, implements the particles dynamics. Handles collisions and bouncing.
 */
class DynamicsSimulation
{

public:
    Parameters   params;                            /*!< Parameters handler instance                                                */
    Walker       walker;                            /*!< Single walker to diffuse                                                   */
    Trajectory   trajectory;                        /*!< Trajectory instance. Handles i/o operations                                */
    std::mt19937 mt;                                /*!< rnd, random generator instance                                             */
    double step_lenght;                             /*!< l, step length                                                             */
    double second_passed;                           /*!< Simulation total time in seconds                                           */
    double max_simulation_time;                     /*!< Maximum simulation time if not passed we carry all the particles           */
    double completed;                               /*!< Auxiliar variable to save the milestone of percentage of completed walkers */
    std::string  ini_pos_file;                      /*!< walkers intitial position file                                             */
    unsigned ini_pos_file_ini_index;                /*!< starting position in the ini walker position file (multicore support)      */
    int id;                                         /*!< Unique id for the dynamic simulation                                       */
    sentinels::Sentinel sentinela;                  /*!< Sentinel initialization to encoutner error in the simulation               */
    std::vector <PLYObstacle>* plyObstacles_list;   /*!< pointer to a vector with all the instances of PLYObstacles                 */
    std::vector <Cylinder>* cylinders_list;         /*!< pointer to a vector with all the isntances of "Cylider" obstacles          */
    std::vector <Sphere>* spheres_list;           /*!< pointer to a vector with all the isntances of "Spheres" obstacles          */
    std::vector<unsigned>  cylinders_deque;         /*!< deque with the indexes of the cylinders (used for optmization)             */
    std::vector<unsigned>  spheres_deque;           /*!< deque with the indexes of the spheres (used for optmization)               */
    std::vector<std::vector<unsigned>> ply_deque;   /*!< deque with the indexes of the triangles of all ply's (used for opt)        */
    std::vector <Voxel> voxels_list;                /*!< vector with all the voxels to be simulated (if any)                        */
    Propagator propagator;                          /*!< Propagator object to compute and save the particles MSD                    */
    double icvf;                                    /*!< Stores the ICVF (1 - Intra-Extra) if needed                                */
    unsigned intra_tries, total_tries;              /*!< Helper variables to compute the estimated ICVF                             */

    /******   Auxiliar variables   ********/
    Eigen::Vector3d step;

    double time_step, time_dt, last_time_dt;        /*!< simulation time steps auxiliar values                                      */

    std::ifstream iniPos;

    time_t start,now;                               /*!< Auxiliar Variable for time recording and estimation for time.              */

    bool print_expected_time;                       /*!< Auxiliar flag for time recording and stimation for time.                   */

    unsigned num_simulated_walkers;                 /*!< Saves the final number of simulated walkers (time limit)                   */

    unsigned aux_walker_index;
    /****** END Auxiliar variables ********/


    /*! \fn  DynamicsSimulation
     *  \brief Default constructor. Initialize everything with 0's and NULL states, object indexes are set to -1.
     */
    DynamicsSimulation();

    /*! \fn     Main constructor.
     *  \param  conf_file: configuration file listing the simulation parameters in a given format.
     *  \brief  Reads all the parameters listed in the params conf_file and stores them in the /t params object.
     */
    DynamicsSimulation(std::string conf_file);

    /*! \fn     Alternative constructor.
     *  \param  params_ Parameters instance.
     *  \brief  Copies all the parameters from /t params_ and stores them in the /t params object.
     */
    DynamicsSimulation(Parameters &params_);

    //! \fn Default destructor.
    /*! \brief Does nothing.
    */
    ~DynamicsSimulation();

    /*! \fn     startSimulation.
     *  \param  dataSynth optional paramter. If this parameter is not given, no signal is computed.
     *  \brief  Starts the dynamics simulation and, if a PGSE sequence is given, computes the DW signal.
     */
    void startSimulation(SimulableSequence* dataSynth = nullptr);

    /*! \fn     readConfigurationFile
     *  \param  conf_file_path paremeters file path.
     *  \brief  Reads all the parameters listed in the param conf_file and stores them in the /t params object.
     */
    void readConfigurationFile(std::string conf_file_path);

    //SET VALUES FUNCTIONS
    /*! \fn  setDuration
     *  \param duration simulation duration.
     *  \brief Sets the simulation duration in milliseconds, this should be synchronized w/r the Time Echo.
     */
    void setDuration(const double& duration);

    /*! \fn    setWalkerNum
     *  \param N number of walkers for the simulation.
     *  \brief set the number of walkers N, for the simulation.
     */
    void setWalkersNum(const unsigned& N);

    /*! \fn    serStepsNum
     *  \param T number of steps to perform for each walker
     *  \brief Set the number of steps to perform for each walker.
     */
    void setStepsNum(const unsigned &T);


    static std::string secondsToMinutes(double);

    /*! \fn    isInIntra
     *  \param position 3d position on space.
     *  \param error minimum distance to be considered "outside" de obstacle (barrier thickness)
     *  \brief return true if the position is inside any of the obstacles. Only obstacles
     *         with a defined "inside region" can be considered. Voxel periodicity is not
     *         considered
     */
    bool isInIntra(Eigen::Vector3d& position, int& cyl_id,  int& ply_id, int& sph_id, double distance_to_be_intra_ply=1e-6);

    /*!
     * \brief   Writes to disk the final propagator matrix.
     */
    void writePropagator(std::string path);

    bool isInsideCylinders(Eigen::Vector3d& position,int& cyl_id,double distance_to_be_inside=1e-6);

    bool isInsidePLY(Eigen::Vector3d& position,int& ply_id,double distance_to_be_inside=1e-6);

    bool isInsideSpheres(Eigen::Vector3d &position, int& sph_id,double distance_to_be_inside);


private:    
    /*! \fn     generateStep
     *  \param  step stores the computed step.
     *  \param  l step size. Can be used to change diffusivity in the medium.
     *  \brief  Computes a random generated orientation in the sphere with given norm.
     *  \todo   Enable the use of pre-computed steps.
     */
    inline void generateStep(Eigen::Vector3d& step , double l);

    /*! \fn     generateDirectedStep
     *  \param  new_step stores the computed step.
     *  \param  direction the new step will be oriented toward this direction
     *  \brief  Computes a random generated orientation and oriented it according to a given direction.
     */
    inline void generateDirectedStep(Eigen::Vector3d& new_step , Eigen::Vector3d &direction);

    /*! \fn     updateWalkerPosition
     *  \param  step to be performed.
     *  \brief  updates the walker position in a step iteration. The methods checks for collisions against all
     *          stored obstacles and voxels and updates the walker position(s)
     *  \return returns false if the was any problem.
     */
    bool updateWalkerPosition(Eigen::Vector3d&step);

    /*! \fn     checkObstacleCollision
     *  \param  amended_step, step to be "amended", this is corrected against bouncing and voxel limits
     *  \param  tmax maximum step size, this value is updated every time a bouncing is performed.
     *  \param  end_points final position where the walker lands.
     *  \param  Collision, Collision instance to store the walker collision (if any).
     *  \brief  Checks for collisions against any obstacle or voxels given a direction and a step size. Only the
     *          more the collision with the higher priority is saved /see #Collision#.
     *  \return returns true if there was any collision.
     */
    inline bool checkObstacleCollision(Eigen::Vector3d& amended_step, double &tmax, Eigen::Vector3d &end_point, Collision &colision);

    /*! \fn     updateWalkerPositionAndHandleBouncing
     *  \param  amended_step, step to be "amended", this is corrected against bouncing and voxel limits
     *  \param  tmax maximum step size, this value is updated every time a bouncing is performed.
     *  \param  collision, Collision instance to store the walker collision (if any).
     *  \brief  Function to follow a collision event. After a successful collision given by /a checkObstacleCollision
     *          this function will handle the collision and decide to ignore it (percolation), bouncing, or label it as
     *          a special case.
     *  \return returns true if the collision was a correct bouncing.
     */
    inline bool updateWalkerPositionAndHandleBouncing(Eigen::Vector3d& amended_step, double& tmax, Collision& colision);

    /*! \fn     handleCollisions
     *  \param  collision A given collision with the highest priority so far.
     *  \param  collision_temp A second collision to compare with.
     *  \param  max_collision_distance: maximum collision distance to be considered successful.
     *  \param  ind unique walker indx identifier. Not used in the version 0.2.
     *  \brief  Warping function to handle the priority between 2 collision in a single step. The method uses the inhered comparisson
     *          in the class \see :Obstacle:.
     */
    inline void handleCollisions(Collision &colision, Collision &colision_tmp,  double &max_collision_distance, unsigned indx);

    /*! \fn     mapWalkerIntoVoxer
     *  \brief  If a voxel is given, maps the walker position back into the voxel, assuming a crossing at
     *          the voxel limits.
     */
    inline void mapWalkerIntoVoxel(Eigen::Vector3d &amended_step, Collision &colision, double barrier_thickness);

    /*! \fn     getTimeDt
     * \param   last_time_dt saves the last time step;
     * \param   time_dt actual time step;
     * \param   new step size if the time was dynamic
     * \param   dataSynt the PGSE sequence, if the steps are dynamic.
     * \param   t number of step.
     * \param   time_step size in milliseconds between steps.
     * \brief   Computes the step time. If the time steps are not dynamic this is just a constant sum.
     */
    inline void getTimeDt(double &last_time_dt, double& time_dt, double& l, SimulableSequence* dataSynth, unsigned t, double time_step);

    /*! \fn     initSimulation
     * \brief   Initialize simulation variables and write (if needed) header files.
     */
    inline void initSimulation();

    /*! \fn     displayExpectedTime
     * \brief   displays the remained (expected) simulation time, based in the time passed so far
     *          and the number of remaining walkers.
     */
    inline bool expectedTimeAndMaxTimeCheck(unsigned w);

    /*! \fn     writeDWSignal
     * \param   dataSynth Simuleable sequence used to the data Sythesis. NULL assumed to skip.
     * \brief   computes and writes the resulting diffusion signal for all the shells.
     */
    inline void writeDWSignal(SimulableSequence *dataSynth);

    /*! \fn     iniWalkerPosition
     * \brief   initialize the first walker position depending if a file was passed, the voxel limits,
     *          ot any other flag (as it can be intra, extra, delta position (not implemented yet)).
     * \todo    Add the flags " onlyIntra", "onlyExtra" and "singlePos".
     */
    inline void iniWalkerPosition();

    /*!
     * \brief   fill the list of indexes in walkers such that the obstacle is close enough for collision.
     * \todo    Implement the function minDistance for PLY's obstacles
     */
    inline void updateWalkerObstacleIndexes(unsigned t_);

    /*!
     * \brief   Initialize the list of obstacles indexes for the collision optimization.
     * \todo    Test the initialization for all the types of obstacles.
     */
    inline void initWalkerObstacleIndexes();

    /*!
     * \brief   Updates the list of indexes inside the inner and outher collision spheres.
     * \param   t number of steps in the simulation. Used to estimate the diffusion coeff.
     */
    inline void updateCollitionSphere(unsigned t);

    /*!
     * \brief   finds an intra celullar 3d position inside the voxel (needs a voxel initialized).
     * \param   intra_pos vector to save the 3d position.
     */
    inline void getAnIntraCellularPosition(Eigen::Vector3d& intra_pos, int &cyl_ind, int &ply_ind, int &sph_ind);

    /*!
     * \brief   finds an extra cellular 3d position inside the voxel (needs a voxel initialized).
     * \param   extra_pos vector to save the 3d position.
     */
    inline void getAnExtraCellularPosition(Eigen::Vector3d& extra_pos);

    /*!
     * \brief   Auxiliary function to checks if a 3d position is still inside the voxel
     * \param   pos to check inside the voxel.
     */
    inline bool checkIfPosInsideVoxel(Eigen::Vector3d& pos);

    /*!
     * \brief   Auxiliary function to initialize the permeability as well as the list of obstacle's
     *          indexes for the collision optimization procedures.
     */
    inline void initObstacleInformation();

    /*!
     * \brief   Function to internally update the log of the propagator.
     */
    inline void updatePropagator(Eigen::Matrix3Xd& log_pos_r);


    /*!
     * \brief   Function to internally normaliza the propagator using the final number of
     *          simualted signals.
     */
    inline void normalizePropagator(float num_samples);

    inline void computeICVF();

    inline bool finalPositionCheck();



};


#endif //DynamicsSimulation_H
