//!  Basic class to store simulation parameters =============================================================/
/*!
*   \details   Basic class to store and handle all the possible simulation parameters.
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   1.5
*===========================================================================================================*/
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include "Eigen/Core"
#include "subdivision.h"
#include <utility>


/*! \class Parameter
 *  \brief Class used to hold and operate all the user and simulation parameters. This is the main class to comunicate
 *         between instances of the simulations and derived classes. So, in a way, it's an interface for the comunication
 *         between component classes in the simulation.
 */
class Parameters
{
public:

    unsigned num_walkers;                           /*!< N, number of walkers                                                       */
    unsigned num_steps;                             /*!< T, number of steps                                                         */
    double diffusivity;                             /*!< D, diffusivity constant                                                    */
    double sim_duration;                            /*!< simulation total time                                                      */
    bool write_traj;                                /*!< flag, write a traj file or not, binary format only                         */
    bool write_txt;                                 /*!< flag, writes DWI output signals in .txt if True                            */
    bool write_bin;                                 /*!< flag, writes the output signal in binary format (True by default)          */
    bool scale_from_stu;                            /*!< flag, true if the scheme file is in standar units m,s                      */
    bool save_phase_shift;                          /*!< flag, saves the phase shift distribution for all particles                 */
    long seed;                                      /*!< Initial seed for the random generator                                      */
    bool verbatim;                                  /*!< False to omit displaying state and warnings                                */
    std::string traj_file;                          /*!< Trajectory file path                                                       */
    std::string output_base_name;                   /*!< output files base name (path + sufix)                                      */
    std::string ini_walkers_file;                   /*!< initial walker position file (if any)                                      */
    unsigned ini_walkers_file_count;                /*!< number of walker positions initialize in the configuration file            */
    std::string ini_walker_flag;                    /*!< where to initialize the walkers                                            */
    std::string scheme_file;                        /*!< signal adquisition scheme file (if any)                                    */
    Eigen::Vector3d min_limits;                     /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits;                     /*!< voxel max limits (if any)                                                  */

    std::vector<std::string> cylinders_files;       /*!< file paths with a list of cilinders obstacles                              */
    std::vector<std::string> PLY_files;             /*!< file paths with PLY obstacle files                                         */
    std::vector<std::string> spheres_files;         /*!< file paths with spheres obstacle files                                     */
    std::vector<double> PLY_scales;                 /*!< Auxiliary vector to save PLY file scales                                   */
    std::vector<double> PLY_percolation;            /*!< Auxiliary vector to save PLY percolation                                   */

    std::vector<float> ini_delta_pos;               /*!< Delta position for the walkers                                             */

    unsigned num_proc;                              /*!< Number of precessors/process to launch in parallel                         */

    std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>>voxels_list; /*!< voxel min and max positions list   (deprecated)        */

    std::vector<Eigen::Vector3f> prop_dirs;         /*!< Saves the directions used to compute the propagator                        */

    std::vector<unsigned> record_pos_times;         /*!< time indexes, used to save the position of all walkers at certain time     */
    std::vector<unsigned> record_phase_times;       /*!< time indexes, used to save the phase shif of all walkers at certain time   */
    std::vector<unsigned> record_prop_times;        /*!< time indexes, used to save the mean propagator of the walkers at c. times  */

    bool   hex_cyl_packing;                         /*!< flag, true if an haxagonal packing should be used                          */
    bool   hex_sphere_packing;                      /*!< flag, true if an haxagonal packing OF SPHERES should be used               */
    double hex_packing_radius;                      /*!< float, constant radius for the cylinders                                   */
    double hex_packing_separation;                  /*!< float, separation distance betwen cylinders (separation > 2*radius)        */
    double hex_packing_icvf;                        /*!< float, ICVF computed or passed as parameter                                */

    bool        gamma_cyl_packing;                  /*!< flag, true if a gamma distribution of cylinders will be initialized        */
    bool        gamma_sph_packing;                  /*!< flag, true if a gamma distribution of  SPHERES will be initialized         */
    bool        gamma_output_conf;
    double      gamma_packing_alpha;
    double      gamma_packing_beta;
    double      gamma_icvf;
    double      gamma_output_configuration;
    unsigned    gamma_num_obstacles;
    float       min_obstacle_radii;                 /*!< Minimum radii (in um) to be sampled                                        */

    bool subdivision_flag           = false;        /*!< flag to check if we have several voxel subdivision to compute the signal   */
    unsigned number_subdivisions    = 0;            /*!< saves the number of subdivisions for an initialzied voxel (needed)         */
    std::string subdivisions_file   = "";           /*!< file with the list of subdivisions coordinates to compute the signal       */
    std::vector<Subdivision> subdivisions;          /*!< saves actual positions of the subdivision to compute the signal            */

    double obstacle_permeability     = 0;           /*!< Obstacles permeability                                                     */
    double collision_sphere_distance = 0;           /*!<  Custiom size for the collision sphere                                     */
    double max_simulation_time       = 0;           /*!< Maximum simulation time for the DYNAMIC SIMULATION                         */

    bool log_phase_shift = false;                   /*!< flag, true to save the final phase shift distribution                      */
    bool log_opp         = false;                   /*!< flag, true to save one per process output                                  */
    bool discard_stucks  = true;                    /*!< flag, true to discard posible stuck particles (max bouncing reached)       */
    bool discard_illegals = true;                   /*!< flag, true to discard possible illegal  crossings, Trump by default.       */

    bool log_propagator = false;                    /*!< flag, true saves the propagator for a given set of directions and times    */

    Eigen::Vector3d min_sampling_area;              /*!< Min defining point to delimiter the uniform sampling of walkers            */
    Eigen::Vector3d max_sampling_area;              /*!< Max defining point to delimiter the uniform sampling of walkers            */
    bool custom_sampling_area;                      /*!< True if a custom sampling area is defined (voxel for default)              */
    bool computeVolume;                             /*!< Forces the volumen computation (slower) even without custom sampling       */
    bool separate_signals;                          /*!< Separate the signals into intra and extra (compute_volume on)              */
    bool img_signal;                                /*!< True to save the img part of the dwi signal (false by default)             */


    /*! \fn Parameters
     *  \brief Default constructor. Sets all the parameters to default and NULL values.
     */
    Parameters();

    /*! \fn  readSchemeFile
     *  \param conf_file
     *  \brief Reads all the parameters from a scheme file in the correct format the function scales them if necessary.
     *         The parameters are passed by listing, first, the parameter name, followed by the value.
     *         The supported parameters are:
     *         number of walkers (N), number of steps (T), duration (duration), PGSE scheme file (scheme_file),
     *         min voxles limits (min limits), max voxel limits (max_limits), diffusivity (diffusivity),
     *         index name for the trajectory and output values (out_traj_file_index), initial walker position file (ini_walkers_file),
     *         write a txt traj flag and header (write_text), write binary traj file and header, write_bin, flag to scale the values
     *         from estandar unit (scale_from_stu), random seed (seed).
     */
    void readSchemeFile(std::string conf_file);

    //Set Methods:

    /*! \fn  setNumberWalkers
     *  \param N number of walkers
     *  \brief Set the number of walkers in the simulation.
    */
    void setNumWalkers(unsigned N);

    /*! \fn  setNumSteps();
     *  \param T number of steps
     *  \brief set the number of steps in the simulation.
    */
    void setNumSteps(unsigned T);

    /*! \fn  setDiffusivity
     *  \param Diff diffusivity value.
     *  \brief set the simulation diffusivity.
    */
    void setDiffusivity(double Diff);

    /*! \fn setSimDuration
     *  \param duration simulation duration.
     *  \brief sets the simulation duration.
    */
    void setSimDuration(double duration);

    /*! \fn setWriteFlag
     *  \param write_bin, boolean flag (0,1).
     *  \brief sets the flag to save a binary traj file.
    */
    void setWriteTrajFlag(bool write_bin);

    /*! \fn setWriteTxtFlag
     *  \param write_bin, boolean flag (0,1).
     *  \brief sets the flag to save a text traj file.
    */
    void setWriteTextFlag(bool write_txt_);

    /*! \fn setMinLimits
     *  \param min_limits_ vector with the minimum voxel limits (bottom left corner).
     *  \brief set the bottom left corner of the voxel to be simulated.
    */
    void setMinLimits(Eigen::Vector3d min_limits_);

    /*! \fn setMaxLimits
     *  \param max_limits_ vector with the maximum voxel limits (bottom right corner).
     *  \brief set the bottom left corner of the voxel to be simulated.
    */
    void setMaxLimits(Eigen::Vector3d max_limits_);

    /*! \fn setTrajFileName
     *  \param traj_file_ prefix of the traj file.
     *  \brief Set the prefix of the name for the traj file (txt and .traj)
    */
    void setTrajFileName(std::string traj_file_);

    /*! \fn setOutputBaseFileName
     *  \param output_base_name prefix for the outputs
     *  \brief Set the prefix of the name for all the outputs in the simulation.
    */
    void setOutputBaseFileName(std::string output_base_name_);

    /*! \fn iniWalkerFileName
     *  \param ini_walker_file_ assigns a file containing the initial position for the walkers.
     *  \brief Set the name of the file where the walkers initial positions should be read.
    */
    void iniWalkersFileName(std::string ini_walkers_file_);

    /*! \fn setSchemeFileName
     *  \param scheme_file_ scheme (PGSE )file name.
     *  \brief Sets the scheme file name to be used for the data synthesis.
    */
    void setSchemeFileName(std::string scheme_file_);

    // Get Methods

    /*! \fn getNumWalkers()
     *  \return Number of walkers N
    */
    unsigned getNumWalkers();

    /*! \fn getNumSteps()
     *  \return Number of Steps
    */
    unsigned getNumSteps();

    /*! \fn getDiffusivity
     *  \return Diffusivity
    */
    double getDiffusivity();

    /*! \fn getWriteTrajFlag
     *  \return flag of the  binary traj file writer
    */
    bool getWriteTrajFlag();

    /*! \fn getWriteTextFlag
     *  \return flag of the text write traj
    */
    bool getWriteTextFlag();

    /*! \fn getMinLimits
     *  \return voxel min limits (left bottom corner)
    */
    Eigen::Vector3d getMinLimits();

    /*! \fn getMaxLimits
     *  \return voxel max limits (right top corner)
    */
    Eigen::Vector3d getMaxLimits();

    /*! \fn getTrajFileName
     *  \return trajectory prefix
    */
    std::string getTrajFileName();

    /*! \fn getOutputBaseFileName
     *  \return Output prefix
    */
    std::string getOutputBaseFileName();

    /*! \fn iniWalkersFileName
     *  \return initial position walkers file name
    */
    std::string getIniWalkersFileName();

    /*! \fn getSchemeFileName
     *  \return name of the scheme file name used (PGSE)
    */
    std::string getSchemeFileName();

    static int str_dist(std::string s, std::string t);

    /*! \fn
     *  \brief adds the number of given subdivisions for the voxel
    */
    void addSubdivisions();

private:

    /*! \fn readObstacles
     *  \param file input iostreams
     *  \brief reads the full list of obstacles on the configuration file.
    */
    void readObstacles(std::ifstream &in);

    /*! \fn readVoxels
     *  \param file input iostreams
     *  \brief reads the full list of voxel on the configuration file.
    */
    void readVoxels(std::ifstream& in);

    /*! \fn r readInfoGatheringParam
     *  \param file input iostreams
     *  \brief reads the the list of times and positions to be recorded and written
     *  in the simulation
    */
    void readInfoGatheringParams(std::ifstream& in);

    /*! \fn readHexagonalParams
     *  \param file input iostreams
     *  \brief reads the parameters needed to define an hexagonal packing of cylinders
    */
    void readHexagonalParams(std::ifstream& in);

    /*! \fn readGammaParams
     *  \param file input iostreams
     *  \brief reads the parameters needed to define an gamma distributed packing of cylinders
    */
    void readGammaParams(std::ifstream& in);

    /*! \fn readHexagonalParams
     *  \param file input iostreams
     *  \brief reads the subdivisions for computing the DW signal
    */
    void readSubdivisionFile();

    /*! \fn
     *  \brief read the directions used to compute the propagator.
    */
    void readPropagatorDirections(std::string dir_path);

    /*! \fn
     *  \brief read a file with one ply scale value and then a list of ply files .
    */
    void readPLYFileList(std::string path);

    /*! \fn
     *  \brief read a file with the path, scale and percolation of each ply
     *
    */
    void readPLYFileListScalePercolation(std::string path);


};

#endif // PARAMETERS_H
