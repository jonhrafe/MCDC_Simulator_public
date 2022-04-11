//!  Basic class to store simulation parameters =============================================================/
/*!
*   \details   Basic class to store and handle all the possible simulation parameters.
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   0.2
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
    double diffusivity_intra;                       /*!< Di, intra-cellular diffusivity constant                                    */
    double diffusivity_extra;                       /*!< De, extra-cellular diffusivity constant                                    */
    double sim_duration;                            /*!< simulation total time                                                      */
    bool write_traj;                                /*!< flag, write a traj file or not, binary format only                         */
    bool write_hit;                                /*!< flag, write a hit file or not, binary format only                         */

    bool write_full_c;

    bool write_txt;                                 /*!< flag, writes DWI output signals in .txt if True                            */
    bool write_bin;                                 /*!< flag, writes the output signal in binary format (True by default)          */
    bool scale_from_stu;                            /*!< flag, true if the scheme file is in standar units m,s                      */
    bool save_phase_shift;                          /*!< flag, saves the phase shift distribution for all particles                 */
    long seed;                                      /*!< Initial seed for the random generator                                      */
    bool verbatim;                                  /*!< False to omit displaying state and warnings                                */
    std::string traj_file;                          /*!< Trajectory file path                                                       */
    std::string hit_file;                           /*!< Hit file path                                                       */
    std::string output_base_name;                   /*!< output files base name (path + sufix)                                      */
    std::string ini_walkers_file;                   /*!< initial walker position file (if any)                                      */
    unsigned ini_walkers_file_count;                /*!< number of walker positions initialize in the configuration file            */
    std::string ini_walker_flag;                    /*!< where to initialize the walkers                                            */
    std::string scheme_file;                        /*!< signal adquisition scheme file (if any)                                    */
    Eigen::Vector3d min_limits;                     /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits;                     /*!< voxel max limits (if any)                                                  */

    std::vector<std::string> spheres_files;         /*!< file paths with a list of spheres obstacles                              */
    std::vector<std::string> cylinders_files;       /*!< file paths with a list of cilinders obstacles                              */
    std::vector<std::string> PLY_files;             /*!< file paths with PLY obstacle files                                         */
    std::vector<double> PLY_scales;                 /*!< Auxiliary vector to save PLY file scales                                   */
    std::vector<float> ini_delta_pos;               /*!< Delta position for the wlakers                                             */

    unsigned num_proc;                              /*!< Number of precessors/process to launch in parallel                         */

    std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>>voxels_list; /*!< voxel min and max positions list   (deprecated)        */

    std::vector<Eigen::Vector3f> prop_dirs;         /*!< Saves the directions used to compute the propagator                        */

    std::vector<unsigned> record_pos_times;         /*!< time indexes, used to save the position of all walkers at certain time     */
    std::vector<unsigned> record_phase_times;       /*!< time indexes, used to save the phase shif of all walkers at certain time   */
    std::vector<unsigned> record_prop_times;        /*!< time indexes, used to save the mean propagator of the walkers at c. times  */

    bool   hex_packing;                             /*!< flag, true if an haxagonal packing should be used                          */
    double hex_packing_radius;                      /*!< float, constant radius for the cylinders                                          */
    double hex_packing_separation;                  /*!< float, separation distance betwen cylinders (separation > 2*radius)        */

    bool        gamma_packing;                      /*!< flag, true if a gamma distribution of cylinders will be initialized        */
    bool        gamma_output_conf;
    double      gamma_packing_alpha;
    double      gamma_packing_beta;
    double      gamma_icvf;
    double      gamma_output_configuration;
    unsigned    gamma_num_cylinders;
    
    /*
     Implementation of gamma sphere packing by Remy -
    */
    bool        gamma_packing_s;                      /*!< flag, true if a gamma distribution of spheres will be initialized        */
    bool        gamma_output_conf_s;
    double      gamma_packing_alpha_s;
    double      gamma_packing_beta_s;
    double      gamma_icvf_s;
    double      gamma_output_configuration_s;
    unsigned    gamma_num_spheres_s;

    /*
     Implementation of multi gamma sphere packing by Remy -
    */
    bool        gamma_packing_smul;                      /*!< flag, true if a multi gamma distribution of spheres will be initialized        */
    bool        gamma_output_conf_smul;
    std::vector<double>      gamma_packing_alpha_smul;
    std::vector<double>      gamma_packing_beta_smul;
    double      gamma_icvf_smul;
    double      gamma_output_configuration_smul;
    std::vector<unsigned>    gamma_num_spheres_smul;


    /*
     Implementation of uniform sphere packing by Remy -
    */
    bool                    uniform_packing_s;                      /*!< flag, true if a gamma distribution of spheres will be initialized        */
    bool                    uniform_packing_output_conf_s;
    double                  uniform_packing_icvf_s;
    double                  uniform_packing_output_configuration_s;
    std::vector<unsigned>   uniform_packing_num_spheres_s;
    std::vector<double>     uniform_packing_radii_s;

    /* 
    Implementation of multiple permeability by Remy
    */ 
    double obstacle_permeability     = -1.0;           /*!< Obstacles permeability if global                                          */
    std::vector<std::string> permeability_files;    /*!< Obstacles permeability file if local                                       */


    bool subdivision_flag           = false;        /*!< flag to check if we have several voxel subdivision to compute the signal   */
    unsigned number_subdivisions    = 0;            /*!< saves the number of subdivisions for an initialzied voxel (needed)         */
    std::string subdivisions_file   = "";           /*!< file with the list of subdivisions coordinates to compute the signal       */
    std::vector<Subdivision> subdivisions;          /*!< saves actual positions of the subdivision to compute the signal            */

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

    /*! \fn  setDiffusivities
     *  \param Diff_i intra-axonal diffusivity value.
     *  \param Diff_e extra-axonal diffusivity value
     *  \brief set the simulation diffusivities.
    */
    void setDiffusivity(double Diff_i, double Diff_e);

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

        /*! \fn setHitFlag
     *  \param write_bin, boolean flag (0,1).
     *  \brief sets the flag to save a binary hit file.
    */
    void setWriteHitFlag(bool write_bin);

        /*! \fn setFullFlag
     *  \param write_full_c_, boolean flag (0,1).
     *  \brief sets the flag to save the full colision file.
    */
    void setWriteFullFlag(bool write_full_c_);


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

    /*! \fn setHitFileName
     *  \param hit_file_ prefix of the hit file.
     *  \brief Set the prefix of the name for the hit file (txt and .hit)
    */
    void setHitFileName(std::string hit_file_);


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

    /*! \fn getDiffusivity_intra
     *  \return intra-axonal Diffusivity
    */
    double getDiffusivity_intra();

    /*! \fn getDiffusivity_extra
     *  \return extra-axonal Diffusivity
    */
    double getDiffusivity_extra();


    /*! \fn getWriteTrajFlag
     *  \return flag of the  binary traj file writer
    */
    bool getWriteTrajFlag();

        /*! \fn getWriteHitFlag
     *  \return flag of the  binary hit file writer
    */
    bool getWriteHitFlag();

    /*! \fn getWriteFullFlag
     *  \return flag of the  binary hit file writer
    */
    bool getWriteFullFlag();

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

    /*! \fn getHitFileName
     *  \return trajectory prefix
    */
    std::string getHitFileName();

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
    
    /*! \fn readGammaParams_s
     *  \param file input iostreams
     *  \brief reads the parameters needed to define an gamma distributed packing of spheres
    */
    void readGammaParams_s(std::ifstream& in);

    /*! \fn readGammaParams_smul
     *  \param file input iostreams
     *  \brief reads the parameters needed to define an multi gamma distributed packing of spheres
    */
    void readGammaParams_smul(std::ifstream& in);


    /*! \fn readUniformParams_s
     *  \param file input iostreams
     *  \brief reads the parameters needed to define a packing of spheres with fixed radius
    */
    void readUniformPackingParams_s(std::ifstream &in);

    /*! \fn readHexagonalParams
     *  \param file input iostreams
     *  \brief reads the subdivisions for computing the DW signal
    */
    void readSubdivisionFile();

    /*! \fn
     *  \brief adds the number of given subdivisions for the voxel
    */
    void addSubdivisions();


    /*! \fn
     *  \brief read the directions used to compute the propagator.
    */
    void readPropagatorDirections(std::string dir_path);

    /*! \fn readPermeability
     *  \param file input iostreams
     *  \brief reads the permeability on the configuration file.
    */
    void readPermeability(std::ifstream &in);


};

#endif // PARAMETERS_H
