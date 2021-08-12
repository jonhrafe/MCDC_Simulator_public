//!  Simulation Input and parameter errors handling class =================================================/
/*!
    \details    Class ot handle the errors in the parameters, logical and on the syntaxis.
    \author     Jonathan Rafael
    \date       March 2016

*=========================================================================================================*/
#ifndef SIMERRNO_H
#define SIMERRNO_H
#include <assert.h>
#include <string>
#include <fstream>
#include "parameters.h"


/*! \class SimErrno
 *  \brief Class to handle the errors in the parameters, logical and on the syntax
 *
 *  This class contains a set of static methods to check that the configuration files exist and
 *  that the parameters are correctly set. This may cause asserts ERRORS or WARNINGS.
 *
 */
class SimErrno
{
public:
    SimErrno();

    //! \fn check if a given file exists
    /*! \param name file path
     *  \brief Return true if the file does exist, false otherwise.
    */
    inline static bool checkFileExist(const std::string name){
        std::ifstream f(name.c_str());
          return f.good();
    }

    //! \fn checks if the given parameters make sense
    /*! \param parameter instance
     *  \brief Return false if any of the parameters are inconsistent or bugged. In may assert the
     *  program.
    */
    static bool checkSimulationParameters(Parameters& params);

    //! \fn checks if the given scheme files exist and make sense
    /*! \param parameter instance
     *  \brief Return false if any of the parameters are inconsistent or bugged. In may assert the
     *  program.
    */
    static bool checkSchemeFile(Parameters& params);

    //! \fn checks if the given PLY file make sense
    /*! \param parameter instance
     *  \brief Return false if any of the PLY files are inconsistent or bugged. In may assert the
     *  program.
    */
    static bool checkPLYFiles(Parameters& params);

    //! \fn checks if the given cylinder list file make sense
    /*! \param parameter instance
     *  \brief Return false if any of the cylinder list files are inconsistent or bugged.
     *  In may assert the program.
    */
    static bool checkCylindersListFile(Parameters &params);

    //! \fn checks if the given cylinder list file make sense
    /*! \param parameter instance
     *  \brief Return false if any of the cylinder list files are inconsistent or bugged.
     *  In may assert the program.
    */
    static bool checkSphereListFile(Parameters &params);

    //! \fn checks if the given initial walker file make sense
    /*! \param parameter instance
     *  \brief Return false if the initial position file is inconsistent or bugged. In may assert the
     *  program.
    */
    static bool checkInitWalkerFile(Parameters &params);

    //! \fn checks if the given voxel limits make sense
    /*! \param parameter instance
     *  \brief Return false if the voxel instances are inconsistent or bugged. In may assert the
     *  program.
    */
    static bool checkVoxelLimits(Parameters& params);

    //! \fn checks if the overall configuration file make sense
    /*! \param parameter instance
     *  \brief Return false if the scheme file does not exist or there are inconsistent or bugs.
     * In may assert the program.
    */
    static bool checkConfigurationFile(const char *configuration_file);

    //! \fn Prints the simulation over all info and parameters.
    /*! \param parameter instance
     *  \param iostream to print to
     *  \param color flag, false if no colour should be display or written
    */
    static void printSimulatinInfo(Parameters &params, std::ostream&, bool color = 1);


    //! \fn checks if the overall output location make sense
    /*! \param parameter instance
     *  \brief Return false if the output location and prefix are inconsistence or bugged
    */
    static void checkOuputPrefixAndWriteInfo(Parameters &params);

    //! \fn checks if the Gamma distribution parameters are correct
    /*! \param parameter instance
     *  \brief Return false if the there are errors or inconsistencies in the gamma distr. parameters
    */
    static bool checkGammaDistributionParamaters(Parameters &params);

    //! \fn prints a warning message to console or a file.
    //! \param message to be displayed or written
    /*! \param iostream where to print
    *  \param colour flag, false if no colour should be display or written
    */
    static void warning(std::string message, std::ostream &, bool color = 1);

    //! \fn prints a info message to console or a file.
    //! \param message to be displayed or written
    /*! \param iostream where to print
    *   \param colour flag, false if no colour should be display or written
    */
    static void info(std::string message, std::ostream &, bool color = 1);

    //! \fn prints a info in menu format message to console or a file.
    //! \param message to be displayed or written
    /*! \param iostream where to print
    *!  \param colour flag, false if no colour should be display or written
    *!  \param spacing at the end of the message
    */
    static void infoMenu(std::string message, std::string value, std::ostream &, bool color = 1, int space = 0);

    //! \fn prints an error message to console or a file.
    //! \param message to be displayed or written
    /*! \param iostream where to print
    *!  \param colour flag, false if no colour should be display or written
    */
    static void error(std::string message, std::ostream &, bool color = 1);

    //! \fn prints the completed percentage and the expected time in a fixed format.
    //! \param completed string with the completed percentage
    //! \param time string with the remaining (expected) time
    /*! \param iostream where to print
    *!  \param colour flag, false if no colour should be display or written
    *!  \param end flag, false if no end of line string should be printing
    */
    static void expectedTime(std::string completed, std::string time, std::ostream &, bool color = 1, std::string steps_second="", std::string endl_str ="");

    //! \fn Prints the current Date and time in the format: DD-MM-YYYY (HH:mm:ss) .
    /*!
    */
    static std::string currentDateTime();

    //! \fn checks if the given subdivision file is in the correct format
    /*! \param parameter instance
     *  \brief Return false if any of the elements in the file are miss configured
    */
    static bool checkSubdivisionsFile(Parameters & params);

    //! \fn appends a repetition label on the prefix command so no results are overwritten
    /*! \param parameter instance
     *  \brief Appends a repetition label on the prefix command so no results are overwritten, helpful if you are
     *         running batch of simulation inside a server.
    */
    static void appendRepetitionLabel(Parameters& params);

};


#endif // SIMERRNO_H
