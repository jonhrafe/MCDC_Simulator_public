//!  Auxiliary class to save scheme_files values  =============================================================/
/*!
*   Helper class to store, handle and read scheme files values .
*   \author Jonathan Rafael
*   \date      November 2016
*   \version   0.2
*=============================================================================================================*/

#ifndef SCHEME_H
#define SCHEME_H
#include <vector>
#include <string>

class Scheme
{
public:
    std::string scheme_file; /*!< Scheme file path                      */
    std::string header;      /*!< Header on the scheme_file             */
    std::string type;        /*!< Sequence type (PGSE only so far)      */
    int num_rep;             /*!< Number of gradients                   */

    float duration;          /*!< time duration (wavefroms)             */
    float T;                 /*!< number of time steps (wavefroms)      */
    bool scale_from_stu;     /*!< True if the input is in standar units */


    std::vector< std::vector<double> > scheme; /*!< Scheme values       */

    /*! \fn  Default constructor.
     *  \brief Intialize everything with 0's
     */
    Scheme();

    /*! \fn  Constructor.
     * \param scheme file location
     *  \brief reads all the parameters listed in the scheme file
     */
    Scheme(std::string scheme_file_);

    ~Scheme(){}

    /*! \fn  read the scheme file list
     *  \brief Intialize and reads all the parameters listed on the scheme file
     */
    void readSchemeFile(std::string scheme_file_,bool scale_from_stu = 0);


private:
    void readPGSE(std::ifstream &in,bool scale_from_stu);
    void readAPGSE(std::ifstream &in,bool scale_from_stu);
    void readWaveForm(std::ifstream &in,bool scale_from_stu);

};

#endif // SCHEME_H
