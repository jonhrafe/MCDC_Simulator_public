//!  SphereGammaDistributionMul Class =============================================================/
/*!
*   \details   Class to construct a substrate taken from Multi Gamma distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      March 2022
*   \version   0.0
=================================================================================================*/

#ifndef SPHEREGAMMADISTRIBUTIONMUL_H
#define SPHEREGAMMADISTRIBUTIONMUL_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "sphere.h"

class SphereGammaDistributionMul
{
public:

    std::vector<unsigned> num_spheres;                             /*!< number of spheres fit inside the substrate                                */
    std::vector<double> alpha;                                   /*!< alpha coefficient of the Gamma distribution                                */
    std::vector<double> beta;                                    /*!< beta coefficient of the gamma distribution                                 */
    double icvf;                                    /*!< Achieved intra-celular volum fraction in the substrate                     */

    Eigen::Vector3d min_limits;                     /*!< voxel min limits (if any) (bottom left corner)                             */
    Eigen::Vector3d max_limits;                     /*!< voxel max limits (if any)                                                  */
    std::vector<Sphere> spheres;                /*!< Spheres vector                                                             */

    /*!
     *  \param P_ Sphere origin
     *  \param radius_ sphere's radius
     *  \param scale scale factor for the values passed. Useful when reading a file.
     *  \brief Initialize everything.
     */
    SphereGammaDistributionMul(std::vector<unsigned> ,std::vector<double>,std::vector<double>, double,Eigen::Vector3d &, Eigen::Vector3d &);


    // void setVariables(unsigned num_sph, double a, double b,double icvf_,Eigen::Vector3d & min_l, Eigen::Vector3d &max_l);


    /*!
     *  \brief Shows a small histogram of the gamma distribution
    */
    void displayGammaDistribution();

    /*!
     *  \brief Samples and constructs a Gamma distribution
    */
    void createGammaSubstrate();

    /*!
     *  \brief Prints the sphere positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out); 

private:

    /*!
     *  \brief Checks for collision between inside a voxel (with periodic boundaries)
     *  \param sph sphere to check collision with
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param spheres_list spheres already added.
     *  \param min_distance that two spheres can be close to.
    */
    bool checkForCollition(Sphere sph, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector<Sphere>& spheres_list, double &min_distance);

    /*!
     *  \brief Auxiliary function to check the BOundary collision
     *  \param sph sphere to check collision with.
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param spheres_list spheres already added.
    */
    void checkBoundaryConditions(Sphere sph, std::vector<Sphere>& spheres_list, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits);

    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added spheres.
     *  \param spheres List of included spheres.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
    */
    double  computeICVF(std::vector<Sphere> &spheres, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);

    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);


};

#endif // SPHEREGAMMADISTRIBUTIONMUL_H
