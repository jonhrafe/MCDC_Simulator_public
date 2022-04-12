//!  CylinderDistribution Class =============================================================/
/*!
*   \details   Base class to construct a substrate made of cylinders from a distribution of radiis placed in
*              a single voxel structure.
*   \author    Remy Gardier
*   \date      April 2022
*   \version   0.0
=================================================================================================*/

#ifndef CYLINDERDISTRIBUTION_H
#define CYLINDERDISTRIBUTION_H

#include "Eigen/Core"
#include <vector>
#include "constants.h"
#include "simerrno.h"
#include <iostream>
#include "cylinder.h"
#include "obstacledistribution.h"

class CylinderDistribution : public ObstacleDistribution
{
public:

    std::vector<Cylinder> cylinders;                    /*!< Spheres vector                                                             */

    /*!
     *  \brief Default constructor. Does nothing
     */
    CylinderDistribution(); 

    /*!
     *  \param icvf_ Intra-cellular volume fraction
     *  \param min_l voxel min limits
     *  \param max_l voxel min limits
     *  \brief Initialize everything.
     */
    CylinderDistribution(double ,Eigen::Vector3d &, Eigen::Vector3d &, std::vector<double> &);


    /*!
     *  \brief Samples and constructs a distribution
    */
    void createSubstrate();

    /*!
     *  \brief Prints the sphere positions in a file or output stream.
     *  \param out ostream where to write the info.
    */
    void printSubstrate(std::ostream& out); 

private:

    /*!
     *  \brief Checks for collision between inside a voxel (with periodic boundaries)
     *  \param cyl cylinder to check collision with
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param cylinder_list cylinders already added.
     *  \param min_distance that two cylinders can be close to.
    */
    bool checkForCollition(Cylinder cyl, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits, std::vector<Cylinder>& cylinder_list, double &min_distance);

    /*!
     *  \brief Auxiliary function to check the BOundary collision
     *  \param cyl cylinder to check collision with.
     *  \param min_limits Voxel min limits.
     *  \param max_limits Voxel max limits.
     *  \param cylinder_list cylinders already added.
    */
    void checkBoundaryConditions(Cylinder cyl, std::vector<Cylinder>& cylinder_list, Eigen::Vector3d min_limits, Eigen::Vector3d max_limits);

    /*!
     *  \brief Computes Intra Celular Volum Fraction given the voxel limits and the list of added cylinders.
     *  \param spheres List of included cylinders.
     *  \param min_limits voxel min limits.
     *  \param max_limits voxel max limits.
    */
    double computeICVF(std::vector<Cylinder> &cylinders, Eigen::Vector3d &min_limits, Eigen::Vector3d &max_limits, int &num_no_repeat);

    /*!
     *  \brief Estimate the minimal voxel side length from Intra Celular Volum Fraction and the radiis of the obstacles
     *  \param radiis List of radiis.
     *  \param icvf_ Intra-cellular volume fraction
     *  \param l voxel side length.
    */
    void computeMinimalSize(std::vector<double> radiis, double icvf_, Eigen::Vector3d& l);

};

#endif // CYLINDERDISTRIBUTION_H
