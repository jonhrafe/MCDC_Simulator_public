//!  Auxiliary Subdivision Class  =============================================================/
/*!
  Auxiliary Class.
  Implementation of the subdivision of a voxel into separate adquisitions

  \date   September 2017
  \author Jonathan Rafael
  \version 0.1.0
*=====================================================================================*/


#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include "Eigen/Core"

class Subdivision
{
public:


    /**
     * @brief Naive constructor
     */
    Subdivision();


    /**
     * @brief Constructor for a defined list of min and max positions.
     */
    Subdivision(Eigen::Vector3f&,Eigen::Vector3f&);

    Eigen::Vector3f min_limits;                         /*!< Vector with the list of min limits points of each subdivisions                 */

    Eigen::Vector3f max_limits;                         /*!< Vector with the list of max limits points of each subdivisions                 */

    int density;                                        /*!< Counter to save the number of particles inside that region                      */

    int density_intra;                                  /*!< Counter to save the number of particles labeled as Intra  in that region        */

    int density_extra;                                  /*!< Counter to save the number of particles labeled as Extra in that region        */

    /**
     * @param pos: 3d position
     * @brief Auxiliary function to check if a 3d position is inside a "subdivision" i.e. defined cube
     */
    bool isInside(Eigen::Vector3d &pos);
};

#endif // SUBDIVISION_H
