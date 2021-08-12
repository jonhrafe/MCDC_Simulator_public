//!  Spin Final class =============================================================/
/*!
  Basic unit of the diffusion process.

  \author Jonathan Rafael
  \date November 2016
  \version 0.2
====================================================================================*/

#ifndef WALKER_H
#define WALKER_H

#include "Eigen/Core"
#include <vector>
#include <deque>
#include "collisionsphere.h"
#include <iostream>

/*! \class Walker
 *  \brief Alias for a particle. Basic unit on the simulation process. Saves all the necessary information to perform
 *         the particles dynamics.
 */
class Walker {
public:
    //! An enum.
    /*! All the possibles states that a walker can be in a given step.
        The next step is perform according to this state
    */
    enum state {on_object, on_edge, on_vertex, on_voxel, free, bouncing};

    //! An enum.
    /*! Possible location of the walker inside the voxel.
        Checks illegal crossings of the barrier (border, lol)
    */
    enum RelativeLocation{unknown,intra,extra};

    Eigen::Vector3d pos_r;                                          /*!< Real walker position for collision, r stands for real                  */

    Eigen::Vector3d pos_v;                                          /*!< Walker current position                                                */

    Eigen::Vector3d last_pos_r;                                     /*!< Walker voxel last position                                             */

    Eigen::Vector3d last_pos_v;                                     /*!< Walker real last position                                              */

    Eigen::Vector3d ini_pos;                                        /*!< Walker intital position                                                */

    Eigen::Vector3d next_direction;                                 /*!< Auxiliar vector for special states cases, decides the next direction   */

    Eigen::Matrix3Xd pos_r_log;                                     /*!< log of the real spin position, used to compute the phase shift         */

    Eigen::Matrix3Xd pos_v_log;                                     /*!< log of the voxel position, used for collision location and bouncing    */

    int in_obj_index;                                               /*!< Auxiliar index to save if the walker was inside a convex object        */

    int in_ply_index;                                               /*!< Auxiliar index to save if the walker was inside a convex ply object    */

    int in_sph_index;                                               /*!< Auxiliar index to save if the walker was inside a sphere               */

    ObstacleCollisionSphere cylinders_collision_sphere;             /*!< Collision sphere for collition against cylidners                       */

    ObstacleCollisionSphere spheres_collision_sphere;               /*!< Collision sphere for collition against cylidners                       */

    PLYCollisionSphere ply_collision_sphere;                        /*!< Collision sphere for collition against PLY meshes                      */

    Eigen::Vector3d initial_sphere_pos_v;                           /*!< Saves the intial positioon of the walker inside the collition sphere   */

    unsigned steps_count;                                           /*!< Counts the number of steps (including bouncings) made.                 */

    state status;                                                   /*!< state memeber */

    RelativeLocation initial_location, location;                    /*!< location on the substrate (if known)*/

    int intra_extra_consensus;                                      /*!< intra o extra position by face collision consensus w/r the normal*/

    unsigned intra_coll_count;                                      /*!< counter of collision in the ïntra-side w/r the normal*/

    unsigned extra_coll_count;                                      /*!< counter of collision in the extra-side w/r the normal*/

    unsigned int index;                                             /*!< Walker identifier (id)*/

    unsigned int rejection_count;                                   /*!< counter of the rejected directions in a single time-step*/

    float steps_per_second;                                         /*!< Particles steps per second speeed.*/

    //! Default constructor.
    /*! Set all variables to cero.*/
    Walker();

    //! Default destructor.
    //!
    /*! Does nothing
    */
    ~Walker() {}

    //! Constructor.
    /*! Initialize the walker position in a random position inside the boundaries
        defined by the limits.
        \param xmin lower x threshold
        \param xmax upper x threshold
        \param ymin lower y threshold
        \param ymax upper y threshold
    */
    Walker(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);

    // Get methods
    void  getRealPosition(double &, double &, double &) const;
    void  getRealPosition(Eigen::Vector3d&) const;
    void  getVoxelPosition(double &, double &, double &) const;
    void  getVoxelPosition(Eigen::Vector3d&) const;
    void  getInitialPosition(double &, double &, double &) const;
    void  getInitialPosition(Eigen::Vector3d&) const;
    void  getNextDirection(Eigen::Vector3d&) const;
    unsigned int getIndex() const;

    // Set methods
    void  setRealPosition(const double &, const double &,const double &);
    void  setRealPosition(const Eigen::Vector3d&);
    void  setVoxelPosition(const double &, const double &,const double &);
    void  setVoxelPosition(const Eigen::Vector3d&);
    void  setInitialPosition(const double &, const double &, const double &);
    void  setInitialPosition(const Eigen::Vector3d &);
    void  setNextDirection(Eigen::Vector3d &);
    void  setRandomInitialPosition(const Eigen::Vector3d &min, const Eigen::Vector3d &max);
    void  setIndex(unsigned int&);

    void setRealPosLog(const Eigen::Vector3d &pos,unsigned t);
    void setRealPosLog(double x, double y, double z, unsigned t);
    void setVoxPosLog(const Eigen::Vector3d &pos,unsigned t);
    void setVoxPosLog(double x, double y, double z, unsigned t);


    void setNumberOfSteps(unsigned T);

};


#endif // WALKER_H
