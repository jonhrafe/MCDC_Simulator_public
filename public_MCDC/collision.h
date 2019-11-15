//!  Collision Final class =====================================================================/
/*!
    \details    Class to save and handle collisions between walkers and objects.
    \author     Jonathan Rafael
    \date       November 2016
===============================================================================================*/

#ifndef COLLISION_H
#define COLLISION_H
#include <Eigen/Core>

/*! \class Collision
 *  \brief Class to save and handle collisions between walkers and objects.
 *
 *  This class should handle all the cases where a collision can happened as well as all the operation
 *  between collision.
 *
 */
class Collision{
public:
    //! \enum  collision_type.
    /*! All the possibles cases or situations where a step can end.
     *  The next step is performed according to this state
    */
    enum collision_type{hit,near,null,boundary,degenerate};

    //! \enum collision_location.
    /*! Only in case of collision (or a very close ending position) this are  the cases where the
     * collision can happened.
    */
    enum collision_location{inside, on_edge, on_vertex, voxel, outside,unknown};

    collision_type type;                /*!< Saves the type of collsion (if any)                    */
    collision_location col_location;    /*!< Save the colocation of the collision over the object   */
    Eigen::Vector3d colision_point;     /*!< Saves the position of colision*/
    Eigen::Vector3d bounced_direction;  /*!< Save the bounced direction for a given obstacle        */

    double rn;                          /*!< saves the local orietnation between the wall and the
                                        particle                                                    */
    double u;                           /*!< u position in baricentric coordinates                  */
    double v;                           /*!< v position in baricentric coordinates                  */
    double t;                           /*!< signed, collision distance                             */

    int triangle_ind;                   /*!< In case of a PLY obstacle saves the triangle index. t
                                        collison distance                                           */
    int obstacle_ind;                   /*!< In case of a generic obstacle saves the obstacle index.*/

    /*! \fn  Default constructor.
     *  \brief Initialize everything with 0's and NULL states, the triangle and object indexes are set
     *  to -1.
     */
    Collision():u(0),v(0),t(1e15),triangle_ind(-1),obstacle_ind(-1){type=null;col_location=unknown;}

    /*! \fn  Default constructor.
     *  \brief Initialize everything with 0's and NULL states, the triangle and object indexes are
     *  set to -1.
     *  \param u is for parametric coordinates of a plane or triangle.
     *  \param u is for parametric coordinates of a plane or triangle.
     *  \param t is the collision distance,
     */
    Collision(double u_,double v_,double t_):u(u_),v(v_),t(t_),triangle_ind(-1),obstacle_ind(-1){}

    //! \fn Default destructor.
    /*! \brief Does nothing.
    */
    ~Collision();

    //! \fn Boolean comparison function
    //! \param coll collision to compare with.
    /*! \brief Compares a second collision to determine which one has more priority.
     *  The comparison is based on the type of collision and distance.
    */
    bool doIHaveMorePiorityThan(Collision &coll);

    //! \fn computeCollisionLocation();
    //! Auxiliar function for Triangular barriers
    /*! \brief Computes, based on the the coordinates u,v,t of the collision, the location relative
     *  to the triangle.
    */
    void computeCollisionLocation();

};


#endif // COLLISION_H
