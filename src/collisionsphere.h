//!  Collision Final class ============================================================================/
/*!
    \details    Class to implement spherical bounding boxes for the WALKER mean diffusion. This class
                provides methods in order to create and update spherical bounding boxes used to compute
                the collisions.
    \author     Jonathan Rafael
    \date       February 2017
======================================================================================================*/

#ifndef COLLISIONSPHERE_H
#define COLLISIONSPHERE_H
#include <vector>

/*! \class Collisionspheren
 * \brief   Father class. this class provides methods in order to create and update
 *          spherical bounding box used to compute the collisions. The implementation
 *          is based on two collision spheres. The inner one (small) and the (outer).
 *          The fist saves the objects where th particle MAY collide in a given time,
 *          While the second saves the full set of obstacles where the particle can
 *          possibly collide, i.e. that are physically achievable for the walker to
 *          collide.
*/
class Collisionsphere
{
public:

    float big_sphere_distance;          /*!< Size of the big (outer) collision sphere                                                                               */
    float small_sphere_distance;        /*!< Size of the small (inner) collision sphere                                                                             */

    unsigned list_size;

    Collisionsphere():big_sphere_distance(0),small_sphere_distance(0),list_size(0){}
};

/*! \class CylinderCollisionSphere
 *  \brief Class to save the cylinderical obstacles that a can collide to a walker.
 */
class ObstacleCollisionSphere: public Collisionsphere{

public:

    unsigned small_sphere_list_end;                         /*!< Index of the LAST element on the list for the small collision sphere                               */
    unsigned big_sphere_list_end;                           /*!< Index of the LAST element on the list for the big collision sphere                                 */

    std::vector<unsigned>* collision_list;                  /*! <Pointer to  List with the cylinders indexes. The indexes are permuted in its position.             */

    ObstacleCollisionSphere();

    //! \fn Removes an index from the inner sphere
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         then decrease the inner sphere end index. This way this index is no longer considered inner collision list.
    */
    void popFromSmallSphere(unsigned i);

    //! \fn Adds one element to the list by moving it in front of the current index and increasing the index.
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         then increase the inner sphere end index. This way this index is now included in the inner collision list.
    */
    void pushToSmallSphere(unsigned i);

    //! \fn Removes one element to the list by moving it in front of the current index and decreasing the index.
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         Then decrease the inner sphere end index. This way this index is now excluded in the outer collision list.
    */
    void popFromBigSphere(unsigned i);

    //! \fn Push one element to the list by moving it in front of the current index and increase the end index.
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         Then increase the inner sphere end index. This way this index is now included in the outer collision list.
    */
    void pushToBigSphere(unsigned i);

    //! \fn Set function to fix the outer sphere size.
    /*! \param size of the list
    */
    void setBigSphereSize(float size);

    //! \fn Set function to fix the inner sphere size.
    /*! \param size of the list
    */
    void setSmallSphereSize(float size);

    //! \fn Push one element to the complete obstacle list
    /*! \param element value to be added to the obstacle list
    */
    void push_index(unsigned int element);

};

/*! \class CylinderCollisionSphere
 *  \brief Class to save the PLY mehses and the subset of triangles that a can collide to a walker
 *
 */
class PLYCollisionSphere: public Collisionsphere{

public:

    std::vector<unsigned> small_sphere_list_end;                /*!< Index vector of the LAST element on the list for the small collision sphere                           */
    std::vector<unsigned> big_sphere_list_end;                  /*!< Index vecotr of the LAST element on the list for the big collision sphere                             */

    std::vector<std::vector<unsigned>>* collision_list;         /*!< Pointer to the list with the triangle indexes for each PLY. The indexes are permuted in its position. */

    PLYCollisionSphere();

    //! \fn Removes an index from the inner sphere
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         then decrease the inner sphere end index. This way this index is no longer considered inner collision list.
    */
    void popFromSmallSphere(unsigned i,unsigned t);             /*!< Removes one index from the list by moving it to the end of the list and decreading the index.         */


    //! \fn Adds one element to the list by moving it in front of the current index and increasing the index.
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         then increse the inner sphere end index. This way this index is now included in the inner collision list.
    */
    void pushToSmallSphere(unsigned i,unsigned t);              /*!< Adds one element to the list by moving it in front of the current index and increasing the index.     */

    //! \fn Removes one element to the list by moving it in front of the current index and decreasing the index.
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         Then decrease the inner sphere end index. This way this index is now excluded in the outer collision list.
    */
    void popFromBigSphere(unsigned i,unsigned t);

    //! \fn Push one element to the list by moving it in front of the current index and increase the end index.
    //! \param index in the collision list. Notice that is the index, no the value.
    /*! \brief This function receives a index from the collision list and moves the value to the last position of the list.
     *         Then increase the inner sphere end index. This way this index is now included in the outer collision list.
    */
    void pushToBigSphere(unsigned i,unsigned t);

    //! \fn Set function to fix the outer sphere size.
    /*! \param size of the list
    */
    void setBigSphereSize(float size);

    //! \fn Set function to fix the inner sphere size.
    /*! \param size of the list
    */
    void setSmallSphereSize(float size);

    //! \fn Push one element to the complete obstacle list
    /*! \param element value to be added to the obstacle list
    */
    void push_ply(std::vector<unsigned> list);

};

#endif // COLLISIONSPHERE_H
