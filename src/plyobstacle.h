//!  PlyObstacle Derived Class ====================================================================/
/*!
*   \details   PLYObstacle derived class. Implements obstacles loaded from pre-defined PY meshes.
*   \author    Jonathan Rafael
*   \date      November 2016
*   \version   0.2
*=================================================================================================*/


#ifndef PLYOBSTACLE_H
#define PLYOBSTACLE_H

#include "obstacle.h"
#include "triangle.h"

/*! \class  PLYObstacle
 *  \brief  Implements obstacles loaded from pre-constructed PLY meshes. The PLY format should be without any other
 *          experiment.
 */
class PLYObstacle : public Obstacle
{
public:

    unsigned    vert_number;
    unsigned    face_number;
    std::string file_path;
    Vertex*     vertices;
    Triangle*   faces;
    double      scale_factor;
    int         id;

    PLYObstacle();
    PLYObstacle(std::string path,double scale_factor_ = 1);
    PLYObstacle(std::string path, std::vector<Eigen::Vector3d> &centers, double max_distance=INFINITY,double scale_factor_ = 1);


    void readPLY_ASCII_triangleFan(std::string ply_file);
    void readPLY_ASCII_triangles(std::string ply_file);
    void readPLY_ASCII_trianglesSubdivitionDistance(std::string ply_file, std::vector<Eigen::Vector3d> &centers, double max_distance);

    void setScaleFactor(double scale){scale_factor = scale;}

//  bool computeStepCollition(Walker &w, double step[3], const double &step_length,double end_point[3], Collision& colision);
    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision);

    bool checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision, std::vector<unsigned> &triangle_list, unsigned list_end);

    double minDistance(Walker& w, unsigned t);

private:
    // Compare 2 coliision and decides wich one has the highest piority and if
    // colision_2 neess to be handled differently
    void handleCollisions(Collision& colision_confirmed, Collision& colision_2, double& max_distance, Eigen::Vector3d &end_point, const unsigned triangle_indx);

    //Function to check if a point is close to a a certain triangle. save the result in a Collision object
    void checkIfItsNearToTriangle(const Eigen::Vector3d end_point, const unsigned triangle_ind, Collision &colision);

    //Given the collision, handles the next walker status and the bouncing, if needed.
    bool updateWalkerStatusAndHandleBouncing(Walker &walker, Eigen::Vector3d &ray_origin, Eigen::Vector3d &step, Collision &colision);

};

#endif // PLYOBSTACLE_H
