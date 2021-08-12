#include "plyobstacle.h"
#include <fstream>
#include <iostream>
#include "constants.h"
#include "Eigen/Dense"

using namespace std;

PLYObstacle::PLYObstacle()
{
    file_path    = "";
    vert_number  = 0;
    face_number  = 0;
    vertices     = nullptr;
    faces        = nullptr;
    scale_factor = 1;
    percolation  = 0;
    count_perc_crossings = 0;
}

PLYObstacle::PLYObstacle(string path, double scale_factor_)
{
    file_path    = "";
    vert_number  = 0;
    face_number  = 0;
    vertices     = nullptr;
    faces        = nullptr;
    scale_factor = scale_factor_;
    percolation  = 0;
    count_perc_crossings = 0;
    readPLY_ASCII_triangles(path);
}

PLYObstacle::PLYObstacle(string path, std::vector<Eigen::Vector3d> &centers, double max_distance, double scale_factor_)
{
    file_path    = "";
    vert_number  = 0;
    face_number  = 0;
    vertices     = nullptr;
    faces        = nullptr;
    scale_factor = scale_factor_;
    percolation  = 0;
    count_perc_crossings = 0;
    readPLY_ASCII_trianglesSubdivitionDistance(path,centers,max_distance);
}


void PLYObstacle::readPLY_ASCII_triangles(std::string ply_file)
{
    if (vertices != nullptr)
        delete[] vertices;
    if (faces != nullptr)
        delete[] faces;

    std::ifstream in(ply_file.c_str(),std::ifstream::in);

    if(!in){
        std::cout << "Error opening file " << ply_file << std::endl;
        assert(1);
        return;
    }

    std::string tmp = "";
    while(tmp.compare("end_header")){
        in >> tmp;

        if(!tmp.compare("vertex")){
            in >> vert_number;
        }
        if(!tmp.compare("face")){
            in >> face_number;
        }
    }

    vertices = new Vertex[vert_number];
    faces = new Triangle[face_number];

    for (unsigned i =0; i< vert_number; i++){
        in >> vertices[i].points[0];
        in >> vertices[i].points[1];
        in >> vertices[i].points[2];

        //cout << vertices[i].points[0] << " " << vertices[i].points[1] << " " << vertices[i].points[2] << " " << endl;
    }


    for (unsigned i =0; i< vert_number; i++){
        vertices[i].points[0]*=scale_factor;
        vertices[i].points[1]*=scale_factor;
        vertices[i].points[2]*=scale_factor;
    }
    int num;
    for (unsigned i = 0; i < face_number; ++i) {
        in >> num;
        //in >> faces[i].index;
        in >> faces[i].indexes[0];
        in >> faces[i].indexes[1];
        in >> faces[i].indexes[2];
        faces[i].vertices = vertices;
        faces[i].saveNormalAndAuxInfo();
        //cout << faces[i].indexes[0] << " " << faces[i].indexes[1] << " "  << faces[i].indexes[2] << endl;
    }

}

void PLYObstacle::readPLY_ASCII_trianglesSubdivitionDistance(string ply_file, vector<Eigen::Vector3d>& centers, double max_distance)
{

    if (vertices != nullptr)
        delete[] vertices;
    if (faces != nullptr)
        delete[] faces;

    std::ifstream in(ply_file.c_str(),std::ifstream::in);

    if(!in){
        std::cout << "Error opening file " << ply_file << std::endl;
        assert(0);
        return;
    }

    std::string tmp = "";
    while(tmp.compare("end_header")){
        in >> tmp;

        if(!tmp.compare("vertex")){
            in >> vert_number;
        }
        if(!tmp.compare("face")){
            in >> face_number;
        }
    }

    vertices = new Vertex[vert_number];
    faces = new Triangle[face_number];

    for (unsigned i =0; i< vert_number; i++){
        in >> vertices[i].points[0];
        in >> vertices[i].points[1];
        in >> vertices[i].points[2];

        //cout << vertices[i].points[0] << " " << vertices[i].points[1] << " " << vertices[i].points[2] << " " << endl;
    }


    for (unsigned i =0; i< vert_number; i++){
        vertices[i].points[0]*=scale_factor;
        vertices[i].points[1]*=scale_factor;
        vertices[i].points[2]*=scale_factor;
    }

    int in_index = 0;
    double  distance;
    int num;
    for (unsigned i = 0; i < face_number; ++i) {
        in >> num;
        //in >> faces[i].index;
        in >> faces[in_index].indexes[0];
        in >> faces[in_index].indexes[1];
        in >> faces[in_index].indexes[2];
        faces[in_index].vertices = vertices;
        faces[in_index].saveNormalAndAuxInfo();

        if(centers.size()>0){
            for (auto c:centers ){
                //auto c= centers[j];
                distance = faces[in_index].minDistance(c);

                if(distance < max_distance){
                    in_index++;
                    break;
                }
            }
        }
        else{
            in_index++;
        }

        //cout << faces[i].indexes[0] << " " << faces[i].indexes[1] << " "  << faces[i].indexes[2] << endl;
    }

    cout << "before " << face_number << endl;
    face_number = in_index;
    cout << "after " << face_number << endl;

}


bool PLYObstacle::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision)
{
    Collision colision_temp;
    colision.type = Collision::null;
    colision.t = INFINITY_VALUE;

    //Origin O
    Eigen::Vector3d ray_origin,end_point;

    walker.getVoxelPosition(ray_origin);

    //To keep track of the closest collision
    double max_collision_distance = step_lenght;

    //En position in case of no collision.
    end_point = ray_origin + max_collision_distance*step;

    //For each triangle on the mesh model
    for (unsigned i=0; i < face_number; i++){
        faces[i].stepIntersects_MT(walker,step,max_collision_distance,colision_temp);
        handleCollisions(colision,colision_temp,max_collision_distance,end_point,i);
    }

    updateWalkerStatusAndHandleBouncing(walker,ray_origin,step,colision);

    if(colision.type == Collision::null){
        return false;
    }

    return true;

}

bool PLYObstacle::checkCollision(Walker &walker, Eigen::Vector3d &step, double &step_lenght, Collision &colision, std::vector<unsigned>& triangle_list, unsigned list_end)
{
    Collision colision_temp;
    colision.type = Collision::null;
    colision.t = INFINITY_VALUE;

    //Origin O
    Eigen::Vector3d ray_origin,end_point;

    walker.getVoxelPosition(ray_origin);

    //To keep track of the closest collision
    double max_collision_distance = step_lenght;

    //En position in case of no collision.
    end_point = ray_origin + max_collision_distance*step;

    //For each triangle on the mesh model
    for (unsigned i=0; i < list_end; i++){
        unsigned triangle_index = triangle_list[i];
        faces[triangle_index].stepIntersects_MT(walker,step,max_collision_distance,colision_temp);
        handleCollisions(colision,colision_temp,max_collision_distance,end_point,triangle_index);
    }

    updateWalkerStatusAndHandleBouncing(walker,ray_origin,step,colision);

    if(colision.type == Collision::null){
        return false;
    }

    return true;
}

void PLYObstacle::handleCollisions(Collision &colision_confirmed, Collision &colision_2, double &max_distance,  Eigen::Vector3d &end_point, const unsigned triangle_indx)
{

    // nothing to do;
    if (colision_2.type == Collision::null)
        return;

    //WARNING: Cuidar este patch
    // Implementa Percolacion
    if(percolation>0.0){
        double _percolation_ ((double)rand()/RAND_MAX);

        if( percolation - _percolation_ > EPS_VAL ){
            count_perc_crossings++;
            return;
        }
    }

    colision_2.triangle_ind = triangle_indx;

    if (colision_confirmed.type == Collision::hit || colision_confirmed.type == Collision::boundary){
        if(colision_2.doIHaveMorePiorityThan(colision_confirmed)){
            colision_confirmed = colision_2;
            colision_confirmed.triangle_ind = triangle_indx;
        }
        return;
    }

    if(colision_confirmed.type == Collision::near ){
        if (colision_2.type == Collision::hit || colision_2.type == Collision::boundary){
            colision_confirmed = colision_2;
            //max_distance = colision_2.t;
            colision_confirmed.triangle_ind = triangle_indx;
        }
        return;
    }

    // if we get here means that colision_confirmed.type = 'null'
    if(colision_2.type == Collision::near){

        checkIfItsNearToTriangle(end_point,triangle_indx,colision_2);

        // if we were near indeed
        if(colision_2.type != Collision::null){
            colision_confirmed = colision_2;
            colision_confirmed.triangle_ind = triangle_indx;
        }
        return;
    }

    colision_confirmed = colision_2;
}

void PLYObstacle::checkIfItsNearToTriangle(const Eigen::Vector3d end_point, const unsigned triangle_ind, Collision &colision)
{
    double EPS = 2e-10;
    double u,v,t;
    bool hit = faces[triangle_ind].rayIntersects_MT(end_point,faces[triangle_ind].normal,u,v,t);

    if( hit && fabs(t) <= EPS){
        colision.t = fabs(t);
        colision.u = u;
        colision.v = v;
        colision.triangle_ind = triangle_ind;
    }
    else{
        colision.type = Collision::null;
    }
}

bool PLYObstacle::updateWalkerStatusAndHandleBouncing(Walker &walker, Eigen::Vector3d &ray_origin, Eigen::Vector3d &step, Collision &colision)
{
    if (colision.type == Collision::null){
        return 0;
    }

    //If the particle bounced
    bool bounced=false;

    // we set the status of the walker ( on_triangle, on_vertex, etc)
    colision.computeCollisionLocation();

    //If was a hit and need to bounce;
    if(colision.type == Collision::hit){
        bounced = true;
        if (colision.col_location == Collision::on_edge || colision.col_location == Collision::on_vertex){
            colision.bounced_direction = -step;
            //WARNING: REMOVE THIS
           // cout << "Edge" << endl;
        }
        else{
            Eigen::Vector3d normal;
            faces[colision.triangle_ind].getNormal(normal);

            Eigen::Vector3d temp_step = step;
            elasticBounceAgainsPlane(ray_origin,normal,colision.t,temp_step);

            colision.bounced_direction = temp_step.normalized();

            //Orientation respect the triangle
            double dot = ((walker.pos_v - faces[colision.triangle_ind].center).normalized()).dot(normal);

            colision.col_location = (dot < -1e-5)?Collision::inside:(dot > 1e-5)?Collision::outside:Collision::unknown;
        }
    }
    else if(colision.type == Collision::near){
        bounced = false;
    }

    return bounced;
}

double PLYObstacle::minDistance(Walker &w, unsigned t)
{
    return faces[t].minDistance(w.pos_v);
}



