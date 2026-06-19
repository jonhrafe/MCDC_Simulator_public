#include "plyobstacle.h"
#include <fstream>
#include <iostream>
#include "Eigen/Dense"
#include "simerrno.h"

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
    
    // Check if the PLY file is binary format
    if (isPLYBinary(path)) {
        readPLY_Binary(path);
    } else {
        readPLY_ASCII_triangles(path);
    }
    
    createAABBs();
    //Todo make a dynamic size for the grid

    double optimal_cell_size = this->AABBgrid.computeOptimalCellSize(this->aabbs,AABB_memory_limit_mb, min_cell_size_um);
    AABBgrid.InitializeGrid(this->aabbs,optimal_cell_size);
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
    
    // Check if the PLY file is binary format
    if (isPLYBinary(path)) {
        readPLY_Binary_trianglesSubdivitionDistance(path, centers, max_distance);
    } else {
        readPLY_ASCII_trianglesSubdivitionDistance(path, centers, max_distance);
    }
    
    createAABBs();

    double optimal_cell_size = this->AABBgrid.computeOptimalCellSize(this->aabbs,AABB_memory_limit_mb, min_cell_size_um);

    if(optimal_cell_size < 100){
        std::string message = "Spheres' grid size: " + std::to_string(optimal_cell_size*1000) + " um";
        SimErrno::info(message,cout);
    }

    AABBgrid.InitializeGrid(this->aabbs,optimal_cell_size);
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
        //cout << faces[i].indexes[0] << " " << faces[i].indexes[1] << "  " << faces[i].indexes[2] << endl;
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

    }
    face_number = in_index;
}

void PLYObstacle::createAABBs()
{
    aabbs.resize(face_number);
    for (unsigned i = 0; i < face_number; i++){
        aabbs[i] = faces[i].computeAABB();
    }
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
        colision.obstacle_id = id;
        bounced = true;
        if (colision.col_location == Collision::on_edge || colision.col_location == Collision::on_vertex){
            colision.bounced_direction = -step;
        }
        else
        {
            Eigen::Vector3d normal;
            faces[colision.triangle_ind].getNormal(normal);

            Eigen::Vector3d temp_step = step;
            elasticBounceAgainsPlane(ray_origin,normal,colision.t,temp_step);

            colision.bounced_direction = temp_step.normalized();

            //Orientation respect the triangle
            double dot = ((walker.pos_v - faces[colision.triangle_ind].center).normalized()).dot(normal);

            colision.col_location = (dot < -1e-5)?Collision::inside:(dot > 1e-5)?Collision::outside:Collision::unknown;

            //WARNING: Cuidar este patch
            // Implementa Percolacion
            if(this->percolation>0.0){
                bool from_intra = (colision.col_location == Collision::inside);
                // Count the hit (membrane reached, a crossing draw is made). Counts both
                // step and bouncing hits, since checkCollision runs for each. P0.3 validation.
                if(from_intra) count_hits_i_e++; else count_hits_e_i++;

                // Seeded, thread-safe per-walker draw (was C rand()/RAND_MAX). P0.1.
                double _percolation_ (walker.rng.uniform());
                double dynamic_percolation = from_intra?this->prob_cross_i_e:this->prob_cross_e_i;

                if( dynamic_percolation - _percolation_ > EPS_VAL ){
                    count_perc_crossings++;
                    if(from_intra) count_cross_i_e++; else count_cross_e_i++;
                    walker.perm_crossed_flag = true;
                    colision.bounced_direction = step;
                    //colision.type = Collision::null;
                    return false;
                }
            }
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

bool PLYObstacle::isPLYBinary(std::string ply_file)
{
    std::ifstream in(ply_file.c_str(), std::ifstream::in);
    if (!in) {
        std::cout << "Error opening file " << ply_file << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.find("binary") != std::string::npos) {
            in.close();
            return true;
        }
        if (line.find("ascii") != std::string::npos) {
            in.close();
            return false;
        }
        if (line == "end_header") {
            break;
        }
    }
    in.close();
    return false; // Default to ASCII if not specified
}

void PLYObstacle::readPLY_Binary(std::string ply_file)
{
    if (vertices != nullptr)
        delete[] vertices;
    if (faces != nullptr)
        delete[] faces;

    // Open file in binary mode
    std::ifstream in(ply_file.c_str(), std::ios::binary);

    if (!in) {
        std::cout << "Error opening file " << ply_file << std::endl;
        assert(1);
        return;
    }

    // Parse header (still in ASCII)
    std::string line;
    
    // Read first line - should be "ply"
    std::getline(in, line);
    if (line.compare("ply") != 0) {
        std::cout << "Not a valid PLY file: missing 'ply' header" << std::endl;
        in.close();
        assert(1);
        return;
    }
    
    // Parse header to get element counts
    vert_number = 0;
    face_number = 0;
    
    while (std::getline(in, line)) {
        // Check if header section is done
        if (line == "end_header") {
            break;
        }
        
        // Extract element counts
        if (line.find("element vertex") != std::string::npos) {
            sscanf(line.c_str(), "element vertex %u", &vert_number);
        } else if (line.find("element face") != std::string::npos) {
            sscanf(line.c_str(), "element face %u", &face_number);
        }
    }
    
    if (vert_number == 0 || face_number == 0) {
        std::cout << "Invalid PLY file: missing vertex or face counts" << std::endl;
        in.close();
        assert(1);
        return;
    }
    
    // Allocate memory for vertices and faces
    vertices = new Vertex[vert_number];
    faces = new Triangle[face_number];
    
    // Read vertex data
    for (unsigned i = 0; i < vert_number; i++) {
        float x, y, z;
        in.read(reinterpret_cast<char*>(&x), sizeof(float));
        in.read(reinterpret_cast<char*>(&y), sizeof(float));
        in.read(reinterpret_cast<char*>(&z), sizeof(float));
        
        vertices[i].points[0] = static_cast<double>(x) * scale_factor;
        vertices[i].points[1] = static_cast<double>(y) * scale_factor;
        vertices[i].points[2] = static_cast<double>(z) * scale_factor;
    }
    
    // Read face data
    for (unsigned i = 0; i < face_number; i++) {
        uint8_t num_vertices;
        in.read(reinterpret_cast<char*>(&num_vertices), sizeof(uint8_t));
        
        if (num_vertices != 3) {
            std::cout << "Non-triangle face detected. Only triangular faces are supported." << std::endl;
            in.close();
            assert(1);
            return;
        }
        
        uint32_t v1, v2, v3;
        in.read(reinterpret_cast<char*>(&v1), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&v2), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&v3), sizeof(uint32_t));
        
        faces[i].indexes[0] = v1;
        faces[i].indexes[1] = v2;
        faces[i].indexes[2] = v3;
        faces[i].vertices = vertices;
        faces[i].saveNormalAndAuxInfo();
    }
    
    in.close();
}

void PLYObstacle::readPLY_Binary_trianglesSubdivitionDistance(std::string ply_file, std::vector<Eigen::Vector3d>& centers, double max_distance)
{
    if (vertices != nullptr)
        delete[] vertices;
    if (faces != nullptr)
        delete[] faces;

    // Open file in binary mode
    std::ifstream in(ply_file.c_str(), std::ios::binary);

    if (!in) {
        std::cout << "Error opening file " << ply_file << std::endl;
        assert(0);
        return;
    }

    // Parse header (still in ASCII)
    std::string line;
    
    // Read first line - should be "ply"
    std::getline(in, line);
    if (line.compare("ply") != 0) {
        std::cout << "Not a valid PLY file: missing 'ply' header" << std::endl;
        in.close();
        assert(0);
        return;
    }
    
    // Parse header to get element counts
    vert_number = 0;
    face_number = 0;
    
    while (std::getline(in, line)) {
        // Check if header section is done
        if (line == "end_header") {
            break;
        }
        
        // Extract element counts
        if (line.find("element vertex") != std::string::npos) {
            sscanf(line.c_str(), "element vertex %u", &vert_number);
        } else if (line.find("element face") != std::string::npos) {
            sscanf(line.c_str(), "element face %u", &face_number);
        }
    }
    
    if (vert_number == 0 || face_number == 0) {
        std::cout << "Invalid PLY file: missing vertex or face counts" << std::endl;
        in.close();
        assert(0);
        return;
    }
    
    // Allocate memory for vertices and faces
    vertices = new Vertex[vert_number];
    faces = new Triangle[face_number];
    
    // Read vertex data
    for (unsigned i = 0; i < vert_number; i++) {
        float x, y, z;
        in.read(reinterpret_cast<char*>(&x), sizeof(float));
        in.read(reinterpret_cast<char*>(&y), sizeof(float));
        in.read(reinterpret_cast<char*>(&z), sizeof(float));
        
        vertices[i].points[0] = static_cast<double>(x) * scale_factor;
        vertices[i].points[1] = static_cast<double>(y) * scale_factor;
        vertices[i].points[2] = static_cast<double>(z) * scale_factor;
    }
    
    // Read face data with subdivision distance filter
    int in_index = 0;
    for (unsigned i = 0; i < face_number; i++) {
        uint8_t num_vertices;
        in.read(reinterpret_cast<char*>(&num_vertices), sizeof(uint8_t));
        
        if (num_vertices != 3) {
            std::cout << "Non-triangle face detected. Only triangular faces are supported." << std::endl;
            in.close();
            assert(0);
            return;
        }
        
        uint32_t v1, v2, v3;
        in.read(reinterpret_cast<char*>(&v1), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&v2), sizeof(uint32_t));
        in.read(reinterpret_cast<char*>(&v3), sizeof(uint32_t));
        
        faces[in_index].indexes[0] = v1;
        faces[in_index].indexes[1] = v2;
        faces[in_index].indexes[2] = v3;
        faces[in_index].vertices = vertices;
        faces[in_index].saveNormalAndAuxInfo();
        
        if (centers.size() > 0) {
            bool include_face = false;
            for (auto c : centers) {
                double distance = faces[in_index].minDistance(c);
                if (distance < max_distance) {
                    include_face = true;
                    break;
                }
            }
            if (include_face) {
                in_index++;
            }
        }
        else {
            in_index++;
        }
    }
    
    face_number = in_index;
    in.close();
}



