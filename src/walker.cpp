/*!
 * Based class Walker.
*/

#include "walker.h"
#include <stdlib.h>     /* srand, rand */
#include <random>       /*random_device for MAC */
#include "cylinder.h"

Walker::Walker()
{
    pos_r = pos_v.setZero(3,1);
    index = 0;
    status = free;
    steps_count = 0;
    initial_location = location = unknown;
    intra_extra_consensus = intra_coll_count = extra_coll_count = rejection_count = steps_count = 0;
    steps_per_second = 0;
    in_ply_index = -1;
    in_obj_index =-1;
}

Walker::Walker(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);


    double t = udist(gen);
    pos_r[0] = last_pos_r[0] = last_pos_v[0] = ini_pos[0] = pos_v[0] = (t+1)*xmin + t*xmax;
    t = udist(gen);
    pos_r[1] = last_pos_r[1] = last_pos_v[1] = ini_pos[1] = pos_v[1] = (t+1)*ymin + t*ymax;
    t = udist(gen);
    pos_r[2] = last_pos_r[2] = last_pos_v[2] = ini_pos[2] = pos_v[2] = (t+1)*zmin + t*zmax;

    steps_count     = 0;
    rejection_count = 0;
    status = free;
    index = 0;
    initial_location  = location = unknown;
    intra_extra_consensus = intra_coll_count = extra_coll_count =0;
    steps_per_second = 0;
    in_ply_index = -1;
    in_obj_index = -1;
}

void Walker::getRealPosition(double &x_, double &y_, double &z_) const
{
    x_ = pos_r[0];
    y_ = pos_r[1];
    z_ = pos_r[2];
}
void Walker::getVoxelPosition(double &x_, double &y_, double &z_) const
{
    x_ = pos_v[0];
    y_ = pos_v[1];
    z_ = pos_v[2];
}

void Walker::getRealPosition(Eigen::Vector3d &_pos_) const
{
    _pos_ = pos_r;
}

void Walker::getVoxelPosition(Eigen::Vector3d &_pos_) const
{
    _pos_ = pos_v;
}

void Walker::getInitialPosition(double &x_, double &y_, double &z_) const
{
    x_ = ini_pos[0];
    y_ = ini_pos[1];
    z_ = ini_pos[2];
}

void Walker::getInitialPosition(Eigen::Vector3d &_ini_pos_) const
{
    _ini_pos_ = ini_pos;
}

void Walker::getNextDirection(Eigen::Vector3d &_direction_) const
{
    _direction_ = next_direction;

}

unsigned int Walker::getIndex() const
{
    return index;
}

void Walker::setRealPosition(const Eigen::Vector3d &_pos_)
{
    last_pos_r = pos_r;
    pos_r = _pos_;
}

void Walker::setVoxelPosition(const Eigen::Vector3d &_pos_)
{
    last_pos_v = pos_v;
    pos_v = _pos_;
}
void Walker::setVoxelPosition(const double &x_, const double &y_, const double &z_)
{
    last_pos_v = pos_v;

    pos_v[0] = x_;
    pos_v[1] = y_;
    pos_v[2] = z_;
}

void Walker::setRealPosition(const double &x_, const double &y_, const double &z_)
{
    last_pos_r = pos_r;

    pos_r[0] = x_;
    pos_r[1] = y_;
    pos_r[2] = z_;
}

void Walker::setInitialPosition(const double &x_,const double &y_,const double &z_)
{
    pos_r[0]  = x_;
    pos_r[1]  = y_;
    pos_r[2]  = z_;

    last_pos_r = last_pos_v = pos_v = ini_pos = pos_r;

    steps_count = 0;
}


void Walker::setInitialPosition(const Eigen::Vector3d &temp)
{
    last_pos_r = last_pos_v = pos_v = ini_pos = pos_r = temp;

    steps_count = 0;
}

void Walker::setNextDirection(Eigen::Vector3d &next_dir)
{
    next_direction = next_dir;
}

void Walker::setIndex(unsigned int& _index)
{
    index = _index;
}

void Walker::setRealPosLog(const Eigen::Vector3d &pos, unsigned t)
{
    this->pos_r_log(0,t)=pos(0);
    this->pos_r_log(1,t)=pos(1);
    this->pos_r_log(2,t)=pos(2);
}

void Walker::setRealPosLog(double x, double y, double z, unsigned t)
{
    this->pos_r_log(0,t)=x;
    this->pos_r_log(1,t)=y;
    this->pos_r_log(2,t)=z;
}

void Walker::setVoxPosLog(const Eigen::Vector3d &pos, unsigned t)
{
    this->pos_v_log(0,t)=pos(0);
    this->pos_v_log(1,t)=pos(1);
    this->pos_v_log(2,t)=pos(2);
}

void Walker::setVoxPosLog(double x, double y, double z, unsigned t)
{
    this->pos_v_log(0,t)=x;
    this->pos_v_log(1,t)=y;
    this->pos_v_log(2,t)=z;
}

void Walker::setNumberOfSteps(unsigned T)
{
    pos_r_log = Eigen::Matrix3Xd::Zero(3,T+1);
    pos_v_log = Eigen::Matrix3Xd::Zero(3,T+1);
}


void Walker::setRandomInitialPosition(const Eigen::Vector3d &_min, const Eigen::Vector3d &_max)
{
    steps_count = 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> udist(0,1);

    double t = udist(gen);
    pos_r[0]  = (1-t)*_min(0)+ t*_max(0);
    t = udist(gen);
    pos_r[1]  = (1-t)*_min(1)+ t*_max(1);
    t = udist(gen);
    pos_r[2]  =  (1-t)*_min(2)+ t*_max(2);

    last_pos_r = last_pos_v = pos_v = ini_pos = pos_r;
    last_pos_r = last_pos_v = pos_r;
}





