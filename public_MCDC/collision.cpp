#include "collision.h"
#include <cmath>
const double EPS_VAL = 1e-15;

Collision::~Collision()
{

}


bool Collision::doIHaveMorePiorityThan(Collision &coll)
{
    // if my type is the same
    if (coll.type == type){
        return t < coll.t;
    }
    else{
        return type < coll.type;
    }
}

// Computes the location on the collision, or near position.
// Referent to be, on the vertex, on the edge, or inside.
void Collision::computeCollisionLocation()
{

    col_location = unknown;

    bool on_edge_flag = ( std::abs(u) < EPS_VAL );  // u=0;
    on_edge_flag     |= ( std::abs(v) < EPS_VAL);  //  v=0;
    on_edge_flag     |= ( std::abs(1.0-u-v) < EPS_VAL);  // u+v = 1 => w=0;

    if(on_edge_flag){
        col_location = on_edge;
        return;
    }
}











