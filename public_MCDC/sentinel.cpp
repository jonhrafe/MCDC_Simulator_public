#include "sentinel.h"
#include "constants.h"
#include <iostream>
using namespace sentinels;

using namespace std;

Sentinel::Sentinel()
{
    error = none;

    bouncings = 0;

    obstacle_id = 0;

    stuck_count = 0;

    rejected_count = 0;

    illegal_count = 0;

    deport_illegals = true; //Trump mode on.
    discard_stucks = true;

    rejected_step = false;
}

void Sentinel::clear(){
    error = none;

    bouncings = 0;

    obstacle_id = 0;
}

void Sentinel::setBouncingError(unsigned bouncings){
    error = stuck;

    this->bouncings = bouncings;
}

void Sentinel::setCrossingError(unsigned id){
    error = crossed;
    this->obstacle_id = id;
}
void Sentinel::setRejectedError(){
    error = rejected;
}


bool Sentinel::checkErrors(Walker &walker, const Parameters &params, bool noPLY, unsigned& bouncing_count)
{

    //If it was rejected to many times we have to take out the walker.
    if(walker.rejection_count>max_rejections){
        setBouncingError(walker.steps_count);
        stuck_count++;
        throw(this->error);
    }

    //Handle "stucked" walkers
    if( (bouncing_count > max_number_bouncings) && discard_stucks){
        setBouncingError(walker.steps_count);
        stuck_count++;
        throw(this->error);
    }

    if( (walker.location != Walker::unknown) && (params.obstacle_permeability <= 0.0) && deport_illegals == true ){
        if(walker.initial_location != walker.location){
            setCrossingError(uint(walker.in_obj_index));
            illegal_count++;
            throw(this->error);
        }
    }

    if(this->rejected_step == true){
        rejected_step = false;
        setRejectedError();
        rejected_count++;
        throw(this->error);
    }
    return false;
}

void Sentinel::deportationProcess(Walker &walker, unsigned& w, unsigned &t, bool &back_tracking, Parameters &params, int id)
{
    if (this->error == Sentinel::ErrorCases::stuck){
        //If the particle got stuck in a corner or bad defined area.
        if(params.verbatim)
            cout << endl <<  SH_FG_GRAY <<  "[INFO]   " << SH_DEFAULT << " Sim: " << id << " " <<
                    "Walker "<< w << " labeled as 'stuck' after " << this->bouncings <<
                    " bouncings.\nBacktraking...\nDone" << endl;
        w--;
        back_tracking = true;
    }

    if(this->error == Sentinel::ErrorCases::crossed){
        //If the particle crosses and object because numerical problems
        if(params.verbatim)
            cout << endl <<  SH_FG_GRAY <<  "[INFO]   " << SH_DEFAULT << " Sim: " << id << " " <<
                    "Walker "<< w << " labeled as 'illegal' after crossing obstacle id: " << this->obstacle_id <<
                    "\nBacktraking...\nDone" << endl;
        w--;
        back_tracking = true;
    }

    if(this->error == Sentinel::ErrorCases::rejected){
        t--;

        if(t>1){
            walker.pos_r = Eigen::Vector3d(walker.pos_r_log(0,t),walker.pos_r_log(1,t),walker.pos_r_log(2,t));
            walker.last_pos_r = Eigen::Vector3d(walker.pos_r_log(0,t-1),walker.pos_r_log(1,t-1),walker.pos_r_log(2,t-1));

            walker.pos_v = Eigen::Vector3d(walker.pos_v_log(0,t),walker.pos_v_log(1,t),walker.pos_v_log(2,t));
            walker.last_pos_v = Eigen::Vector3d(walker.pos_v_log(0,t-1),walker.pos_v_log(1,t-1),walker.pos_v_log(2,t-1));
        }
        walker.rejection_count++;
    }
}
