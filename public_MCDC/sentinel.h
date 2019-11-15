//!  Collision Final class ====================================================================/
/*!
    \details    Auxiliar class to check error and misbehaviours during the particle dynamics.
    \author     Jonathan Rafael
    \date       Junes 2017

*=============================================================================================*/
#ifndef SENTINEL_H
#define SENTINEL_H
#include "walker.h"
#include "parameters.h"

/*! \class Sentinels
 *  \brief Class used to check the possible numerical errors or un-handed cases inside the
 *  dynamic simulation.
 */

namespace sentinels{;

class Sentinel
{
public:

    unsigned stuck_count;
    unsigned illegal_count;

    enum ErrorCases{none,stuck,crossed,rejected,rejected_initial_pos};

    unsigned bouncings;
    unsigned obstacle_id;
    unsigned rejected_count;
    bool rejected_step;
    bool deport_illegals;
    bool discard_stucks ;

    ErrorCases error;

    Sentinel();

    void clear();

    void setBouncingError(unsigned bouncings);

    void setCrossingError(unsigned);

    void setRejectedError();

    bool checkErrors(Walker &w, const Parameters &params, bool noPLY, unsigned &bouncing_count);

    void deportationProcess(Walker &walker, unsigned &w, unsigned& t , bool& back_tracking, Parameters &params, int id);

};

}
#endif // SENTINEL_H
