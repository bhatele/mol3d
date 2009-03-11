/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file Compute.h
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */

#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "defs.h"

struct force {
  BigReal x;
  BigReal y;
  BigReal z;
};

class ParticleForceMsg : public CMessage_ParticleForceMsg {
  public:
    int lengthX;
    int lengthY;
    int lengthZ;
   
    force* forces;
    int lengthUpdates;
};

// Class representing the interaction agents between a couple of cells
class Compute : public CBase_Compute {
  private:
    int cellCount;  // to count the number of interact() calls
    ParticleDataMsg *bufferedMsg;
    CkVec<int> *pairList;

  public:
    Compute();
    Compute(CkMigrateMessage *msg);

    void interact(ParticleDataMsg *msg);
    //void calcPairForcesPL(ParticleDataMsg* first, ParticleDataMsg* second);
    //void calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second);
    //void calcInternalForces(ParticleDataMsg* first);
           
};

#endif
