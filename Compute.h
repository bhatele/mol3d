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

#include "common.h"

// Class representing the interaction agents between a couple of cells
class Compute : public CBase_Compute {
  private:
    int cellCount;  // to count the number of interact() calls
    CkVec<Particle> bufferedParticles;
    int bufferedX;
    int bufferedY;
    int bufferedZ;

    void interact(CkVec<Particle> &first, CkVec<Particle> &second);
    void interact(Particle &first, Particle &second);

  public:
    Compute();
    Compute(CkMigrateMessage *msg);

    void interact(CkVec<Particle> particles, int x, int y, int z);
    void calcForces(CkVec<Particle> &first, CkVec<Particle> &second);
};

#endif
