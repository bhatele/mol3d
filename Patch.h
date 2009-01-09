/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file Patch.h
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */


#ifndef __PATCH_H__
#define __PATCH_H__

/** \class Main
 *
 */
class Main : public CBase_Main {
  private:
    int phase;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    void allDone();
    void startUpDone();
};

class ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  public:
    int lengthX;
    int lengthY;
    int lengthZ;
    double* particleLocX;
    double* particleLocY;
    double* particleLocZ;
    int lengthAll;
    int x;
    int y;
    int z;
};

/** \class Patch
 *  Class representing a cell in the grid. 
 *  We consider each cell as a square of LxL units
 */
class Patch : public CBase_Patch {
  private:
    CkVec<Particle> particles;
    CkVec<Particle> incomingParticles;
    int forceCount;		// to count the returns from interactions
    int stepCount;		// to count the number of steps, and decide when to stop
    int updateCount;
    bool updateFlag;
    bool incomingFlag;
    int computesList[NUM_NEIGHBORS][6];

    void migrateToPatch(Particle p, int &px, int &py, int &pz);
    void updateProperties();	// updates properties after receiving forces from computes
    void limitVelocity(Particle &p);
    Particle& wrapAround(Particle &p);
    void print();		// prints all its particles
    CProxySection_Compute mCastSecProxy;

  public:
    Patch();
    Patch(CkMigrateMessage *msg);
    ~Patch();

    void start();
    void createComputes();
    void createSection();
    void receiveParticles(CkVec<Particle> &);
    void receiveForces(ParticleForceMsg *updates);
    void checkNextStep();	// checks whether to continue with next step
#ifdef RUN_LIVEVIZ
    void requestNextFrame(liveVizRequestMsg *m);
#endif
};

#endif
