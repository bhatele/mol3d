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
  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    void allDone();
    void computeCreationDone();
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
    void checkNextStep();	// checks whether to continue with next step
    void print();		// prints all its particles

  public:
    Patch();
    Patch(CkMigrateMessage *msg);
    ~Patch();

    void start();
    void createComputes();
    void receiveParticles(CkVec<Particle> &);
    void receiveForces(CkVec<Particle> &);
#ifdef RUN_LIVEVIZ
    void requestNextFrame(liveVizRequestMsg *m);
#endif
};

#endif
