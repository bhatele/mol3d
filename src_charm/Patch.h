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

struct loc{
  BigReal x;
  BigReal y;
  BigReal z;
};

class ParticleDataMsg : public CkMcastBaseMsg, public CMessage_ParticleDataMsg {
  public:
    //BigReal* particleLocX;
    //BigReal* particleLocY;
    //BigReal* particleLocZ;
    loc* coords;
    BigReal* charge;
    int lengthAll;
    int x;
    int y;
    int z;
    bool updateList;
    bool deleteList;
};

class FileDataMsg : public CMessage_FileDataMsg {
  public:
    BigReal* charge;
    BigReal* mass;
    loc* coords; //encoded as x1 y1 z1 x2 y2 z2...
    int length;
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
    int myNumParts;
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
    Patch(FileDataMsg *fdmsg);
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
