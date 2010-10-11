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

class sqrtPars {
  public:
    BigReal a; //constant
    BigReal b;
    BigReal c;
    BigReal d; //x^3

    void pup(PUP::er &p) {
      p | a; p | b;
      p | c;  p | d;
    }
};


class sqrtTable : public CMessage_sqrtTable {
  public:
    sqrtPars* pars;
    int length;
    BigReal delta;
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
    int numLists;
    int bmsgLenAll;
    ParticleDataMsg *bufferedMsg;
    CkVec<int> *pairList;
    CkSectionInfo cookie1;
    CkSectionInfo cookie2;

  public:
    Compute();
    Compute(CkMigrateMessage *msg);

    void interact(ParticleDataMsg *msg);
    //void calcPairForcesPL(ParticleDataMsg* first, ParticleDataMsg* second);
    //void calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second);
    //void calcInternalForces(ParticleDataMsg* first);

    void pup(PUP::er &p) {
      CBase_Compute::pup(p);
      p | cookie1;
      p | cookie2;
#ifdef USE_SECTION_MULTICAST
      if (p.isUnpacking() && CkInRestarting()) {
        cookie1.get_redNo() = 0;
        if (!(thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2))
        cookie2.get_redNo() = 0;
      }
#endif
      p | cellCount;
      p | numLists;
      p | bmsgLenAll;
      int hasMsg = (bmsgLenAll >= 0); // only pup if msg will be used
      p | hasMsg;
      if (hasMsg){
	//CkPrintf("HERE?\n");
	if (p.isUnpacking())
	  bufferedMsg = new (bmsgLenAll) ParticleDataMsg;
	p | *bufferedMsg;
      }
      else
        bufferedMsg = NULL;

      int hasList = (numLists >= 0  && pairList != NULL);
      p | hasList;
      if (hasList){
	//CkPrintf("NUMLISTS = %d\n", numLists);
	if (p.isUnpacking())
	  pairList = new CkVec<int>[numLists];
	PUParray(p, pairList, numLists);
      }
      else{
//	CkPrintf("DID NOT NEED TO PUP PAIRLIST\n");
	pairList = NULL;
      }
      //CkPrintf("done pupping a compute \n");
    }
    void ResumeFromSync();           
};

#endif
