/** \file Compute.h
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */

#include "defs.h"
#include "mol3d.decl.h"
#include "Patch.h"
#include "Compute.h"
#include "ckmulticast.h"
#include "nonbonded.h"
#include "ConfigList.h"
//#include "sqrtTable.h"


extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ sqrtTable *rootTable;

extern /* readonly */ bool usePairLists;
extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;	// Number of Chare Rows
extern /* readonly */ int patchArrayDimY;	// Number of Chare Columns
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int patchSizeX;
extern /* readonly */ int patchSizeY;
extern /* readonly */ int patchSizeZ;
extern /* readonly */ int ptpCutOff;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ BigReal stepTime; 

extern /* readonly */ double A;			// Force Calculation parameter 1
extern /* readonly */ double B;			// Force Calculation parameter 2

// Compute - Default constructor
Compute::Compute() {
  LBTurnInstrumentOff();
  cellCount = 0;
  numLists = -1;
  bmsgLenAll = -1;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg): CBase_Compute(msg)  { 
  usesAtSync = CmiTrue;
  delete msg;
}
  
// Entry method to receive vector of particles
void Compute::interact(ParticleDataMsg *msg){
  int i;

  // self interaction check
  // FIX THIS: :self compute could be a wrap
  if (thisIndex.x1 ==thisIndex.x2 && thisIndex.y1 ==thisIndex.y2 && thisIndex.z1 ==thisIndex.z2) {
    bool doatSync = false;
    bmsgLenAll = -1;
    if (msg->doAtSync){
     // LBTurnInstrumentOff();
      //AtSync();
      doatSync = true;
      //LBTurnInstrumentOff();
    }
    if (msg->lbOn)
      LBTurnInstrumentOn();
    CkGetSectionInfo(cookie1,msg);
    calcInternalForces(msg, &cookie1);
    if(doatSync)
      AtSync();
    //if (msg->lbOn)
      //LBTurnInstrumentOn();
  } else {
    if (cellCount == 0) {
      bufferedMsg = msg;
      bmsgLenAll = bufferedMsg->lengthAll;
      cellCount++;
    } else if (cellCount == 1) {
      // if both particle sets are received, compute interaction
      cellCount = 0;
      bool doatSync = false;
      bmsgLenAll = -1;
      if (msg->doAtSync){
	//LBTurnInstrumentOff();
	//AtSync();
	doatSync = true;
      }
      if (msg->lbOn)
	LBTurnInstrumentOn();
      if (usePairLists){
//	if (bufferedMsg->lengthAll == msg->lengthAll){
	  if (bufferedMsg->x*patchArrayDimY*patchArrayDimZ + bufferedMsg->y*patchArrayDimZ + bufferedMsg->z < msg->x*patchArrayDimY*patchArrayDimZ + msg->y*patchArrayDimZ + msg->z){ 
	    if (bufferedMsg->lengthAll <= msg->lengthAll)
	      pairList = calcPairForcesPL(bufferedMsg, msg, pairList, &numLists, &cookie1, &cookie2);
	    else
	      pairList = calcPairForcesPL(msg, bufferedMsg, pairList, &numLists, &cookie2, &cookie1);
	  }
	  else{
	    if (bufferedMsg->lengthAll <= msg->lengthAll)
	      pairList = calcPairForcesPL(bufferedMsg, msg, pairList, &numLists, &cookie2, &cookie1);
	    else
	      pairList = calcPairForcesPL(msg, bufferedMsg, pairList, &numLists, &cookie1, &cookie2);
	  }
      }
      else {
	if (bufferedMsg->x*patchArrayDimY*patchArrayDimZ + bufferedMsg->y*patchArrayDimZ + bufferedMsg->z < msg->x*patchArrayDimY*patchArrayDimZ + msg->y*patchArrayDimZ + msg->z){ 
	  if (bufferedMsg->lengthAll <= msg->lengthAll)
	    calcPairForces(bufferedMsg, msg, &cookie1, &cookie2);
	  else
	    calcPairForces(msg, bufferedMsg, &cookie2, &cookie1);
	}
	else{
	  if (bufferedMsg->lengthAll <= msg->lengthAll)
	    calcPairForces(bufferedMsg, msg, &cookie2, &cookie1);
	  else
	    calcPairForces(msg, bufferedMsg, &cookie1, &cookie2);
	}
      }
      //if (msg->lbOn)
	//LBTurnInstrumentOn();
      bufferedMsg = NULL;
      if(doatSync)
	AtSync();
    }
  }
}


void Compute::ResumeFromSync(){
  LBTurnInstrumentOff();
  //CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  /*if (thisIndex.x1==thisIndex.x2 && thisIndex.y1==thisIndex.y2 && thisIndex.z1==thisIndex.z2) {
    patchArray(thisIndex.x1, thisIndex.y1, thisIndex.z1).resume();
    //mCastGrp->rebuild(cookie1);
  }
  else{
    patchArray(thisIndex.x1, thisIndex.y1, thisIndex.z1).resume();
    patchArray(thisIndex.x2, thisIndex.y2, thisIndex.z2).resume();
    //mCastGrp->rebuild(cookie1);
    //mCastGrp->rebuild(cookie2);
  }*/
}
