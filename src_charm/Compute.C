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

#include "defs.h"
#ifdef RUN_LIVEVIZ
  #include "liveViz.h"
#endif
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
  cellCount = 0;
  numLists = -1;
  bmsgLenAll = -1;
  usesAtSync = CmiTrue;
}

Compute::Compute(CkMigrateMessage *msg) { }
  
// Entry method to receive vector of particles
void Compute::interact(ParticleDataMsg *msg){
  int i;

  // self interaction check
  if (thisIndex.x1==thisIndex.x2 && thisIndex.y1==thisIndex.y2 && thisIndex.z1==thisIndex.z2) {
    if (msg->doAtSync)
      AtSync();
    bmsgLenAll = -1;
    CkGetSectionInfo(cookie1,msg);
    calcInternalForces(msg, &cookie1);
  } else {
    if (cellCount == 0) {
      bufferedMsg = msg;
      bmsgLenAll = bufferedMsg->lengthAll;
      cellCount++;
    } else if (cellCount == 1) {
      // if both particle sets are received, compute interaction
      cellCount = 0;
      if (msg->doAtSync)
	AtSync();
      bmsgLenAll = -1;
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
	if (bufferedMsg->lengthAll <= msg->lengthAll)
          calcPairForces(bufferedMsg, msg);
	else
	  calcPairForces(msg, bufferedMsg);
      }
      bufferedMsg = NULL;
    }
  }
}
