/** \file Main.C
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */

#include <time.h>
#include "ckmulticast.h"

#include "defs.h"
#include "mol3d.decl.h"
#include "Reader.h"

#include "Main.h"
#include "Patch.h"
#include "Compute.h"

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Patch patchArray;
/* readonly */ CProxy_Compute computeArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ int numParts;
/* readonly */ int patchArrayDimX;	// Number of Chares in X
/* readonly */ int patchArrayDimY;	// Number of Chares in Y
/* readonly */ int patchArrayDimZ;	// Number of Chares in Z
/* readonly */ int patchSizeX;
/* readonly */ int patchSizeY;
/* readonly */ int patchSizeZ;

/* readonly */ BigReal stepTime; 
/* readonly */ bool twoAwayX;
/* readonly */ bool twoAwayY;
/* readonly */ bool twoAwayZ;
/* readonly */ int numNbrs;
/* readonly */ int nbrsX;
/* readonly */ int nbrsY;
/* readonly */ int nbrsZ;

// double A = 1.60694452*pow(10.0, -134);
// double B = 1.031093844*pow(10.0, -77);

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {
  LBTurnInstrumentOff();
  stepTime = CmiWallTimer();
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS START UP ...\n");
  numParts = 0;

  /*patchArrayDimX = PATCHARRAY_DIM_X;
  patchArrayDimY = PATCHARRAY_DIM_Y;
  patchArrayDimZ = PATCHARRAY_DIM_Z;
  patchSizeX = PATCH_SIZE_X;
  patchSizeY = PATCH_SIZE_Y;
  patchSizeZ = PATCH_SIZE_Z;*/
  numNbrs = NUM_NEIGHBORS;
  nbrsX = NBRS_X;
  nbrsY = NBRS_Y;
  nbrsZ = NBRS_Z;

  mainProxy = thisProxy;
  phase = 0;
  int bFactor = 2;

  if (msg->argc < 2) {
    CkAbort("Atleast specify the config file to run this programn\n\n");
  }
  if (msg->argc > 2) bFactor = atoi(msg->argv[2]);

#ifdef USE_SECTION_MULTICAST
  // initializing the CkMulticastMgr
  mCastGrpID = CProxy_CkMulticastMgr::ckNew(bFactor);
#endif
}

void Main::startUpDone() {
  switch(phase) {
    case 0:
      CkPrintf("\nNUMBER OF PATCHES: %d X %d X %d .... CREATED\n", patchArrayDimX, patchArrayDimY, patchArrayDimZ);
      phase++;

      // initializing the 6D compute array
      for (int x=0; x<patchArrayDimX; x++)
	for (int y=0; y<patchArrayDimY; y++)
	  for (int z=0; z<patchArrayDimZ; z++)
	    patchArray(x, y, z).createComputes();
      break;
 
    case 1:
      computeArray.doneInserting();
      CkPrintf("NUMBER OF COMPUTES: %d .... CREATED\n", (numNbrs/2+1) * patchArrayDimX * patchArrayDimY * patchArrayDimZ);
#ifdef USE_SECTION_MULTICAST
      phase++;
      patchArray.createSection();
      break;

    case 2:
      CkPrintf("MULTICAST SECTIONS .... CREATED\n");
#endif

      CkPrintf("STARTING SIMULATION .... \n\n");
      patchArray.start();
      break;
  }
}

/* Constructor for chare object migration */
Main::Main(CkMigrateMessage* msg): CBase_Main(msg) { }

void Main::lbBarrier(){
  // CkPrintf("got to lbBarrier at %f\n", CmiWallTimer());
  patchArray.resume();
}

void Main::ftBarrier(){
  // CkPrintf("got to ftBarrier at %f\n", CmiWallTimer());
  CkCallback cb(CkIndex_Patch::ftresume(), patchArray);
  CkStartMemCheckpoint(cb);
}

void Main::allDone() {
  CkPrintf("SIMULATION COMPLETE. \n\n");  CkExit();
}

#include "mol3d.def.h"
