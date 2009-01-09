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

#include "common.h"
#ifdef RUN_LIVEVIZ
  #include "liveViz.h"
#endif
#include "Patch.decl.h"
#include "Compute.h"
#include "Patch.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;

extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;	// Number of Chare Rows
extern /* readonly */ int patchArrayDimY;	// Number of Chare Columns
extern /* readonly */ int patchSize;
extern /* readonly */ double radius;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ double stepTime; 

extern double A;			// Force Calculation parameter 1
extern double B;			// Force Calculation parameter 2

// Compute - Default constructor
Compute::Compute() {
  cellCount = 0;
}

Compute::Compute(CkMigrateMessage *msg) { }
  
// Entry method to receive vector of particles
void Compute::interact(ParticleDataMsg *msg){
  int i;

  // self interaction check
  if (thisIndex.x1==thisIndex.x2 && thisIndex.y1==thisIndex.y2 && thisIndex.z1==thisIndex.z2) {
    //for (int i = 0; i < msg->lengthAll; i++){
    //  CkPrintf("post loc = %f ",msg->particleLocX[i]);
    //}
    calcForces(msg, msg, true);
  } else {
    if (cellCount == 0) {
      bufferedMsg = msg;
      cellCount++;
    } else if (cellCount == 1) {
      // if both particle sets are received, compute interaction
      cellCount = 0;
      calcForces(bufferedMsg, msg, false);
    }
  }
}

// Local function to compute all the interactions between pairs of particles in two sets
void Compute::calcForces(ParticleDataMsg* first, ParticleDataMsg* second, bool bSame) {
  int i, j;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  double rx, ry, rz, r, fx, fy, fz, f;
  ParticleForceMsg *firstmsg = new (firstLen, firstLen, firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg;
  if (bSame)
    secondmsg = firstmsg;
  else
    secondmsg = new (secondLen, secondLen, secondLen) ParticleForceMsg;
 // CkPrintf("firstlen = %d, secondlen = %d \n", firstLen, secondLen);
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  for(i = 0; i < firstLen; i++){
    for(j = 0; j < secondLen; j++) {
      //CkPrintf("(%lf)(%lf) ",first->particleLocX[i],second->particleLocX[j]); 
      // computing base values
      rx = first->particleLocX[i] - second->particleLocX[j];
      ry = first->particleLocY[i] - second->particleLocY[j];
      rz = first->particleLocZ[i] - second->particleLocZ[j];
      r = sqrt(rx*rx + ry*ry + rz*rz);
      // We include 0.000001 to ensure that r doesn't tend to zero in the force calculation
      //if(r < 0.000001 || r >= patchSize) CkPrintf(" radii was OB ");
      if(r >= 0.000001 && r < patchSize){
        f = A / pow(r,12) - B / pow(r,6);
	//CkPrintf("force = %f ", f);
        fx = f * rx / r;
        fy = f * ry / r;
        fz = f * rz / r;
        if (i == 0){
  	  secondmsg->forcesX[j] = -fx;
  	  secondmsg->forcesY[j] = -fy;
	  secondmsg->forcesZ[j] = -fz;
        }
        else{
	  secondmsg->forcesX[j] -= fx;
	  secondmsg->forcesY[j] -= fy;
	  secondmsg->forcesZ[j] -= fz;
        }

        if (j == 0){
	  firstmsg->forcesX[i] = fx;
	  firstmsg->forcesY[i] = fy;
	  firstmsg->forcesZ[i] = fz;
        }
        else{
	  firstmsg->forcesX[i] += fx;
	  firstmsg->forcesY[i] += fy;
	  firstmsg->forcesZ[i] += fz;
        }
      }
      else {
        if (i == 0){
          secondmsg->forcesX[j] = 0;
          secondmsg->forcesY[j] = 0;
          secondmsg->forcesZ[j] = 0;
        }
        if (j == 0){
          firstmsg->forcesX[i] = 0;
          firstmsg->forcesY[i] = 0;
          firstmsg->forcesZ[i] = 0;
        }
      }
    }
  }
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  if (!bSame)
    patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
}

