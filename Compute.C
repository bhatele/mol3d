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
#include "Main.decl.h"
#include "Compute.h"
#include "Patch.h"
#include "ConfigList.h"

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;

extern /* readonly */ bool usePairLists;
extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;	// Number of Chare Rows
extern /* readonly */ int patchArrayDimY;	// Number of Chare Columns
extern /* readonly */ int patchArrayDimZ;
extern /* readonly */ int patchSize;
extern /* readonly */ int ptpCutOff;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ BigReal stepTime; 

extern /* readonly */ double A;			// Force Calculation parameter 1
extern /* readonly */ double B;			// Force Calculation parameter 2

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
    calcInternalForces(msg);
  } else {
    if (cellCount == 0) {
      bufferedMsg = msg;
      cellCount++;
    } else if (cellCount == 1) {
      // if both particle sets are received, compute interaction
      cellCount = 0;
      if (usePairLists){
	if (bufferedMsg->lengthAll <= msg->lengthAll)
	  calcPairForcesPL(bufferedMsg, msg);
	else
	  calcPairForcesPL(msg, bufferedMsg);
      }
      else {
	if (bufferedMsg->lengthAll <= msg->lengthAll)
          calcPairForces(bufferedMsg, msg);
	else
	  calcPairForces(msg, bufferedMsg);
      }
    }
  }
}

//calculate pair forces using pairlists
void Compute::calcPairForcesPL(ParticleDataMsg* first, ParticleDataMsg* second) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  ParticleForceMsg *firstmsg = new (firstLen, firstLen, firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen, secondLen, secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = patchSize * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->particleLocX[i] += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = patchSize * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->particleLocY[i] += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = patchSize * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->particleLocZ[i] += diff;
  } 
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10, -20);

  //check if pairlist needs to be updated
  if (first->updateList){
    pairList = new CkVec<int>[firstLen];
    for (i = 0; i < firstLen; i++){
      for (j = 0; j < secondLen; j++){
        rx = first->particleLocX[i] - second->particleLocX[j];
        ry = first->particleLocY[i] - second->particleLocY[j];
        rz = first->particleLocZ[i] - second->particleLocZ[j];
	rsqd = rx*rx + ry*ry + rz*rz;
	if (rsqd >= 0.001 && rsqd < ptpCutOffSqd)
	  pairList[i].push_back(j);
      }
    }
  }

  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
 
  memset(firstmsg->forcesX, 0, firstLen * sizeof(BigReal));
  memset(firstmsg->forcesY, 0, firstLen * sizeof(BigReal));
  memset(firstmsg->forcesZ, 0, firstLen * sizeof(BigReal));
  memset(secondmsg->forcesX, 0, secondLen * sizeof(BigReal));
  memset(secondmsg->forcesY, 0, secondLen * sizeof(BigReal));
  memset(secondmsg->forcesZ, 0, secondLen * sizeof(BigReal));
  for(i = 0; i < firstLen; i++){
    eField = first->charge[i];
    for(j = 0; j < pairList[i].length(); j++) {
      jpart = pairList[i][j];
      rx = first->particleLocX[i] - second->particleLocX[jpart];
      ry = first->particleLocY[i] - second->particleLocY[jpart];
      rz = first->particleLocZ[i] - second->particleLocZ[jpart];
      rsqd = rx*rx + ry*ry + rz*rz;
      if (rsqd >= 0.001){
	rsqd = rsqd * powTwenty;
	r = sqrt(rsqd);
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
	//CkPrintf("rsqd = %E, rsix = %E, rtwelve = %E\n", rsqd, rSix, rTwelve);
        f = (BigReal)(A / rTwelve - B / rSix);
	//CkPrintf("%E \n\n", f);
	f -= eField * constants * second->charge[jpart] / rsqd;
	//CkPrintf("%E \n", f);
	fr = f /r;
	fx = rx * fr;
	fy = ry * fr;
	fz = rz * fr;
	secondmsg->forcesX[jpart] -= fx;
	secondmsg->forcesY[jpart] -= fy;
	secondmsg->forcesZ[jpart] -= fz;
	firstmsg->forcesX[i] += fx;
	firstmsg->forcesY[i] += fy;
	firstmsg->forcesZ[i] += fz;
      }
    }
  }
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
  if (first->deleteList){
    for (i = 0; i < firstLen; i++){
      pairList[i].removeAll();
    }
    delete [] pairList;
  }
  delete first;
  delete second;
}

//calculate pair forces without using pairlists
void Compute::calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  ParticleForceMsg *firstmsg = new (firstLen, firstLen, firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen, secondLen, secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = patchSize * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->particleLocX[i] += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = patchSize * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->particleLocY[i] += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = patchSize * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->particleLocZ[i] += diff;
  } 
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10, -20);
  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
 
  memset(firstmsg->forcesX, 0, firstLen * sizeof(BigReal));
  memset(firstmsg->forcesY, 0, firstLen * sizeof(BigReal));
  memset(firstmsg->forcesZ, 0, firstLen * sizeof(BigReal));
  memset(secondmsg->forcesX, 0, secondLen * sizeof(BigReal));
  memset(secondmsg->forcesY, 0, secondLen * sizeof(BigReal));
  memset(secondmsg->forcesZ, 0, secondLen * sizeof(BigReal));
  for(i = 0; i < firstLen; i++){
    eField = first->charge[i];
    for(jpart = 0; jpart < secondLen; jpart++){
      rx = first->particleLocX[i] - second->particleLocX[jpart];
      ry = first->particleLocY[i] - second->particleLocY[jpart];
      rz = first->particleLocZ[i] - second->particleLocZ[jpart];
      rsqd = rx*rx + ry*ry + rz*rz;
      if (rsqd >= 0.001 && rsqd < ptpCutOffSqd){
	rsqd = rsqd * powTwenty;
	r = sqrt(rsqd);
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * second->charge[jpart] / rsqd;
	fr = f /r;
	fx = rx * fr;
	fy = ry * fr;
	fz = rz * fr;
	secondmsg->forcesX[jpart] -= fx;
	secondmsg->forcesY[jpart] -= fy;
	secondmsg->forcesZ[jpart] -= fz;
	firstmsg->forcesX[i] += fx;
	firstmsg->forcesY[i] += fy;
	firstmsg->forcesZ[i] += fz;
      }
    }
  }
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
  delete first;
  delete second;
}

// Local function to compute all the interactions between pairs of particles in two sets
void Compute::calcInternalForces(ParticleDataMsg* first) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  BigReal powTwenty, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  ParticleForceMsg *firstmsg = new (firstLen, firstLen, firstLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
    
  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
  
  memset(firstmsg->forcesX, 0, firstLen * sizeof(BigReal));
  memset(firstmsg->forcesY, 0, firstLen * sizeof(BigReal));
  memset(firstmsg->forcesZ, 0, firstLen * sizeof(BigReal));
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10, -20);
  for(i = 0; i < firstLen; i++){
    eField = first->charge[i];
    for(j = i+1; j < firstLen; j++) {
      // computing base values
      rx = first->particleLocX[i] - first->particleLocX[j];
      ry = first->particleLocY[i] - first->particleLocY[j];
      rz = first->particleLocZ[i] - first->particleLocZ[j];
      rsqd = rx*rx + ry*ry + rz*rz;
      //check if r >= .000001 to make sure force calc doesnt tend to 0
      if(rsqd >= 0.001 && rsqd < ptpCutOffSqd){
	rsqd = rsqd * powTwenty;
	r = sqrt(rsqd);
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * first->charge[j] / (rsqd);
        
	fr = f /r;
	fx = rx * fr;
        fy = ry * fr;
        fz = rz * fr;
	firstmsg->forcesX[j] -= fx;
	firstmsg->forcesY[j] -= fy;
	firstmsg->forcesZ[j] -= fz;
	firstmsg->forcesX[i] += fx;
	firstmsg->forcesY[i] += fy;
	firstmsg->forcesZ[i] += fz;
      }
    }
  }
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  delete first;
}

