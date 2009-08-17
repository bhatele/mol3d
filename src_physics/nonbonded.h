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

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ vdwParams *vdwTable;

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

//extern /* readonly */ double A;			// Force Calculation parameter 1
//extern /* readonly */ double B;			// Force Calculation parameter 2

//calculate pair forces using pairlists
inline CkVec<int>* calcPairForcesPL(ParticleDataMsg* first, ParticleDataMsg* second, CkVec<int> *pairList, int *numLists, CkSectionInfo* cookie1, CkSectionInfo* cookie2) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  double A, B;
  vdwPars *vdwp;
  sqrtPars pars;
  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = patchSizeX * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->coords[i].x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = patchSizeY * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->coords[i].y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = patchSizeZ * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->coords[i].z += diff;
  } 
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10, -20);

  //check if pairlist needs to be updated
  if (first->updateList){
    pairList = new CkVec<int>[firstLen];
    *numLists = firstLen;
    for (i = 0; i < firstLen; i++){
      firstx = first->coords[i].x;
      firsty = first->coords[i].y;
      firstz = first->coords[i].z;
      for (j = 0; j < secondLen; j++){
        rx = firstx - second->coords[j].x;
        ry = firsty - second->coords[j].y;
        rz = firstz - second->coords[j].z;
	rsqd = rx*rx + ry*ry + rz*rz;
	if (rsqd >= 0.001 && rsqd < ptpCutOffSqd)
	  pairList[i].push_back(j);
      }
    }
  }

  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
 
  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  memset(secondmsg->forces, 0, secondLen * 3*sizeof(BigReal));
  for(i = 0; i < firstLen; i++){
    eField = first->charge[i];
    firstx = first->coords[i].x;
    firsty = first->coords[i].y;
    firstz = first->coords[i].z;
    for(j = 0; j < pairList[i].length(); j++) {
      jpart = pairList[i][j];
      rx = firstx - second->coords[jpart].x;
      ry = firsty - second->coords[jpart].y;
      rz = firstz - second->coords[jpart].z;
      rsqd = rx*rx + ry*ry + rz*rz;
      if (rsqd >= 0.001){
	rsqd = rsqd * powTwenty;
	//r = sqrt(rsqd);
	pars = rootTable->pars[(int)(rsqd/rootTable->delta)];
	r = pars.a + rsqd*(pars.b+rsqd*(pars.c+ rsqd*pars.d));
	vdwp = &(vdwTable->params[first->vdwIndex[i]*vdwTable->numParams + second->vdwIndex[jpart]]);
	A = vdwp->A;
	B = vdwp->B;
	//r = r * 10000000000;
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * second->charge[jpart] / rsqd;  //positive force should be attractive
	fr = f /r;
	fx = rx * fr;
	fy = ry * fr;
	fz = rz * fr;
	secondmsg->forces[jpart].x += fx;
	secondmsg->forces[jpart].y += fy;
	secondmsg->forces[jpart].z += fz;
	firstmsg->forces[i].x -= fx;
	firstmsg->forces[i].y -= fy;
	firstmsg->forces[i].z -= fz;
      }
    }
  }
#ifdef USE_SECTION_MULTICAST
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*cookie1, first);
  mCastGrp->contribute(sizeof(BigReal)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  CkGetSectionInfo(*cookie2, second);
  mCastGrp->contribute(sizeof(BigReal)*3*secondmsg->lengthUpdates, secondmsg->forces, CkReduction::sum_double, *cookie2);
  delete firstmsg;
  delete secondmsg;
#else
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
#endif
  if (first->deleteList){
    for (i = 0; i < firstLen; i++){
      pairList[i].removeAll();
    }
    delete [] pairList;
    *numLists = -1;
    pairList = NULL;
  }
  delete first;
  delete second;
  return pairList;
}

//calculate pair forces without using pairlists
inline void calcPairForces(ParticleDataMsg* first, ParticleDataMsg* second, CkSectionInfo* cookie1, CkSectionInfo* cookie2) {
  int i, j, jpart, ptpCutOffSqd, diff;
  int firstLen = first->lengthAll;
  int secondLen = second->lengthAll;
  BigReal powTwenty, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  double A, B;
  sqrtPars pars;
  vdwPars *vdwp;

  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  ParticleForceMsg *secondmsg = new (secondLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
  secondmsg->lengthUpdates = secondLen;
  //check for wrap around and adjust locations accordingly
  if (abs(first->x - second->x) > 1){
    diff = patchSizeX * patchArrayDimX;
    if (second->x < first->x)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->coords[i].x += diff;
  }
  if (abs(first->y - second->y) > 1){
    diff = patchSizeY * patchArrayDimY;
    if (second->y < first->y)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->coords[i].y += diff;
  }
  if (abs(first->z - second->z) > 1){
    diff = patchSizeZ * patchArrayDimZ;
    if (second->z < first->z)
      diff = -1 * diff; 
    for (i = 0; i < firstLen; i++)
      first->coords[i].z += diff;
  } 
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10, -20);
  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
 
  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  memset(secondmsg->forces, 0, secondLen * 3*sizeof(BigReal));
  for(i = 0; i < firstLen; i++){
    eField = first->charge[i];
    for(jpart = 0; jpart < secondLen; jpart++){
      rx = first->coords[i].x - second->coords[jpart].x;
      ry = first->coords[i].y - second->coords[jpart].y;
      rz = first->coords[i].z - second->coords[jpart].z;
      rsqd = rx*rx + ry*ry + rz*rz;
      if (rsqd >= 0.001 && rsqd < ptpCutOffSqd){
	rsqd = rsqd * powTwenty;
	//r = rsqd/1000000000;
	//r = sqrt(rsqd);
	pars = rootTable->pars[(int)(rsqd/rootTable->delta)];
	r = pars.a + rsqd*(pars.b+rsqd*(pars.c+ rsqd*pars.d));
	vdwp = &vdwTable->params[first->vdwIndex[i]*vdwTable->numParams + second->vdwIndex[jpart]];
	A = vdwp->A;
	B = vdwp->B;
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * second->charge[jpart] / rsqd;
	fr = f /r;
	fx = rx * fr;
	fy = ry * fr;
	fz = rz * fr;
	secondmsg->forces[jpart].x -= fx;
	secondmsg->forces[jpart].y -= fy;
	secondmsg->forces[jpart].z -= fz;
	firstmsg->forces[i].x += fx;
	firstmsg->forces[i].y += fy;
	firstmsg->forces[i].z += fz;
      }
    }
  }
#ifdef USE_SECTION_MULTICAST
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  CkGetSectionInfo(*cookie1, first);
  mCastGrp->contribute(sizeof(BigReal)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  CkGetSectionInfo(*cookie2, second);
  mCastGrp->contribute(sizeof(BigReal)*3*secondmsg->lengthUpdates, secondmsg->forces, CkReduction::sum_double, *cookie2);
  delete firstmsg;
  delete secondmsg;
#else
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
  patchArray(second->x, second->y, second->z).receiveForces(secondmsg);
#endif
  delete first;
  delete second;
}

// Local function to compute all the interactions between pairs of particles in two sets
inline void calcInternalForces(ParticleDataMsg* first, CkSectionInfo *cookie1) {
  int i, j, ptpCutOffSqd;
  int firstLen = first->lengthAll;
  BigReal powTwenty, firstx, firsty, firstz, rx, ry, rz, r, rsqd, fx, fy, fz, f, fr, eField, constants;
  double rSix, rTwelve;
  double A, B;
  vdwPars *vdwp;
  sqrtPars pars;

  ParticleForceMsg *firstmsg = new (firstLen) ParticleForceMsg;
  firstmsg->lengthUpdates = firstLen;
    
  constants = COULOMBS_CONSTANT * ELECTRON_CHARGE * ELECTRON_CHARGE;
  
  memset(firstmsg->forces, 0, firstLen * 3*sizeof(BigReal));
  ptpCutOffSqd = ptpCutOff * ptpCutOff;
  powTwenty = pow(10, -20);
  for(i = 0; i < firstLen; i++){
    eField = first->charge[i];
    firstx = first->coords[i].x;
    firsty = first->coords[i].y;
    firstz = first->coords[i].z;
    for(j = i+1; j < firstLen; j++) {
      // computing base values
      rx = firstx - first->coords[j].x;
      ry = firsty - first->coords[j].y;
      rz = firstz - first->coords[j].z;
      rsqd = rx*rx + ry*ry + rz*rz;
      //check if r >= .000001 to make sure force calc doesnt tend to 0
      if(rsqd >= 0.001 && rsqd < ptpCutOffSqd){
	rsqd = rsqd * powTwenty;
	//r = sqrt(rsqd);
	pars = rootTable->pars[(int)(rsqd/rootTable->delta)];
	r = pars.a + rsqd*(pars.b+rsqd*(pars.c+ rsqd*pars.d));
	vdwp = &vdwTable->params[first->vdwIndex[i]*vdwTable->numParams + first->vdwIndex[j]];
	A = vdwp->A;
	B = vdwp->B;
	rSix = ((double)rsqd) * rsqd * rsqd;
	rTwelve = rSix * rSix;
        f = (BigReal)(A / rTwelve - B / rSix);
	f -= eField * constants * first->charge[j] / (rsqd);
        
	fr = f /r;
	fx = rx * fr;
        fy = ry * fr;
        fz = rz * fr;
	firstmsg->forces[j].x += fx;
	firstmsg->forces[j].y += fy;
	firstmsg->forces[j].z += fz;
	firstmsg->forces[i].x -= fx;
	firstmsg->forces[i].y -= fy;
	firstmsg->forces[i].z -= fz;
      }
    }
  }
#ifdef USE_SECTION_MULTICAST
  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  //CkGetSectionInfo(*cookie1, first);
  //double might be incorrect here
  //double *forceArr = new double[3*firstmsg->lengthUpdates];
  //for (int i = 0; i < firstmsg->lengthUpdates; i++){
    //forceArr[3*i] = firstmsg->forces[i].x;
    //forceArr[3*i+1] = firstmsg->forces[i].y;
    //forceArr[3*i+2] = firstmsg->forces[i].z;
  //}
  //CkPrintf("lengthupdates = %d\n", firstmsg->lengthUpdates);
  mCastGrp->contribute(sizeof(BigReal)*3*firstmsg->lengthUpdates, firstmsg->forces, CkReduction::sum_double, *cookie1);
  delete firstmsg;
#else
  patchArray(first->x, first->y, first->z).receiveForces(firstmsg);
#endif
  delete first;
}

