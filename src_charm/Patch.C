/** \file Patch.C
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */

#include "time.h"
#include "defs.h"
#include "mol3d.decl.h"
#include "Main.h"
#include "Patch.h"
#include "Compute.h"
#include "ConfigList.h"
#include "molfile_plugin.h"
#include "PluginIOMgr.h"
#ifdef USE_SECTION_MULTICAST
  #include "ckmulticast.h"
#endif

  
extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;	// Number of Chares in X
extern /* readonly */ int patchArrayDimY;	// Number of Chares in Y
extern /* readonly */ int patchArrayDimZ;	// Number of Chares in Z
extern /* readonly */ int patchSizeX;
extern /* readonly */ int patchSizeY;
extern /* readonly */ int patchSizeZ;
extern /* readonly */ int ptpCutOff;
extern /* readonly */ int patchMargin;
extern /* readonly */ int patchOriginX;
extern /* readonly */ int patchOriginY;
extern /* readonly */ int patchOriginZ;
extern /* readonly */ int migrateStepCount;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ int firstLdbStep; 
extern /* readonly */ int ldbPeriod; 
extern /* readonly */ int ftPeriod; 
extern /* readonly */ BigReal stepTime; 
extern /* readonly */ BigReal timeDelta;
extern /* readonly */ bool usePairLists;
extern /* readonly */ bool twoAwayX;
extern /* readonly */ int numNbrs;
extern /* readonly */ int nbrsX;
extern /* readonly */ int nbrsY;
extern /* readonly */ int nbrsZ;

extern /* readonly */ double A;
extern /* readonly */ double B;

// Default constructor
Patch::Patch(FileDataMsg* fdmsg) {
  LBTurnInstrumentOff();
  int i;
  inbrs = numNbrs;
  usesAtSync = CmiTrue;
  //double d1 = CmiWallTimer();
  // Particle initialization
  myNumParts = 0;
  for(i=0; i < fdmsg->length; i++) {
    particles.push_back(Particle());
    particles[myNumParts].charge = fdmsg->charge[i];
    particles[myNumParts].mass = fdmsg->mass[i];

    particles[myNumParts].x = fdmsg->coords[i].x;
    particles[myNumParts].y = fdmsg->coords[i].y;
    particles[myNumParts].z = fdmsg->coords[i].z;

    particles[myNumParts].vx = 0;
    particles[myNumParts].vy = 0;
    particles[myNumParts].vz = 0;
    particles[myNumParts].fx = 0;
    particles[myNumParts].fy = 0;
    particles[myNumParts].fz = 0;
     
    particles[myNumParts].id = (thisIndex.x*patchArrayDimX + thisIndex.y) * numParts / (patchArrayDimX*patchArrayDimY)  + i;
    
    particles[myNumParts].vdw_type = fdmsg->vdw_type[i];
    myNumParts++;
  }	
  //CkPrintf("Creating %d particles on Patch [%d][%d][%d]\n", myNumParts, thisIndex.x, thisIndex.y, thisIndex.z);


  updateCount = 0;
  forceCount = 0;
  stepCount = 0;
  resumeCount = 0;
  updateFlag = false;
  incomingFlag = false;
  pause = false;
  incomingParticles.resize(0);
 // setMigratable(CmiFalse);
  delete fdmsg;
  //loadTime = CmiWallTimer() - d1;
}

// Constructor for chare object migration
Patch::Patch(CkMigrateMessage *msg): CBase_Patch(msg) { 
  usesAtSync = CmiTrue;
  delete msg;
}  
                                       
Patch::~Patch() {}


void Patch::createComputes() {
  //double d1 = CmiWallTimer();
  int num;  
  
  int x = thisIndex.x;
  int y = thisIndex.y;
  int z = thisIndex.z;
  int px1, py1, pz1, dx, dy, dz, px2, py2, pz2;

  // For Round Robin insertion
  int numPes = CkNumPes();
  int currPe = CkMyPe();

  computesList = new int*[numNbrs];
  for (int i =0; i < numNbrs; i++){
    computesList[i] = new int[6];
  }
 
  /*  The computes X are inserted by a given patch:
   *
   *	^  X  X  X
   *	|  0  X  X
   *	y  0  0  0
   *	   x ---->
   */

  // these computes will be created by other patches
  for (num=0; num<numNbrs; num++) {
    dx = num / (nbrsY * nbrsZ)            - nbrsX/2;
    dy = (num % (nbrsY * nbrsZ)) / nbrsZ - nbrsY/2;
    dz = num % nbrsZ                       - nbrsZ/2;

    if (num >= numNbrs/2){
      px1 = x + 2;
      px2 = x+dx+2;
      py1 = y + 2;
      py2 = y+dy+2;
      pz1 = z + 2;
      pz2 = z+dz+2;
      computeArray(px1, py1, pz1, px2, py2, pz2).insert((++currPe)%numPes);
      computesList[num][0] = px1; computesList[num][1] = py1; computesList[num][2] = pz1; 
      computesList[num][3] = px2; computesList[num][4] = py2; computesList[num][5] = pz2;
    }
    else {
      px2 = WRAP_X(x+dx);
      py2 = WRAP_Y(y+dy);
      pz2 = WRAP_Z(z+dz);
      px1 = x;
      py1 = y;
      pz1 = z; 
      px1 = px2 - dx + 2;
      px2 = px2+2;
      py1 = py2 - dy + 2;
      py2 = py2+2;
      pz1 = pz2 - dz + 2;
      pz2 = pz2+2;
      computesList[num][0] = px2; computesList[num][1] = py2; computesList[num][2] = pz2; 
      computesList[num][3] = px1; computesList[num][4] = py1; computesList[num][5] = pz1;
    }

    //insert only the upper right half computes
  } // end of for loop

  contribute(CkCallback(CkIndex_Main::startUpDone(), mainProxy));
  //loadTime += CmiWallTimer()-d1;
}

void Patch::createSection() {
  localCreateSection();
  contribute(CkCallback(CkIndex_Main::startUpDone(), mainProxy));
}

void Patch::localCreateSection() {
#ifdef USE_SECTION_MULTICAST
  CkVec<CkArrayIndex6D> elems;
  for (int num=0; num<numNbrs; num++)
    elems.push_back(CkArrayIndex6D(computesList[num][0], computesList[num][1], computesList[num][2], computesList[num][3], computesList[num][4], computesList[num][5]));

  CkArrayID computeArrayID = computeArray.ckGetArrayID();
  mCastSecProxy = CProxySection_Compute::ckNew(computeArrayID, elems.getVec(), elems.size()); 

  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);
  mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkIndex_Patch::reduceForces(NULL), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));

#endif
}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Patch::start() {
  int x = thisIndex.x;
  int y = thisIndex.y;
  int z = thisIndex.z;
  int len = particles.length();
  
  if (stepCount == 0 && x+y+z ==0)
    stepTime = CmiWallTimer();

  ParticleDataMsg* msg = new (len) ParticleDataMsg;
  msg->x = x;
  msg->y = y;
  msg->z = z;
  msg->lengthAll = len;
  msg->deleteList = false;
  msg->updateList = false;
  msg->doAtSync = false;
  // If using pairlists determine whether or not its time to update the pairlist
  if (usePairLists){
    // set up pairlsit at startup
    if (stepCount == 0)
      msg->updateList = true;
    if (stepCount > 1) {
      // delete pairlist if about to migrate
      if (stepCount % migrateStepCount == 0)
	msg->deleteList = true;
      // rebild pairlist once migrated
      if (stepCount % migrateStepCount == 1)
	msg->updateList = true;
      // delete pairlist if load balancing is about to be done
      if ((stepCount - firstLdbStep) % ldbPeriod == 0)
	msg->deleteList = true;
      // rebuild pairlist if load balancing was just done
      if ((stepCount - firstLdbStep) % ldbPeriod == 1)
	msg->updateList = true;
    }
  }
#ifdef USE_SECTION_MULTICAST
  // if we are using section mutlicast and we just did migration we need to rebuild the section
  if (stepCount > 1 && (stepCount - firstLdbStep) % ldbPeriod == 1){
    CkVec<CkArrayIndex6D> elems;
    for (int num=0; num<numNbrs; num++)
      elems.push_back(CkArrayIndex6D(computesList[num][0], computesList[num][1], computesList[num][2], computesList[num][3], computesList[num][4], computesList[num][5]));

    CkArrayID computeArrayID = computeArray.ckGetArrayID();
    mCastSecProxy = CProxySection_Compute::ckNew(computeArrayID, elems.getVec(), elems.size());
    
//	CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
    //mCastSecProxy.ckSectionDelegate(mCastGrp);
//	mCastGrp->resetSection(mCastSecProxy);
//	mCastGrp->setReductionClient(mCastSecProxy, new CkCallback(CkIndex_Patch::reduceForces(NULL), thisProxy(thisIndex.x, thisIndex.y, thisIndex.z)));
  
  }
#endif
  msg->lbOn = false;
  if (((stepCount > firstLdbStep - 1) && stepCount % ldbPeriod == 1) || stepCount == 0){
    // if (x + y + z == 0) CkPrintf("Starting Load Balancer Instrumentation at %f\n", CmiWallTimer());
    msg->lbOn = true;
  }
  if ((stepCount - firstLdbStep) % ldbPeriod == 0){
    msg->doAtSync = true;
    pause = true;
  }
  for (int i = 0; i < len; i++){
    msg->part[i].coord.x = particles[i].x;
    msg->part[i].coord.y = particles[i].y;
    msg->part[i].coord.z = particles[i].z;
    msg->part[i].charge = particles[i].charge;
    msg->part[i].vdwIndex = particles[i].vdw_type;
  }
  //if (x+y+z < 1)
    //  CkPrintf("patch %d %d %d, updatelist = %d deletelist = %d\n",x,y,z,msg->updateList, msg->deleteList);
#ifdef USE_SECTION_MULTICAST
  mCastSecProxy.interact(msg);
#else
  int px1, py1, pz1, px2, py2, pz2;
  
  for(int num=0; num<numNbrs; num++) {
    px1 = computesList[num][0];
    py1 = computesList[num][1];
    pz1 = computesList[num][2];
    px2 = computesList[num][3];
    py2 = computesList[num][4];
    pz2 = computesList[num][5];
    if (num == numNbrs-1)
      computeArray(px1, py1, pz1, px2, py2, pz2).interact(msg);
    else {
      ParticleDataMsg* newMsg = new (len) ParticleDataMsg;
      newMsg->x = x;
      newMsg->y = y;
      newMsg->z = z;
      newMsg->lengthAll = len;
      if (usePairLists){
	newMsg->updateList = msg->updateList;
	newMsg->deleteList = msg->deleteList;
      }
      newMsg->doAtSync = msg->doAtSync;
      newMsg->lbOn = msg->lbOn;
      memcpy(newMsg->part, msg->part, len*sizeof(partData));
      //memcpy(newMsg->charge, msg->charge, len*sizeof(BigReal));
      //memcpy(newMsg->vdwIndex, msg->vdwIndex, len*sizeof(int));
      computeArray(px1, py1, pz1, px2, py2, pz2).interact(newMsg);
    } 
  }
#endif
}

//reduction to update forces coming from a compute
void Patch::reduceForces(CkReductionMsg *msg) {
  //double d1 = CmiWallTimer();
  int i, lengthUp;
  forceCount=numNbrs;
  int* forces = (int*)msg->getData();
  lengthUp = msg->getSize()/sizeof(BigReal);
  //CkPrintf("lengthup = %d numparts = %d\n",lengthUp, particles.length());
  for(i = 0; i < lengthUp; i+=3){
    particles[i/3].fx += forces[i];
    particles[i/3].fy += forces[i+1];
    particles[i/3].fz += forces[i+2];
  }
  applyForces();
  delete msg;
  //loadTime += CmiWallTimer()-d1;
}



// Function to update forces coming from a compute
void Patch::receiveForces(ParticleForceMsg *updates) {
  int i;
  // incrementing the counter for receiving updates
  forceCount++;

  // updating force information
  for(i = 0; i < updates->lengthUpdates; i++){
    particles[i].fx += updates->forces[i].x;
    particles[i].fy += updates->forces[i].y;
    particles[i].fz += updates->forces[i].z;
  }
  delete updates;
  applyForces();
}


void Patch::applyForces(){
  int i, x, y, z, x1, y1, z1;
  // if all forces are received, then it must recompute particles location
  if (forceCount == numNbrs) {
    CkVec<Particle> *outgoing = new CkVec<Particle>[numNbrs];

    // Received all it's forces from the interactions.
    forceCount = 0;
  
    // Update properties on own particles
    updateProperties();

    // sending particles to neighboring cells
    x = thisIndex.x;
    y = thisIndex.y;
    z = thisIndex.z;
    if (stepCount > 0 && (stepCount % migrateStepCount) == 0){
      for(i=0; i<particles.length(); i++) {
	migrateToPatch(particles[i], x1, y1, z1);
	if(x1 !=0 || y1!=0 || z1 !=0) {
	  //CkPrintf("PARTICLE MIGRATING!\n");
	  outgoing[(x1+1)*nbrsY*nbrsZ + (y1+1)*nbrsZ + (z1+1)].push_back(wrapAround(particles[i]));
	  particles.remove(i);
	}
      }
    
   
      for(int num=0; num<numNbrs; num++) {
	x1 = num / (nbrsY * nbrsZ)            - nbrsX/2;
	y1 = (num % (nbrsY * nbrsZ)) / nbrsZ - nbrsY/2;
	z1 = num % nbrsZ                       - nbrsZ/2;

	patchArray(WRAP_X(x+x1), WRAP_Y(y+y1), WRAP_Z(z+z1)).receiveParticles(outgoing[num]);
      }
    }
    else
      incomingFlag = true;

    updateFlag = true;
	      
    // checking whether to proceed with next step
    thisProxy(x, y, z).checkNextStep();
    //checkNextStep();
    delete [] outgoing;
  //  thisProxy(x, y, z).checkNextStep();
    //checkNextStep();
  }
//  else { CkPrintf("forcecount = %d/%d on patch %d %d %d\n", forceCount, numNbrs, thisIndex.x, thisIndex.y, thisIndex.z); }

}

void Patch::migrateToPatch(Particle p, int &px, int &py, int &pz) {
  // currently this is assuming that particles are
  // migrating only to the immediate neighbors
  int x = thisIndex.x * patchSizeX + patchOriginX;
  int y = thisIndex.y * patchSizeY + patchOriginY;
  int z = thisIndex.z * patchSizeZ + patchOriginZ;

  if (p.x < x) px = -1;
  else if (p.x > x+patchSizeX) px = 1;
  else px = 0;

  if (p.y < y) py = -1;
  else if (p.y > y+patchSizeY) py = 1;
  else py = 0;

  if (p.z < z) pz = -1;
  else if (p.z > z+patchSizeZ) pz = 1;
  else pz = 0;
}

// Function that checks whether it must start the following step or wait until other messages are received
void Patch::checkNextStep(){
  int i;
  double timer;

  if (updateFlag && incomingFlag) {
    // resetting flags
    updateFlag = false;
    incomingFlag = false;
    stepCount++;

    // adding new elements
    for (i = 0; i < incomingParticles.length(); i++)
      particles.push_back(incomingParticles[i]);
    incomingParticles.removeAll();

    if (thisIndex.x==0 && thisIndex.y==0 && thisIndex.z==0 && stepCount%20==0) {
      timer = CmiWallTimer();
      CkPrintf("Step %d Benchmark Time %f ms/step, Total Time Elapsed %f ms\n", stepCount, ((timer - stepTime)/20)*1000, timer);
      stepTime = timer;
//      if (stepCount == 300)
//	traceBegin();
  //    if (stepCount == 400)
//	traceEnd();

    }
 //   if (stepCount == 300 && thisIndex.x*patchArrayDimY*patchArrayDimZ + thisIndex.y*patchArrayDimZ + thisIndex.z < 8)
 //     traceBegin();
 //   if (stepCount == 301 && thisIndex.x*patchArrayDimY*patchArrayDimZ + thisIndex.y*patchArrayDimZ + thisIndex.z < 8)
 //     traceEnd();

    // checking for next step
    if (stepCount >= finalStepCount) {
     // CkPrintf("Final number of particles is %d on Patch [%d][%d][%d]\n", particles.length(), thisIndex.x, thisIndex.y, thisIndex.z);
      print();
      contribute(CkCallback(CkIndex_Main::allDone(), mainProxy)); 
    } else {
      if (!pause){
	if ((stepCount > firstLdbStep - 1) && stepCount % ldbPeriod == 0){
	  // if (x + y + z == 0) CkPrintf("Starting Load Balancer Instrumentation at %f\n", CmiWallTimer());
	  contribute(CkCallback(CkIndex_Main::lbBarrier(),mainProxy));
	  return;
	}
        if (stepCount % ftPeriod == 1 && stepCount > 1) {
	  contribute(CkCallback(CkIndex_Main::ftBarrier(),mainProxy));
	  return;
        }
	thisProxy(thisIndex.x, thisIndex.y, thisIndex.z).start();
      }
      else{
	AtSync();
//	loadTime += CmiWallTimer()-d1;
	//CkPrintf("Patch %d on processor %d had load %f\n", thisIndex.x*patchArrayDimY*patchArrayDimZ + thisIndex.y*patchArrayDimZ + thisIndex.z, CkMyPe(), loadTime);
//	loadTime = 0;
	contribute(CkCallback(CkIndex_Main::lbBarrier(),mainProxy));
      }
	//loadTime = 0;
    }
  }
}

void Patch::ResumeFromSync(){
    // if (thisIndex.x+thisIndex.y+thisIndex.z == 0) CkPrintf("patch 0 ResumeFromSync at %f\n",CmiWallTimer());
    pause = false;
    stepTime = CmiWallTimer();
    thisProxy(thisIndex.x, thisIndex.y, thisIndex.z).start();
    LBTurnInstrumentOff();
    //resumeCount = 0;
}

void Patch::resume(){
 /* if (++resumeCount == numNbrs){
    pause = false;
    thisProxy(thisIndex.x, thisIndex.y, thisIndex.z).start();
    resumeCount = 0;
  }*/
  if (!pause){
    // if (thisIndex.x+thisIndex.y+thisIndex.z == 0) CkPrintf("patch 0 calling LBInstrumentation on at %f\n",CmiWallTimer());
    LBTurnInstrumentOn();
    loadTime = 0;
    start();
  }
  else {
    // if (thisIndex.x+thisIndex.y+thisIndex.z == 0) CkPrintf("patch 0 calling AtSync at %f\n",CmiWallTimer());
    AtSync();
  }
}

void Patch::ftresume(){
  if (thisIndex.x==0 && thisIndex.y==0 && thisIndex.z ==0)
      CkPrintf("patch 0 calling ftresume at %f\n",CmiWallTimer());
  start();
}

// Function that receives a set of particles and updates the 
// forces of them into the local set
void Patch::receiveParticles(CkVec<Particle> &updates) {
  updateCount++;

  for( int i=0; i < updates.length(); i++) {
    incomingParticles.push_back(updates[i]);
  }

  // if all the incoming particle updates have been received, we must check 
  // whether to proceed with next step
  if(updateCount == numNbrs ) {
    updateCount = 0;
    incomingFlag = true;
    checkNextStep();
  }
}

// Function to update properties (i.e. acceleration, velocity and position) in particles
void Patch::updateProperties() {
  int i;
  BigReal powTen, powFteen, realTimeDelta, invMassParticle;
  powTen = pow(10.0, -10);
  powFteen = pow(10.0, -15);
  realTimeDelta = timeDelta * powFteen;
  for(i = 0; i < particles.length(); i++) {
    // applying kinetic equations
    invMassParticle = (AVAGADROS_NUMBER / (particles[i].mass * powTen));
    particles[i].ax = particles[i].fx * invMassParticle;
    particles[i].ay = particles[i].fy * invMassParticle;
    particles[i].az = particles[i].fz * invMassParticle;
    particles[i].vx = particles[i].vx + particles[i].ax * realTimeDelta;
    particles[i].vy = particles[i].vy + particles[i].ay * realTimeDelta;
    particles[i].vz = particles[i].vz + particles[i].az * realTimeDelta;

    limitVelocity(particles[i]);

    particles[i].x = particles[i].x + particles[i].vx * realTimeDelta;
    particles[i].y = particles[i].y + particles[i].vy * realTimeDelta;
    particles[i].z = particles[i].z + particles[i].vz * realTimeDelta;

    particles[i].fx = 0.0;
    particles[i].fy = 0.0;
    particles[i].fz = 0.0;
  }
}

void Patch::limitVelocity(Particle &p) {
  if( fabs( p.vx ) > MAX_VELOCITY ) {
    if( p.vx < 0.0 )
      p.vx = -MAX_VELOCITY;
    else
      p.vx = MAX_VELOCITY;
  }

  if( fabs(p.vy) > MAX_VELOCITY ) {
    if( p.vy < 0.0 )
      p.vy = -MAX_VELOCITY;
    else
      p.vy = MAX_VELOCITY;
  }

  if( fabs(p.vz) > MAX_VELOCITY ) {
    if( p.vz < 0.0 )
      p.vz = -MAX_VELOCITY;
    else
      p.vz = MAX_VELOCITY;
  }
}

Particle& Patch::wrapAround(Particle &p) {
  if(p.x < patchOriginX) p.x += patchSizeX*patchArrayDimX;
  if(p.y < patchOriginY) p.y += patchSizeY*patchArrayDimY;
  if(p.z < patchOriginZ) p.z += patchSizeZ*patchArrayDimZ;
  if(p.x > patchOriginX + patchSizeX*patchArrayDimX) p.x -= patchSizeX*patchArrayDimX;
  if(p.y > patchOriginY + patchSizeY*patchArrayDimY) p.y -= patchSizeY*patchArrayDimY;
  if(p.z > patchOriginZ + patchSizeZ*patchArrayDimZ) p.z -= patchSizeZ*patchArrayDimZ;

  return p;
}

// Prints all particles 
void Patch::print(){
#ifdef PRINT
  int i;
  CkPrintf("*****************************************************\n");
  CkPrintf("Patch (%d, %d)\n", thisIndex.x, thisIndex.y);

  for(i=0; i < particles.length(); i++)
    CkPrintf("Patch (%d,%d) %-5d %7.4f %7.4f \n", thisIndex.x, thisIndex.y, i, particles[i].x, particles[i].y);
  CkPrintf("*****************************************************\n");
#endif
}

