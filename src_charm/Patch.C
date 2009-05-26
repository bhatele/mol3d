/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file Patch.C
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */

#include "time.h"
#include "defs.h"
#ifdef RUN_LIVEVIZ
  #include "liveViz.h"
#endif
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
//extern /* readonly */ CProxy_GridCompute gridComputeArray;
//extern /* readonly */ CProxy_PMECompositor PMECompArray;
//extern /* readonly */ CProxy_PMEDecompositor PMEDecompArray;
extern /* readonly */ CkGroupID mCastGrpID;


extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;	// Number of Chares in X
extern /* readonly */ int patchArrayDimY;	// Number of Chares in Y
extern /* readonly */ int patchArrayDimZ;	// Number of Chares in Z
extern /* readonly */ int patchSize;
extern /* readonly */ int ptpCutOff;
extern /* readonly */ int patchMargin;
extern /* readonly */ int patchOriginX;
extern /* readonly */ int patchOriginY;
extern /* readonly */ int patchOriginZ;
extern /* readonly */ int compArrayLenX;
extern /* readonly */ int compArrayLenY;
extern /* readonly */ int compArrayLenZ;
extern /* readonly */ int pmeGridDimX;
extern /* readonly */ int pmeGridDimY;
extern /* readonly */ int pmeGridDimZ;
extern /* readonly */ BigReal pmeCellSize;
extern /* readonly */ int pmeCutOff;
extern /* readonly */ int migrateStepCount;
extern /* readonly */ int finalStepCount; 
extern /* readonly */ BigReal stepTime; 
extern /* readonly */ BigReal timeDelta;
extern /* readonly */ bool usePairLists;

extern /* readonly */ double A;
extern /* readonly */ double B;

// Default constructor
Patch::Patch(){}
//Patch::Patch(FileDataMsg* fdmsg) {
void Patch::Iinsert(FileDataMsg* fdmsg){
  int i;

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
  updateFlag = false;
  incomingFlag = false;
  incomingParticles.resize(0);
  setMigratable(CmiFalse);
  delete fdmsg;
}

// Constructor for chare object migration
Patch::Patch(CkMigrateMessage *msg) { }  
                                       
Patch::~Patch() {}


void Patch::createComputes() {
  int num;  
  
  int x = thisIndex.x;
  int y = thisIndex.y;
  int z = thisIndex.z;
  int px1, py1, pz1, dx, dy, dz, px2, py2, pz2;

  // For Round Robin insertion
  int numPes = CkNumPes();
  int currPe = CkMyPe();
 
  /*  The computes X are inserted by a given patch:
   *
   *	^  X  X  X
   *	|  0  X  X
   *	y  0  0  0
   *	   x ---->
   */

  // these computes will be created by other patches
  for (num=0; num<NUM_NEIGHBORS; num++) {
    dx = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
    dy = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z - NBRS_Y/2;
    dz = num % NBRS_Z                       - NBRS_Z/2;

    if (dx == 0) {
      px1 = px2 = x;

      if (dy == 0) { 
	py1 = py2 = y;
	if (dz == 0) { pz1 = pz2 = z; }
	if (dz > 0) { (z >= patchArrayDimZ - NBRS_Z/2) ? ( pz1 = WRAP_Z(z+dz), pz2 = z ) : ( pz1 = z, pz2 = z+dz ); }
	if (dz < 0) { (z < NBRS_Z/2) ? ( pz1 = z, pz2 = WRAP_Z(z+dz) ) : ( pz1 = z+dz, pz2 = z ); }
      }

      if (dy > 0) { 
	(y >= patchArrayDimY - NBRS_Y/2) ? 
	( py1 = WRAP_Y(y+dy), pz1 = WRAP_Z(z+dz), py2 = y, pz2 = z ) : 
	( py1 = y, pz1 = z, py2 = y+dy, pz2 = WRAP_Z(z+dz) );
      }

      if (dy < 0) { 
	(y < NBRS_Y/2) ? 
	( py1 = y, pz1 = z, py2 = WRAP_Y(y+dy), pz2 = WRAP_Z(z+dz) ) : 
	( py1 = y+dy, pz1 = WRAP_Z(z+dz), py2 = y, pz2 = z ); 
      }
    } // dx == 0

    if (dx > 0) {
      (x >= patchArrayDimX - NBRS_X/2) ? 
      ( px1 = WRAP_X(x+dx), py1 = WRAP_Y(y+dy), pz1 = WRAP_Z(z+dz), px2 = x, py2 = y, pz2 = z ) : 
      ( px1 = x, py1 = y, pz1 = z, px2 = WRAP_X(x+dx), py2 = WRAP_Y(y+dy), pz2 = WRAP_Z(z+dz) ) ;
    }

    if (dx < 0) {
      (x < NBRS_X/2) ? 
      ( px1 = x, py1 = y, pz1 = z, px2 = WRAP_X(x+dx), py2 = WRAP_Y(y+dy), pz2 = WRAP_Z(z+dz) ) :
      ( px1 = WRAP_X(x+dx), py1 = WRAP_Y(y+dy), pz1 = WRAP_Z(z+dz), px2 = x, py2 = y, pz2 = z ) ;
    }

    computesList[num][0] = px1; computesList[num][1] = py1; computesList[num][2] = pz1; 
    computesList[num][3] = px2; computesList[num][4] = py2; computesList[num][5] = pz2;

    //insert only the upper right half computes
    if (num >= NUM_NEIGHBORS/2)
      computeArray(px1, py1, pz1, px2, py2, pz2).insert((++currPe) % numPes);
  } // end of for loop

  contribute(CkCallback(CkIndex_Main::startUpDone(), mainProxy));
}

void Patch::createSection() {
#ifdef USE_SECTION_MULTICAST
  CkVec<CkArrayIndex6D> elems;
  for (int num=0; num<NUM_NEIGHBORS; num++)
    elems.push_back(CkArrayIndex6D(computesList[num][0], computesList[num][1], computesList[num][2], computesList[num][3], computesList[num][4], computesList[num][5]));

  CkArrayID computeArrayID = computeArray.ckGetArrayID();
  mCastSecProxy = CProxySection_Compute::ckNew(computeArrayID, elems.getVec(), elems.size()); 

  CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(mCastGrpID).ckLocalBranch();
  mCastSecProxy.ckSectionDelegate(mCastGrp);

  contribute(CkCallback(CkIndex_Main::startUpDone(), mainProxy));
	//CkPrintf("%E \n", f);
#endif
}

// Function to start interaction among particles in neighboring cells as well as its own particles
void Patch::start() {
  int x = thisIndex.x;
  int y = thisIndex.y;
  int z = thisIndex.z;
  int len = particles.length();
  
  ParticleDataMsg* msg = new (len, len, len) ParticleDataMsg;
  msg->x = x;
  msg->y = y;
  msg->z = z;
  msg->lengthAll = len;
  // If using pairlists determine whether or not its time to update the pairlist
  if (usePairLists){
    msg->deleteList = false;
    msg->doAtSync = false;
    if (stepCount == 0 || ((stepCount % migrateStepCount == 1) && stepCount > 1)){
      msg->updateList = true;
    }
    else{
      msg->updateList = false;
      if (stepCount % migrateStepCount == 0){
	msg->deleteList = true;
	if (migrateStepCount*1024 % stepCount == 0)
	  msg->doAtSync = true;
      }
    }
  }
  for (int i = 0; i < len; i++){
    msg->coords[i].x = particles[i].x;
    msg->coords[i].y = particles[i].y;
    msg->coords[i].z = particles[i].z;
    msg->charge[i] = particles[i].charge;
    msg->vdwIndex[i] = particles[i].vdw_type;
  }
#ifdef USE_SECTION_MULTICAST
  mCastSecProxy.interact(msg);
#else
  int px1, py1, pz1, px2, py2, pz2;
  
  for(int num=0; num<NUM_NEIGHBORS; num++) {
    px1 = computesList[num][0];
    py1 = computesList[num][1];
    pz1 = computesList[num][2];
    px2 = computesList[num][3];
    py2 = computesList[num][4];
    pz2 = computesList[num][5];
    if (num == NUM_NEIGHBORS-1)
      computeArray(px1, py1, pz1, px2, py2, pz2).interact(msg);
    else {
      ParticleDataMsg* newMsg = new (len, len, len) ParticleDataMsg;
      newMsg->x = x;
      newMsg->y = y;
      newMsg->z = z;
      newMsg->lengthAll = len;
      if (usePairLists){
	newMsg->updateList = msg->updateList;
	newMsg->deleteList = msg->deleteList;
	newMsg->doAtSync = msg->doAtSync;
      }
      memcpy(newMsg->coords, msg->coords, 3*len*sizeof(BigReal));
      memcpy(newMsg->charge, msg->charge, len*sizeof(BigReal));
      memcpy(newMsg->vdwIndex, msg->vdwIndex, len*sizeof(int));
      computeArray(px1, py1, pz1, px2, py2, pz2).interact(newMsg);
    } 
  }
#endif
}

// Function to update forces coming from a compute
void Patch::receiveForces(ParticleForceMsg *updates) {
  int i, x, y, z, x1, y1, z1;
  // incrementing the counter for receiving updates
  forceCount++;

  // updating force information
  for(i = 0; i < updates->lengthUpdates; i++){
    particles[i].fx += updates->forces[i].x;
    particles[i].fy += updates->forces[i].y;
    particles[i].fz += updates->forces[i].z;
  }

  // if all forces are received, then it must recompute particles location
  if (forceCount == NUM_NEIGHBORS) {
    CkVec<Particle> outgoing[NUM_NEIGHBORS];

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
	  outgoing[(x1+1)*NBRS_Y*NBRS_Z + (y1+1)*NBRS_Z + (z1+1)].push_back(wrapAround(particles[i]));
	  particles.remove(i);
	}
      }
    
   
      for(int num=0; num<NUM_NEIGHBORS; num++) {
	x1 = num / (NBRS_Y * NBRS_Z)            - NBRS_X/2;
	y1 = (num % (NBRS_Y * NBRS_Z)) / NBRS_Z - NBRS_Y/2;
	z1 = num % NBRS_Z                       - NBRS_Z/2;

	patchArray(WRAP_X(x+x1), WRAP_Y(y+y1), WRAP_Z(z+z1)).receiveParticles(outgoing[num]);
      }
    }
    else
      incomingFlag = true;

    updateFlag = true;
	      
    // checking whether to proceed with next step
    thisProxy(x, y, z).checkNextStep();
  }

  delete updates;
}

void Patch::migrateToPatch(Particle p, int &px, int &py, int &pz) {
  // currently this is assuming that particles are
  // migrating only to the immediate neighbors
  int x = thisIndex.x * patchSize + patchOriginX;
  int y = thisIndex.y * patchSize + patchOriginY;
  int z = thisIndex.z * patchSize + patchOriginZ;

  if (p.x < x) px = -1;
  else if (p.x > x+patchSize) px = 1;
  else px = 0;

  if (p.y < y) py = -1;
  else if (p.y > y+patchSize) py = 1;
  else py = 0;

  if (p.z < z) pz = -1;
  else if (p.z > z+patchSize) pz = 1;
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

    if (thisIndex.x==0 && thisIndex.y==0 && thisIndex.z==0 && stepCount%100==0) {
      timer = CmiWallTimer();
      CkPrintf("Step %d Benchmark Time %f ms/step\n", stepCount, ((timer - stepTime)/100)*1000);
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
      thisProxy(thisIndex.x, thisIndex.y, thisIndex.z).start();
    }
  }
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
  if(updateCount == NUM_NEIGHBORS ) {
    updateCount = 0;
    incomingFlag = true;
    checkNextStep();
  }
}

// Function to update properties (i.e. acceleration, velocity and position) in particles
void Patch::updateProperties() {
  int i;
  BigReal powTen, powFteen, realTimeDelta, invMassParticle;
  powTen = pow(10, -10);
  powFteen = pow(10, -15);
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
  if(p.x < patchOriginX) p.x += patchSize*patchArrayDimX;
  if(p.y < patchOriginY) p.y += patchSize*patchArrayDimY;
  if(p.z < patchOriginZ) p.z += patchSize*patchArrayDimZ;
  if(p.x > patchOriginX + patchSize*patchArrayDimX) p.x -= patchSize*patchArrayDimX;
  if(p.y > patchOriginY + patchSize*patchArrayDimY) p.y -= patchSize*patchArrayDimY;
  if(p.z > patchOriginZ + patchSize*patchArrayDimZ) p.z -= patchSize*patchArrayDimZ;

  return p;
}

// Helper function to help with LiveViz
void color_pixel(unsigned char*buf,int width, int height, int xpos,int ypos,
                             unsigned char R,unsigned char G,unsigned char B) {
  if(xpos>=0 && xpos<width && ypos>=0 && ypos<height) {
    buf[3*(ypos*width+xpos)+0] = R; 
    buf[3*(ypos*width+xpos)+1] = G; 
    buf[3*(ypos*width+xpos)+2] = B; 
  }
}
    
#ifdef RUN_LIVEVIZ
// Each chare provides its particle data to LiveViz
void Patch::requestNextFrame(liveVizRequestMsg *lvmsg) {
  // These specify the desired total image size requested by the client viewer
  int wdes = lvmsg->req.wid;
  int hdes = lvmsg->req.ht;
   
  int myWidthPx = wdes / patchArrayDimX;
  int myHeightPx = hdes / patchArrayDimY;
  int sx=thisIndex.x*myWidthPx;
  int sy=thisIndex.y*myHeightPx; 

  // set the output pixel values for rectangle
  // Each component is a char which can have 256 possible values
  unsigned char *intensity= new unsigned char[3*myWidthPx*myHeightPx];
  for(int i=0; i<myHeightPx; ++i)
    for(int j=0; j<myWidthPx; ++j)
      color_pixel(intensity,myWidthPx,myHeightPx,j,i,0,0,0);	// black background

  for (int i=0; i < particles.length(); i++ ) {
    int xpos = (int)((particles[i].x /(BigReal) (patchSize*patchArrayDimX)) * wdes) - sx;
    int ypos = (int)((particles[i].y /(BigReal) (patchSize*patchArrayDimY)) * hdes) - sy;

    Color c(particles[i].id);
    color_pixel(intensity,myWidthPx,myHeightPx,xpos+1,ypos,c.R,c.B,c.G);
    color_pixel(intensity,myWidthPx,myHeightPx,xpos-1,ypos,c.R,c.B,c.G);
    color_pixel(intensity,myWidthPx,myHeightPx,xpos,ypos+1,c.R,c.B,c.G);
    color_pixel(intensity,myWidthPx,myHeightPx,xpos,ypos-1,c.R,c.B,c.G);
    color_pixel(intensity,myWidthPx,myHeightPx,xpos,ypos,c.R,c.B,c.G);
  }
        
  for(int i=0; i<myHeightPx; ++i)
    for(int j=0; j<myWidthPx; ++j) {
      // Draw red lines
      if(i==0 || j==0) {
	color_pixel(intensity,myWidthPx,myHeightPx,j,i,128,0,0);
      }
    }
        
  liveVizDeposit(lvmsg, sx,sy, myWidthPx,myHeightPx, intensity, this, max_image_data);
  delete[] intensity;
}
#endif

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

