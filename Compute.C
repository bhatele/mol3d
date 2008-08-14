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
  bufferedX = 0;
  bufferedY = 0;
  bufferedZ = 0;
}

Compute::Compute(CkMigrateMessage *msg) { }
  
// Entry method to receive vector of particles
void Compute::interact(CkVec<Particle> particles, int x, int y, int z) {

  int i;

  // self interaction check
  if( thisIndex.x1 == thisIndex.x2 && thisIndex.y1 == thisIndex.y2 && thisIndex.z1 == thisIndex.z2) {
    interact(particles, particles);
    patchArray(x, y, z).receiveForces(particles);
  } else {
    if(cellCount == 0) {
      bufferedX = x;
      bufferedY = y;
      bufferedZ = z;
      bufferedParticles = particles;
      cellCount++;
    } else if (cellCount == 1) {
      // if both particle sets are received, compute interaction
      cellCount = 0;
      interact(bufferedParticles, particles);
      patchArray(bufferedX, bufferedY, bufferedZ).receiveForces(bufferedParticles);
      patchArray(x, y, z).receiveForces(particles);
    }
  }
}

// Local function to compute all the interactions between pairs of particles in two sets
void Compute::interact(CkVec<Particle> &first, CkVec<Particle> &second){
  int i, j;
  for(i = 0; i < first.length(); i++)
    for(j = 0; j < second.length(); j++)
      interact(first[i], second[j]);
}

// Local function for computing interaction among two particles
// There is an extra test for interaction of identical particles, in which case there is no effect
void Compute::interact(Particle &first, Particle &second){
  float rx, ry, rz, r, fx, fy, fz, f;

  // computing base values
  rx = first.x - second.x;
  ry = first.y - second.y;
  rz = first.z - second.z;
  r = sqrt(rx*rx + ry*ry + rz*rz);

  // We include 0.000001 to ensure that r doesn't tend to zero in the force calculation
  // if(r < 0.000001 || r >= DEFAULT_RADIUS)
  if(r < 0.000001 || r >= patchSize)
    return;

  f = A / pow(r,12) - B / pow(r,6);
  fx = f * rx / r;
  fy = f * ry / r;
  fz = f * rz / r;

  // updating particle properties
  second.fx -= fx;
  second.fy -= fy;
  second.fz -= fz;
  first.fx += fx;
  first.fy += fy;
  first.fz += fz;
}

