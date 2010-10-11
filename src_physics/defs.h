/** \file common.h
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include "pup.h"

typedef double BigReal;

#define AVAGADROS_NUMBER        (6.022141 * pow(10.0,23))
#define COULOMBS_CONSTANT       (8.987551 * pow(10.0,-9))
#define ELECTRON_CHARGE         (1.602176 * pow(10.0,-19))

#define USE_PAIRLISTS		true	// generally faster if true

#define DEFAULT_DELTA		1	// in femtoseconds

#define DEFAULT_FIRST_LDB	101
#define DEFAULT_LDB_PERIOD	500
#define DEFAULT_FT_PERIOD	100000

#define PATCHARRAY_DIM_X	3
#define PATCHARRAY_DIM_Y	3
#define PATCHARRAY_DIM_Z	3
#define PTP_CUT_OFF		13	// Rc in NAMD, cut off for atom to atom interactions
#define PATCH_MARGIN		0 	// constant difference between cut off and patch size
#define PATCH_SIZE_X		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_SIZE_Y		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_SIZE_Z		(PTP_CUT_OFF + PATCH_MARGIN)
#define PATCH_ORIGIN_X		0
#define PATCH_ORIGIN_Y		0
#define PATCH_ORIGIN_Z		0

#define MIGRATE_STEPCOUNT	5
#define DEFAULT_FINALSTEPCOUNT	101
#define MAX_VELOCITY		30.0

#define KAWAY_X			1
#define KAWAY_Y			1
#define KAWAY_Z			1
#define NBRS_X			(2*KAWAY_X+1)
#define NBRS_Y			(2*KAWAY_Y+1)
#define NBRS_Z			(2*KAWAY_Z+1)
#define NUM_NEIGHBORS		(NBRS_X * NBRS_Y * NBRS_Z)

#define WRAP_X(a)		(((a)+patchArrayDimX)%patchArrayDimX)
#define WRAP_Y(a)		(((a)+patchArrayDimY)%patchArrayDimY)
#define WRAP_Z(a)		(((a)+patchArrayDimZ)%patchArrayDimZ)

// Struct for keeping track of vdw parameters
class vdwPars{
  public:
    BigReal A;
    BigReal B;
    BigReal A14;
    BigReal B14;

    vdwPars(){
      A = B = A14 = B14 = 0;
    }
    
    void pup(PUP::er &p) {
      p | A; p | B; p | A14; p | B14;
    }

};

// Class for keeping track of the properties for a particle
class Particle {
  public:
    int vdw_type;

    int id;
    BigReal mass;	// mass of the particle
    BigReal charge;     // charge of particle
    BigReal x;		// position in x axis
    BigReal y;		// position in y axis
    BigReal z;		// position in z axis

    BigReal fx;		// total forces on x axis
    BigReal fy;		// total forces on y axis
    BigReal fz;		// total forces on z axis

    BigReal ax;		// acceleration on x axis
    BigReal ay;		// acceleration on y axis
    BigReal az;		// acceleration on z axis

    BigReal vx;		// velocity on x axis
    BigReal vy;		// velocity on y axis
    BigReal vz;		// velocity on z axis

    // Default constructor
    Particle() {
      fx = fy = fz = 0.0;
    }

    // Function for pupping properties
    void pup(PUP::er &p) {
      p | id; p | mass; p | charge;
      p | x;  p | y;  p | z;
      p | fx; p | fy; p | fz;
      p | ax; p | ay; p | az;
      p | vx; p | vy; p | vz;
      p | vdw_type;
    }
};

class Color {
  public:
    unsigned char R, G, B;

    // Generate a unique color for each index from 0 to total-1
    Color(int index){
      int total = 8;
      if(index % total == 0) {
	R = 255;
	G = 100;
	B = 100;
      } else if(index % total == 1) {
	R = 100;
	G = 255;
	B = 100;
      } else if(index % total == 2) {
	R = 100;
	G = 100;
	B = 255;
      } else if(index % total == 3) {
	R = 100;
	G = 255;
	B = 255;
      } else if(index % total == 4) {
	R = 100;
	G = 255;
	B = 255;
      } else if(index % total == 5) {
	R = 255;
	G = 255;
	B = 100;
      } else if(index % total == 6) {
	R = 255;
	G = 100;
	B = 255;
      } else {
	R = 170;
	G = 170;
	B = 170;
      }
    }	
};

#endif
