#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

#include "common.h"
#include "Vector.h"
#include "Lattice.h"

#include "MGridforceParams.h"

class SimParameters {
  public:
    bool paraTypeXplorOn;           //  FLAG TRUE-> parametrs are XPLOR format (default)
    bool paraTypeCharmmOn;          //  FLAG TRUE-> parametrs are CHARMM format

    bool cosAngles;
};
#endif
