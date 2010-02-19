/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file Main.h
 *  Author: Abhinav S Bhatele
 *  Date Created: August 11th, 2008
 */
#include "SimParameters.h"
#include "Parameters.h"

#ifndef __MAIN_H__
#define __MAIN_H__


/** \class Main
 *
 */
class Main : public CBase_Main {
  private:
    int phase;
    const char* structureFilename;
    const StringList* sl_parameters;
    SimParameters *simParams; //not actually the full simparameters class from namd
    Parameters *fileParams;    

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    void readConfigFile(const char* filename);
    void readParameterFile(const StringList* sl_params, SimParameters* sParams);
    FileDataMsg* readParticleData();
    void allDone();
    void lbBarrier();
    void ftBarrier();
    void startUpDone();
};
#endif
