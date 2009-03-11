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


#ifndef __MAIN_H__
#define __MAIN_H__

/** \class Main
 *
 */
class Main : public CBase_Main {
  private:
    int phase;
    const char* structureFilename;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    void readConfigFile(const char* filename);
    FileDataMsg* readParticleData();
    void allDone();
    void startUpDone();
};
#endif
