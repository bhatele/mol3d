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
    CProxy_Reader rd;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);

    void pup(PUP::er &p) {
      Chare::pup(p);
      p|phase;
    }

    void allDone();
    void lbBarrier();
    void ftBarrier();
    void startUpDone();
};
#endif
