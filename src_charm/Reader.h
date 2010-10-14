/** \file Reader.C
 *  Authors: Abhinav S Bhatele
 *  Date Created: October 14th, 2010
 */

#ifndef __READER_H__
#define __READER_H__


class vdwParams : public CMessage_vdwParams {
  public:
    vdwPars *params;
    int numParams;
};


class sqrtPars {
  public:
    BigReal a; //constant
    BigReal b;
    BigReal c;
    BigReal d; //x^3

    void pup(PUP::er &p) {
      p | a; p | b;
      p | c;  p | d;
    }
};


class sqrtTable : public CMessage_sqrtTable {
  public:
    sqrtPars* pars;
    int length;
    BigReal delta;
};


class loc{
  public:
    BigReal x;
    BigReal y;
    BigReal z;

    void pup(PUP::er &p){
      p|x; p|y; p|z;
    }
};


struct force {
  BigReal x;
  BigReal y;
  BigReal z;
};


class partData{
  public:
    loc coord;
    BigReal charge;
    int vdwIndex;
    
    void pup(PUP::er &p){
      p|coord; p|charge; p|vdwIndex;
    }
};


class FileDataMsg : public CMessage_FileDataMsg {
  public:
    BigReal* charge;
    BigReal* mass;
    loc* coords; //encoded as x1 y1 z1 x2 y2 z2...
    int* vdw_type;
    int length;
};

#endif
