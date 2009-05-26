#ifndef __TESTINSERT_H__
#define __TESTINSERT_H__

class dataMsg : public CMessage_dataMsg {
  public:
    int *stuff;
    int length;
};

class Main : public CBase_Main {
  private:
    int dimX, dimY, dimZ;
    int counter;

  public:
    Main(CkArgMsg* msg);
    Main(CkMigrateMessage* msg);
    void checkIn();
};

class Elements : public CBase_Elements {
  private:
    int *myData;

  public:
    Elements(CkMigrateMessage* msg);
    Elements(dataMsg *msg);
};
#endif
