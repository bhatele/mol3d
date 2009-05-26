#include "testInsert.decl.h"
#include "testInsert.h"

/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ CProxy_Elements elementsArray;

Main::Main(CkMigrateMessage* msg) { }
Main::Main(CkArgMsg* msg){
  counter = 0;
  CkPrintf("beginning insert test\n");
  dimX = 10;
  dimY = 10;
  dimZ = 10;

  mainProxy = thisProxy;

  if (msg->argc > 1) dimX = atoi(msg->argv[1]);
  if (msg->argc > 2) dimY = atoi(msg->argv[2]);
  if (msg->argc > 3) dimZ = atoi(msg->argv[3]);

  elementsArray = CProxy_Elements::ckNew();  
  
  srand(100);
  int msgSize = 100;

  int numPes = CkNumPes();
  int currPe = -1;
  int pe;

  for (int x = 0; x < dimX; x++){
    for (int y = 0; y < dimY; y++){
      for (int z = 0; z < dimZ; z++){
	dataMsg* dmsg = new (msgSize) dataMsg;
	for (int i = 0; i < msgSize; i++){
	  dmsg->stuff[i] = rand();
	}
	dmsg->length = msgSize;
	pe = (++currPe) % numPes;
	elementsArray(x,y,z).insert(dmsg, pe);
      }
    }
  }
  elementsArray.doneInserting();
  CkPrintf("elements created\n");
}

void Main::checkIn(){
  counter++;
  if (counter == dimX*dimY*dimZ){
    CkPrintf("test successful\n");
    CkExit();
  }
}

Elements::Elements(CkMigrateMessage* msg) {}
Elements::Elements(dataMsg *msg){
  myData = msg->stuff;
  mainProxy.checkIn();
}

#include "testInsert.def.h"
