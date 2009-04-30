/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/** \file Main.C
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
#include "sqrtTable.h"
#include "ConfigList.h"
#include "molfile_plugin.h"
#include "PluginIOMgr.h"
#ifdef USE_SECTION_MULTICAST
  #include "ckmulticast.h"
#endif

/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Patch patchArray;
/* readonly */ CProxy_Compute computeArray;
///* readonly */ CProxy_GridCompute gridComputeArray;
///* readonly */ CProxy_PMECompositor PMECompArray;
///* readonly */ CProxy_PMEDecompositor PMEDecompArray;
/* readonly */ CkGroupID mCastGrpID;

/* readonly */ vdwParams *vdwTable;
/* readonly */ sqrtTable *rootTable;

/* readonly */ int numParts;
/* readonly */ int patchArrayDimX;	// Number of Chares in X
/* readonly */ int patchArrayDimY;	// Number of Chares in Y
/* readonly */ int patchArrayDimZ;	// Number of Chares in Z
/* readonly */ int patchSize;
/* readonly */ int ptpCutOff;
/* readonly */ int patchMargin;
/* readonly */ int patchOriginX;
/* readonly */ int patchOriginY;
/* readonly */ int patchOriginZ;
/* readonly */ int compArrayLenX;
/* readonly */ int compArrayLenY;
/* readonly */ int compArrayLenZ;
/* readonly */ int pmeGridDimX;
/* readonly */ int pmeGridDimY;
/* readonly */ int pmeGridDimZ;
/* readonly */ BigReal pmeCellSize;
/* readonly */ int pmeCutOff;
/* readonly */ int migrateStepCount;
/* readonly */ int finalStepCount; 
/* readonly */ BigReal stepTime; 
/* readonly */ BigReal timeDelta;
/* readonly */ bool usePairLists;

/* readonly */ double A = 1.60694452*pow(10, -134);			// Force Calculation parameter 1
/* readonly */ double B = 1.031093844*pow(10, -77);			// Force Calculation parameter 2

// Entry point of Charm++ application
Main::Main(CkArgMsg* msg) {
  stepTime = CmiWallTimer();
  CkPrintf("\nLENNARD JONES MOLECULAR DYNAMICS RUNNING ...\n");
  numParts = 0;

  usePairLists = USE_PAIRLISTS;
  patchArrayDimX = PATCHARRAY_DIM_X;
  patchArrayDimY = PATCHARRAY_DIM_Y;
  patchArrayDimZ = PATCHARRAY_DIM_Z;
  ptpCutOff = PTP_CUT_OFF;
  patchMargin = PATCH_MARGIN;
  patchSize = PATCH_SIZE;
  patchOriginX = PATCH_ORIGIN_X;
  patchOriginY = PATCH_ORIGIN_Y;
  patchOriginZ = PATCH_ORIGIN_Z;
  compArrayLenX = COMPARRAY_LEN_X;
  compArrayLenY = COMPARRAY_LEN_Y;
  compArrayLenZ = COMPARRAY_LEN_Z;
  pmeGridDimX = PMEGRID_DIM_X;
  pmeGridDimY = PMEGRID_DIM_Y;
  pmeGridDimZ = PMEGRID_DIM_Z;
  pmeCellSize = PME_CELL_SIZE;
  pmeCutOff = PME_CUT_OFF;
  migrateStepCount = MIGRATE_STEPCOUNT;
  finalStepCount = DEFAULT_FINALSTEPCOUNT;
  structureFilename = STRUCTURE_FILENAME;
  //paramsFileName = PARAMS_FILENAME;
  timeDelta = DEFAULT_DELTA;

  simParams = new SimParameters();  
  simParams->paraTypeCharmmOn = false;
  simParams->paraTypeXplorOn = true;
  simParams->cosAngles = false;  

  mainProxy = thisProxy;
  phase = 0;


  int numPes = CkNumPes();
  int currPe = -1, pe;

  //get square root interpolation table
  rootTable = fillTable(10000, ((BigReal)ptpCutOff*ptpCutOff)*pow(10,-18)/10000);

  //read config file
  if (msg->argc > 1) {
    readConfigFile(msg->argv[1]);
    if (simParams->paraTypeCharmmOn || simParams->paraTypeXplorOn)
      readParameterFile(sl_parameters, simParams);
  }


  //reading data
  FileDataMsg *fdmsg = readParticleData();

  // initializing the 3D patch array

  // patchArray = CProxy_Patch::ckNew(fdmsg, patchArrayDimX, patchArrayDimY, patchArrayDimZ);

  patchArray = CProxy_Patch::ckNew();

  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++) {
	FileDataMsg *fdmsgcopy;
	if (x == patchArrayDimX -1 && y == patchArrayDimY -1 && z == patchArrayDimZ -1)
	  fdmsgcopy = fdmsg;
	else {
	  fdmsgcopy = new (numParts, numParts, numParts*3, numParts) FileDataMsg;
	  memcpy(fdmsgcopy->mass, fdmsg->mass, numParts*sizeof(BigReal));
	  memcpy(fdmsgcopy->charge, fdmsg->charge, numParts*sizeof(BigReal));
	  memcpy(fdmsgcopy->coords, fdmsg->coords, numParts*3*sizeof(BigReal));
	  memcpy(fdmsgcopy->vdw_type, fdmsg->vdw_type, numParts*sizeof(int));
	  fdmsgcopy->length = fdmsg->length;
	}
	pe = (++currPe) % numPes;
	patchArray(x, y, z).insert(fdmsgcopy, pe);
      }
  patchArray.doneInserting();
  CkPrintf("%d PATCHES CREATED\n", patchArrayDimX * patchArrayDimY * patchArrayDimZ);

#ifdef USE_SECTION_MULTICAST
  // initializing the CkMulticastMgr
  mCastGrpID = CProxy_CkMulticastMgr::ckNew();
#endif

  // initializing the 6D compute array
  computeArray = CProxy_Compute::ckNew();
  for (int x=0; x<patchArrayDimX; x++)
    for (int y=0; y<patchArrayDimY; y++)
      for (int z=0; z<patchArrayDimZ; z++)
	patchArray(x, y, z).createComputes();


  delete msg;
  // initializing GridCompute array
//  gridComputeArray = CProxy_GridCompute::ckNew(patchArrayDimX, patchArrayDimY, patchArrayDimZ);

  // initializing PME comp / decomp arrays
//  CkArrayOptions opts(patchArrayDimX, patchArrayDimY, patchArrayDimZ);
//  opts.bindTo(gridComputeArray);
//  PMECompArray = CProxy_PMECompositor::ckNew(opts);
//  PMEDecompArray = CProxy_PMEDecompositor::ckNew(pmeGridDimX / compArrayLenX, pmeGridDimY / compArrayLenY, pmeGridLenZ / compArrayLenZ);
}

// Constructor for chare object migration
Main::Main(CkMigrateMessage* msg) { }

void Main::allDone() {
  CkPrintf("SIMULATION COMPLETE. \n\n");  CkExit();
}

void Main::startUpDone() {
  switch(phase) {
    case 0:
      computeArray.doneInserting();
      CkPrintf("%d COMPUTES CREATED\n", (NUM_NEIGHBORS/2+1) * patchArrayDimX * patchArrayDimY * patchArrayDimZ);
#ifdef USE_SECTION_MULTICAST
      phase++;
      patchArray.createSection();
      break;

    case 1:
      CkPrintf("MULTICAST SECTIONS CREATED\n");
#endif

#ifdef RUN_LIVEVIZ
      // setup liveviz
      CkCallback c(CkIndex_Patch::requestNextFrame(0), patchArray);
      liveVizConfig cfg(liveVizConfig::pix_color,true);
      liveVizInit(cfg,patchArray,c);
#endif

      patchArray.start();
      break;
  }
}

// retrieve necessary information out of a .namd file
void Main::readConfigFile(const char* filename){
  ConfigList *cfg;
  string sON = "on";
  cfg = new ConfigList(filename);
  const StringList* sl_structure = cfg->find("structure");
  if (sl_structure != NULL)
    structureFilename = sl_structure->data;
  
  const StringList* sl_numsteps = cfg->find("numsteps");
  if (sl_numsteps != NULL)
    finalStepCount = atoi(sl_numsteps->data);
  
  const StringList* sl_stepspercycle = cfg->find("stepspercycle");
  if (sl_stepspercycle != NULL)
    migrateStepCount = atoi(sl_stepspercycle->data);
   
  const StringList* sl_cutoff = cfg->find("cutoff");
  if (sl_cutoff != NULL)
    ptpCutOff = atoi(sl_cutoff->data);
  
  const StringList* sl_margin = cfg->find("margin");
  if (sl_margin != NULL)
    patchMargin = atoi(sl_margin->data);
  
  patchSize = patchMargin + ptpCutOff;

  const StringList* sl_timestep = cfg->find("timestep");
  if (sl_timestep != NULL)
    timeDelta = atof(sl_timestep->data);

  sl_parameters = cfg->find("parameters");

  const StringList* sl_paraTypeXplor = cfg->find("paraTypeXplor");
  if (sl_paraTypeXplor != NULL && sON.compare(sl_paraTypeXplor->data) == 0){
    simParams->paraTypeXplorOn = true;
    simParams->paraTypeCharmmOn = false;
  }

  const StringList* sl_paraTypeCharmm = cfg->find("paraTypeCharmm");
  if (sl_paraTypeCharmm != NULL && sON.compare(sl_paraTypeCharmm->data) == 0){
    simParams->paraTypeCharmmOn = true;
    simParams->paraTypeXplorOn = false;
  }
}

// retrieves particle data out of a .js file
FileDataMsg* Main::readParticleData() {
  int numAtoms;
  PluginIOMgr *pIOMgr = new PluginIOMgr();
  molfile_plugin_t *pIOHdl = pIOMgr->getPlugin();

  if (pIOHdl == NULL) {
    CkPrintf("ERROR: Failed to match requested plugin type\n");
  }


  void *pIOFileHdl = pIOHdl->open_file_read(structureFilename,
                                              pIOHdl->name, &numAtoms);
  numParts = numAtoms;
  if(pIOFileHdl ==  NULL) {
    CkPrintf("ERROR: Opening structure file failed!\n");
    CkExit();
  }
  CkPrintf("opened file to read, numatoms = %d \n", numAtoms);

  //get charge and mass
  int optflags = MOLFILE_BADOPTIONS;
  molfile_atom_t *atomarray = (molfile_atom_t *) malloc(numAtoms*sizeof(molfile_atom_t));
  memset(atomarray, 0, numAtoms*sizeof(molfile_atom_t));
  
  int rc = pIOHdl->read_structure(pIOFileHdl, &optflags, atomarray);
  
  if (rc != MOLFILE_SUCCESS && rc != MOLFILE_NOSTRUCTUREDATA) {
    free(atomarray);
    CkPrintf("ERROR: plugin failed reading structure data\n");
    CkExit();
  }
  if(optflags == MOLFILE_BADOPTIONS) {
    free(atomarray);
    CkPrintf("ERROR: plugin didn't initialize optional data flags\n");
    CkExit();
  }
  
  FileDataMsg *fdmsg = new (numAtoms, numAtoms, numAtoms, numAtoms) FileDataMsg;
  Atom fakeAtom;
  for (int i = 0; i < numAtoms; i++){
    fdmsg->mass[i] = atomarray[i].mass; 
    fdmsg->charge[i] = atomarray[i].charge; //charge may be being received with wrong units!!! previously had factor of 10 / 4
    //CkPrintf("atom name = %s atom type = %s\n",atomarray[i].name, atomarray[i].type);
    fileParams->assign_vdw_index(atomarray[i].type, &fakeAtom);
    fdmsg->vdw_type[i] = fakeAtom.vdw_type;
  }

  fdmsg->length = numAtoms;
  //get coordinates
  molfile_timestep_t ts;
  float *atomcoords;
  memset(&ts, 0, sizeof(molfile_timestep_t));
  
  atomcoords = (float *) malloc(3*numAtoms*sizeof(float));
  memset(atomcoords, 0, 3*numAtoms*sizeof(float));
  ts.coords = atomcoords;
  if (pIOHdl->read_next_timestep(pIOFileHdl, numAtoms, &ts)) {
    free(atomcoords);
    pIOHdl->close_file_read(pIOFileHdl);
    CkPrintf("ERROR: failed reading atom coordinates");
    CkExit();
  }
  float minX, maxX, minY, maxY, minZ, maxZ;
  minX = maxX = minY = maxY = minZ = maxZ = 0;
  for (int i = 0; i < numAtoms; i++){
    if (atomcoords[0] > maxX)
      maxX = atomcoords[0];
    if (atomcoords[0] < minX)
      minX = atomcoords[0];
    if (atomcoords[1] > maxY)
      maxY = atomcoords[1];
    if (atomcoords[1] < minY)
      minY = atomcoords[1];
    if (atomcoords[2] > maxZ)
      maxZ = atomcoords[2];
    if (atomcoords[2] < minZ)
      minZ = atomcoords[2];
    fdmsg->coords[i].x = atomcoords[0];
    fdmsg->coords[i].y = atomcoords[1];
    fdmsg->coords[i].z = atomcoords[2];
    atomcoords += 3;
  }
  // determine appropriate patch dimensions
  patchArrayDimX = ceil((maxX - minX) / patchSize);
  patchArrayDimY = ceil((maxY - minY) / patchSize);
  patchArrayDimZ = ceil((maxZ - minZ) / patchSize);
  patchOriginX = (int)minX;
  patchOriginY = (int)minY;
  patchOriginZ = (int)minZ;
  CkPrintf("origin is [%d][%d][%d]\n", patchOriginX, patchOriginY, patchOriginZ);
  return fdmsg;

  // do we need to read the velocities?
}

// reads parameter file
void Main::readParameterFile(const StringList* sl_params, SimParameters* sParams){
  int numParams;
  vdwPars *vdwp;
  Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
  Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;
  Real rA, rB, rA14, rB14;
  CkPrintf("Reading parameter files\n");
  fileParams = new Parameters(sParams, (StringList*)sl_params);
  numParams = fileParams->get_num_vdw_params();
  CkPrintf("Assigned %d vdw_params\n", numParams);
  vdwTable = new (numParams*numParams) vdwParams;
  vdwTable->numParams = numParams;
  for (int i = 0; i < numParams; i++){
    for (int j = 0; j < numParams; j++){
      vdwp = &vdwTable->params[i*numParams +j];
      if (fileParams->get_vdw_pair_params((Index)i, (Index)j, &rA, &rB, &rA14, &rB14) == 1){
	CkPrintf("retreived params for (%d, %d)\n", i,j);
	vdwp->A = (BigReal)rA;
	vdwp->B = (BigReal)rB;
	vdwp->A14 = (BigReal)rA14;
	vdwp->B14 = (BigReal)rB14;
      }
      else {
	fileParams->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
                                       &epsilon_i14, (Index)i);
	fileParams->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14,
                                       &epsilon_j14, (Index)j);
	BigReal sigma_ij = sqrt(sigma_i*sigma_j);
	BigReal sigma_ij14 = sqrt(sigma_i14*sigma_j14);
	BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);
	BigReal epsilon_ij14 = sqrt(epsilon_i14*epsilon_j14);

	sigma_ij *= sigma_ij*sigma_ij;
	sigma_ij *= sigma_ij;
	sigma_ij14 *= sigma_ij14*sigma_ij14;
	sigma_ij14 *= sigma_ij14;

	//  Calculate LJ constants A & B
	vdwp->B = 4.0 * sigma_ij * epsilon_ij;
	vdwp->A = vdwp->B * sigma_ij;
	vdwp->B14 = 4.0 * sigma_ij14 * epsilon_ij14;
        vdwp->A14 = vdwp->B14 * sigma_ij14;
      }
    }
  }
}
#include "mol3d.def.h"
