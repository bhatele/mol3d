/** \file Reader.C
 *  Authors: Abhinav S Bhatele
 *  Date Created: October 14th, 2010
 */

#include "SimParameters.h"
#include "Parameters.h"

#include "defs.h"
#include "mol3d.decl.h"
#include "Reader.h"

#include "sqrtTable.h"
#include "ConfigList.h"
#include "molfile_plugin.h"
#include "PluginIOMgr.h"

PUPbytes(SimParameters)

/* readonly */ vdwParams *vdwTable;
/* readonly */ sqrtTable *rootTable;

/* readonly */ int ptpCutOff;
/* readonly */ int patchMargin;
/* readonly */ int patchOriginX;
/* readonly */ int patchOriginY;
/* readonly */ int patchOriginZ;
/* readonly */ int migrateStepCount;
/* readonly */ int finalStepCount;
/* readonly */ int firstLdbStep;
/* readonly */ int ldbPeriod;
/* readonly */ int ftPeriod;
/* readonly */ BigReal timeDelta;
/* readonly */ bool usePairLists;

extern /* readonly */ CProxy_Main mainProxy;
extern /* readonly */ CProxy_Patch patchArray;
extern /* readonly */ CProxy_Compute computeArray;
extern /* readonly */ CkGroupID mCastGrpID;

extern /* readonly */ int numParts;
extern /* readonly */ int patchArrayDimX;       // Number of Chares in X
extern /* readonly */ int patchArrayDimY;       // Number of Chares in Y
extern /* readonly */ int patchArrayDimZ;       // Number of Chares in Z
extern /* readonly */ int patchSizeX;
extern /* readonly */ int patchSizeY;
extern /* readonly */ int patchSizeZ;
extern /* readonly */ bool twoAwayX;
extern /* readonly */ bool twoAwayY;
extern /* readonly */ bool twoAwayZ;
extern /* readonly */ int numNbrs;
extern /* readonly */ int nbrsX;
extern /* readonly */ int nbrsY;
extern /* readonly */ int nbrsZ;

class Reader : public CBase_Reader {
  private:
    const char* structureFilename;
    const StringList* sl_parameters;
    SimParameters *simParams;	// not the full simparameters class from namd
    Parameters *fileParams;

  public:

    Reader(CkArgMsg *msg) {
      usePairLists = USE_PAIRLISTS;

      ptpCutOff = PTP_CUT_OFF;
      patchMargin = PATCH_MARGIN; 
      patchOriginX = PATCH_ORIGIN_X;
      patchOriginY = PATCH_ORIGIN_Y;
      patchOriginZ = PATCH_ORIGIN_Z;
  
      migrateStepCount = MIGRATE_STEPCOUNT;
      finalStepCount = DEFAULT_FINALSTEPCOUNT;
      firstLdbStep = DEFAULT_FIRST_LDB;
      ldbPeriod = DEFAULT_LDB_PERIOD;
      ftPeriod = DEFAULT_FT_PERIOD;
      timeDelta = DEFAULT_DELTA;

      twoAwayX = false;
      twoAwayY = false;
      twoAwayZ = false;

      simParams = new SimParameters();
      simParams->paraTypeCharmmOn = false;
      simParams->paraTypeXplorOn = true;
      simParams->cosAngles = false;

      // read config file
      readConfigFile(msg->argv[1]);
      if (simParams->paraTypeCharmmOn || simParams->paraTypeXplorOn)
	readParameterFile(sl_parameters, simParams);

      // get square root interpolation table
      rootTable = fillTable(10000, ((BigReal)ptpCutOff*ptpCutOff)*pow(10.0,-20)/10000);

      // initializing the 3D patch array
      patchArray = CProxy_Patch::ckNew();
      computeArray = CProxy_Compute::ckNew();

      //reading data
      FileDataMsg *fdmsg = readParticleData();

      FileDataMsg **fdmsgs;
      fdmsgs = new FileDataMsg*[patchArrayDimZ*patchArrayDimY*patchArrayDimX];
      CkVec<int> partsToSend;
      int numSend;
      int numPes = CkNumPes(), currPe = -1, pe;

      for (int x=0; x<patchArrayDimX; x++)
	for (int y=0; y<patchArrayDimY; y++)
	  for (int z=0; z<patchArrayDimZ; z++) {
	    partsToSend.removeAll();
	    for (int i = 0; i < numParts; i++) {
	      if (((int)((fdmsg->coords[i].x - patchOriginX) / patchSizeX)) == x && ((int)((fdmsg->coords[i].y-patchOriginY) / patchSizeY)) == y
		 && ((int)((fdmsg->coords[i].z-patchOriginZ) / patchSizeZ)) == z)
		partsToSend.push_back(i);
	    }
	    numSend = partsToSend.size();
	    fdmsgs[x*patchArrayDimZ*patchArrayDimY+y*patchArrayDimZ+z] = new (numSend, numSend, numSend, numSend) FileDataMsg;
	    FileDataMsg *fdmsgcopy = fdmsgs[x*patchArrayDimZ*patchArrayDimY+y*patchArrayDimZ+z];
	    //fdmsgcopy = new (numSend, numSend, numSend, numSend) FileDataMsg;
	    for (int i = 0; i < numSend; i++){
	      fdmsgcopy->mass[i] = fdmsg->mass[partsToSend[i]];
	      fdmsgcopy->charge[i] = fdmsg->charge[partsToSend[i]];
	      fdmsgcopy->coords[i] = fdmsg->coords[partsToSend[i]];
	      fdmsgcopy->vdw_type[i] = fdmsg->vdw_type[partsToSend[i]];
	    }
	    fdmsgcopy->length = numSend;
	    pe = (++currPe) % numPes;
	    patchArray(x, y, z).insert(fdmsgcopy, pe);
	  }
      patchArray.doneInserting();
    }

    Reader(CkMigrateMessage* msg): CBase_Reader(msg) { }

    void pup(PUP::er &p) {
      Chare::pup(p);
      if (p.isUnpacking())  simParams = new SimParameters;
      p|*simParams;
    }

    /* retrieve necessary information out of a .namd file */
    void readConfigFile(const char* filename) {
      ConfigList *cfg;
      string sON = "on";
      cfg = new ConfigList(filename);

      const StringList* sl_structure = cfg->find("structure");
      if (sl_structure != NULL)
	structureFilename = sl_structure->data;
     
      const StringList* sl_numsteps = cfg->find("numsteps");
      if (sl_numsteps != NULL)
	finalStepCount = atoi(sl_numsteps->data);
     
      const StringList* sl_twoawayx = cfg->find("twoAwayX");
      if (sl_twoawayx != NULL)
	twoAwayX = sl_twoawayx->data[0] == 'y';
      if (twoAwayX)
	CkPrintf("performing 2-away X decomposition\n");
     
      const StringList* sl_twoawayy = cfg->find("twoAwayY");
      if (sl_twoawayy != NULL)
	twoAwayY = sl_twoawayy->data[0] == 'y';
      if (twoAwayY)
	CkPrintf("performing 2-away Y decomposition\n");

      const StringList* sl_twoawayz = cfg->find("twoAwayZ");
      if (sl_twoawayz != NULL)
	twoAwayZ = sl_twoawayz->data[0] == 'y';
      if (twoAwayZ)
	CkPrintf("performing 2-away Z decomposition\n");

      const StringList* sl_stepspercycle = cfg->find("stepspercycle");
      if (sl_stepspercycle != NULL)
	migrateStepCount = atoi(sl_stepspercycle->data);

      const StringList* sl_cutoff = cfg->find("cutoff");
      if (sl_cutoff != NULL)
	ptpCutOff = atoi(sl_cutoff->data);

      const StringList* sl_margin = cfg->find("margin");
      if (sl_margin != NULL)
	patchMargin = atoi(sl_margin->data);

      if (twoAwayX) {
	patchSizeX = (2*patchMargin + ptpCutOff)/2;
	nbrsX = 5;
      } else
	patchSizeX = 2*patchMargin + ptpCutOff;

      if (twoAwayY) {
	patchSizeY = (2*patchMargin + ptpCutOff)/2;
	nbrsY = 5;
      } else
	patchSizeY = 2*patchMargin + ptpCutOff;

      if (twoAwayZ) {
	patchSizeZ = (2*patchMargin + ptpCutOff)/2;
	nbrsZ = 5;
      } else
	patchSizeZ = 2*patchMargin + ptpCutOff;

      numNbrs = nbrsX * nbrsY * nbrsZ;

      const StringList* sl_timestep = cfg->find("timestep");
      if (sl_timestep != NULL)
	timeDelta = atof(sl_timestep->data);

      const StringList* sl_firstLdbStep = cfg->find("firstLdbStep");
      if (sl_firstLdbStep != NULL)
	firstLdbStep = atof(sl_firstLdbStep->data);

      const StringList* sl_ldbPeriod = cfg->find("ldbPeriod");
      if (sl_ldbPeriod != NULL)
	ldbPeriod = atof(sl_ldbPeriod->data);

      const StringList* sl_ftPeriod = cfg->find("ftPeriod");
      if (sl_ftPeriod != NULL)
	ftPeriod = atof(sl_ftPeriod->data);

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

    /* retrieves particle data out of a .js file */
    FileDataMsg* readParticleData() {
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
	fdmsg->charge[i] = atomarray[i].charge; //charge may be being received with
	// wrong units!!! previously had factor of 10 / 4
	fileParams->assign_vdw_index(atomarray[i].type, &fakeAtom);
	fdmsg->vdw_type[i] = fakeAtom.vdw_type;
      }

      fdmsg->length = numAtoms;
      // get coordinates
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

      for (int i = 0; i < numAtoms; i++) {
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
      patchArrayDimX = floor((maxX - minX) / patchSizeX);
      patchSizeX = ceil((maxX - minX) / patchArrayDimX);
      patchArrayDimY = floor((maxY - minY) / patchSizeY);
      patchSizeY = ceil((maxY - minY) / patchArrayDimY);
      patchArrayDimZ = floor((maxZ - minZ) / patchSizeZ);
      patchSizeZ = ceil((maxZ - minZ) / patchArrayDimZ);
      patchOriginX = (int)minX;
      patchOriginY = (int)minY;
      patchOriginZ = (int)minZ;

      // do we need to read the velocities?
      return fdmsg;
    }

    /* reads parameter file */
    void readParameterFile(const StringList* sl_params, SimParameters* sParams) {
      int numParams;
      vdwPars *vdwp;
      Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
      Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;
      Real rA, rB, rA14, rB14;
      CkPrintf("Reading parameter files\n");
      fileParams = new Parameters(sParams, (StringList*)sl_params);
      numParams = fileParams->get_num_vdw_params();
      // CkPrintf("Assigned %d vdw_params\n", numParams);
      vdwTable = new (numParams*numParams) vdwParams;
      vdwTable->numParams = numParams;
      for (int i = 0; i < numParams; i++) {
	for (int j = 0; j < numParams; j++) {
	  vdwp = &vdwTable->params[i*numParams +j];
	  if(fileParams->get_vdw_pair_params((Index)i, (Index)j, &rA, &rB, &rA14, &rB14) == 1) {
	    // CkPrintf("retreived params for (%d, %d)\n", i,j);
	    vdwp->A = (BigReal)rA;
	    vdwp->B = (BigReal)rB;
	    vdwp->A14 = (BigReal)rA14;
	    vdwp->B14 = (BigReal)rB14;
	  } else {
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
};

#include "reader.def.h"
