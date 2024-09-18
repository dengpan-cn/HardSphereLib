#ifndef __SUBFUNC__
#define __SUBFUNC__

#include "StructEventDriven_ndim.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*FuncPtrToolkitWriteConf)(FILE* fout, Box* box, Particle* particle, Update* update);

void readCmdLineVar(Variable* var, int argc, char const* argv[]);

//====read and write configurations======
void readConf(Box* box, Particle* particle, Update* update, Variable* var);
void writeConf(Box* box, Particle* particle, Update* update, Variable* var);
int emergWriteConfHalt(Box* box, Particle* particle, Update* update, char* fname);//write configurations then halt.

// turn on/off memory reordering.
// Memory reordering (sorting) is performed at the fisrt time of builing neighbour list and every 200 times of neighbour list rebuilding.
// The sorting is performed when building neighbour list.
int turnOffMemSort(Box *box, Particle *particle, Update *update);// Abort if failed!
int turnOnMemSort(Box *box, Particle *particle, Update *update);// Abort if failed!
int checkMemSortFlag(Box *box, Particle *particle, Update *update);// Abort if failed!

//=======box deformation========
void setBoxPara(Box* box);
void calcBoxLenAngle(Box* box, uptriMat para);
void calcDistBoxPlane(doubleVector distPlane, doubleVector boxEdge[DIM]);
void setUnits(Update* update, double distanceUnits);
void reInitSim(Box* box, Particle* particle, Update* update);

void syncAll(Box* box, Particle* particle, Update* update);
void adjustImg(Box* box, Particle* particle);

void calcKinTensor(Box* box, Particle* particle, Update* update);
void calcPressure(Box* box, Particle* particle, Update* update);
void calcPoten(Box* box, Particle* particle, Update* update);

void thermostat(Box* box, Particle* particle, Update* update);
void zeroMomentum(Box* box, Particle* particle, Update* update);
void moveBarycenter(Box* box, Particle* particle, Update* update);

//=======generating velocities================
void genGaussianVeloc(Box* box, Particle* particle, Update* update, double targT);

void instant_inflate(Box* box, Particle* particle, Update* update, double deltaVF);
//=============================================
int checkOverlap(Box* box, Particle* particle, Update* update);

#if (DIM == 2 || DIM == 3)
void writeLammpsTopo(Box* box, Particle* particle, Update* update,  char *fname);
#endif

//============write and read dump File==========================
typedef struct writeDumpFile {
    int revNum;
    FILE* fdump;
} writeDumpFile;
writeDumpFile* getWriteDumpFile(Update* update);
writeDumpFile* addWriteDumpFile(Box* box, Particle* particle, Update* update, Variable* var);
int writeDump(Box* box, Particle* particle, Update* update);
int delWriteDumpFile(Update* update);

typedef struct mmapBinFile {
    int refCnt;
    int fd;               // the descriptor the mapped file
    struct stat binStat;  // the stat of the opened file
    void* dataSection;    // pointer to the mapped memory
    int nStep;            // the number of timestep
    int stepSize;
    int headerSize;
    int revNum;
} mmapBinFile;  // binary dump file
mmapBinFile* openBinFile(char* fname);
int readSimInfo(Box* box, Particle* particle, Update* update, mmapBinFile* binFile);
int readDump(Box* box, Particle* particle, Update* update, mmapBinFile* binFile, int whichStep);
int closeBinFile(mmapBinFile** binFilePtr);

void generateUnwarpPos(Box* box, Particle* particle);

typedef struct SWAP {
    double max_dsig;// sigA - sigB >= max_dsig is rejected.
    double swap_frac;// fraction over collisions.
    
    //====bin list====
    doubleVector binLen;
    intVector nbin;
    int totBin;
    
    idPosRadius *binHead, *binList;
    int maxBinHead;
    intVector* binIdx;
    
    //===adjacent bin===
    intVector* deltaAdjBin;
    int nAdjBin;
    //================
    int stopBySwap;
    doubleVector* xyz0;
    double k_0, selfISF, msd;
    double max_fs, min_msd;
    //====record last state===
    int lastColCnt, lastStepCol;
    
    //===statistics====
    long long int ntrial, naccept;
} SWAP;
SWAP* getSwap(Update* update);
SWAP* addSwap(Box* box, Particle* particle, Update* update, Variable* var);
ReturnType swapMCmove(Box* box, Particle* particle, Update* update);
int delSwap(Update* update);

#ifdef __cplusplus
}
#endif

#endif
