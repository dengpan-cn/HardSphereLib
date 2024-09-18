#ifndef __STRUCTSIM__
#define __STRUCTSIM__

#include "VecMathEventDriven_ndim.h"

typedef struct Box {
    int dim;
    uptriMat boxH, invBoxH;
    double volume;
    
    doubleVector boxEdge[DIM];
    
    doubleVector cornerLo;
    doubleVector cornerHi;
    
    bool isShapeFixed;
} Box;
typedef struct Particle {
    int nAtom, nAtomType;
    
    doubleVector* pos;
    doubleVector* upos;
    doubleVector* veloc;
    intVector* img;
    double* timeStamp;
    
    // the diameter in code is diameterScale[iatom] * meandiameter.
    double *diameterScale, meanDiameter;
    
    //===just for future feature===
    int* type;
    double* mass;
    double* massPerType;
    //===just for future feature===
    
    int *id2tag, *tag2id;  // tag is the label, id is the index;
    int sortFlag;          // increased when the index is reordered.

    // The library MUST NOT rewrite isSortForbidden.
    bool isSortForbidden;  // default: false. True: Never doing sorting and (label == index).
    
    bool isSizeFixed;
    bool isSync;
} Particle;

#define __minSkinSet__ 5E-2
#define __maxSkinSet__ 1.0
#define AveNrebuild 10
#define colPerAtomBetweenRebuild 2.0
typedef struct idPosRadius {
    int id, parent;
    doubleVector pos;
    double radius;
    doubleVector veloc;
} idPosRadius;
typedef struct NebrList {
    double maxDiameterScale, minDiameterScale;
    double skinSet, rskin, maxRskinSet;
    
    doubleVector binLen;
    intVector nbin;
    int totBin, allocBin;
    
    doubleVector* xyzHold;
    
    long long int accColRebuild;   // record collision used for tune parameter;
    long long int colRebuild;      // # of collisions;
    //bool nptisoFlag;  // flag for nptiso
    //====bin list====
    idPosRadius *binHead, *binList;
    int maxBinHead;
    //====Nebr list====
    int2* list;  // first: nebr; second: colPair;
    int* nNebr;
    int maxNebrPerAtom;
    
    double* colPairTable;
    int allocColTable;
    
    double* tnebr;  // predicted time of tnebr
    //===adjacent bin===
    intVector* deltaAdjBin;
    int nAdjBin;
    //====
    
    int *binHead4sort, *binList4sort;
    int totBin4sort, allocBin4sort;
    int* oid2nid;
    void* buffer;
    
    long int nBuild;
    bool isValid, doSort;
    bool compelInit;
} NebrList;

#define JammingZmin 1E12
//#define NoseHoverChain 20
typedef enum ReturnType {// the state is synchronized except TrivialReturn.
    TrivialReturn = 1 << 0,       // trivial Return;
    DensityReturn = 1 << 1,       // reaching target density;
    PressReturn = 1 << 2,         // reaching target pressure;
    JamReturn = 1 << 3,           // Jammed before reaching target density or Z;
    ThermoReturn = 1 << 4,        // return due to output cmd;
    ErrorReturn = 1 << 5,         // return due to error;
    HaltReturn = 1 << 6,          // Halt running [currently trigered by SWAP];
} ReturnType;

typedef struct Update {
    double massUnits, energyUnits, distanceUnits;
    double timeUnits, forceUnits, velocityUnits;
    double pressureUnits, volumeUnits;  // units
    
    double Epair, volFrac, Tint, Z, dof;
    uptriMat Ztensor, Kintensor;
    uptriMat sColVirialTensor;
    double accTimePeriod;//operated by calcPressure().
    int outputPeriod, nextOutputStep;
    bool Edone, Zdone, Tdone;
    
    NebrList nebrList;
    MinEventHeap eventList;
    
    double kTt;    // temperature 1.0;
    double rrateSet; // \dot{X}/X = rrate; X is diameter or radius;
    
    double currentStamp;//operated by syncAll().
    double runtimeReal;
    
    int ZMperiod, ZCperiod, Ttperiod;
    int nextZMstep, nextZCstep, nextTtstep;
    
    int colCnt, errCnt;
    int stepCol, stepErr;
    // when colCnt == nAtom; do 0 -> colCnt and stepCol++;
    
    ReturnType rtype;
    Toolkit toolkit;
} Update;


#endif
