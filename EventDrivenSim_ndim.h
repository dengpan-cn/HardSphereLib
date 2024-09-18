#ifndef __HARDCOLFORWARD__
#define __HARDCOLFORWARD__

#include "SimEventDrivenSub_ndim.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SimEDMD {
    // performing velocity scaling thermostat whenever building event list.
    double rtimeStepReal;
} SimEDMD;
typedef struct ExpandEDMD {
    // performing velocity scaling thermostat whenever building event list.
    double rrateSet;// \dot{X}/X = rrate; X: diameter or radius;
    
    bool isVolCtrl;
    double targetVF, targetZ;

    double rtimeStepReal;
} ExpandEDMD;

typedef struct NptIsoEDMD {
    // performing velocity scaling thermostat every step.
    double predtStep;
    double rrate;
    bool zeroRate;
    
    double Zstart, Zstop, tauZ, Ztarg, Zrate;
} NptIsoEDMD;

//===Event Driven Molecular Dynamics===
SimEDMD* getSimEdmd(Box* box, Particle* particle, Update* update);
SimEDMD* addSimEdmd(Box* box, Particle* particle, Update* update);
int delSimEdmd(Box* box, Particle* particle, Update* update);
ReturnType colForward_edmd(Box* box, Particle* particle, Update* update, double runTimeReal, int maxStep);

#define __rrateThreshold__ 1E-12
#define __maxRrate__ 2E-4

//===Event Driven Molecular Dynamics Expand diameter===
ExpandEDMD* addExpandEdmd(Box* box, Particle* particle, Update* update, Variable* var);
ExpandEDMD* getExpandEdmd(Box* box, Particle* particle, Update* update);
int delExpandEdmd(Box* box, Particle* particle, Update* update);
ReturnType colForward_expand(Box* box, Particle* particle, Update* update, int maxStep);

//===Const-Z (both P and T) Event Driven Molecular Dynamics===
NptIsoEDMD* addNptEdmd(Box* box, Particle* particle, Update* update, Variable* var);
NptIsoEDMD* getNptEdmd(Box* box, Particle* particle, Update* update);
int delNptEdmd(Box* box, Particle* particle, Update* update);
ReturnType stepForward_npt(Box* box, Particle* particle, Update* update, int maxStep);


#ifdef __cplusplus
}
#endif

#endif
