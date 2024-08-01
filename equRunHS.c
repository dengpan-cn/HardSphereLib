//
//  equRunHS.c
//  eventSim
//
//  Created by Deng Pan on 2024/8/1.
//

#include "EventDrivenSim_ndim.h"
#include "SimEventDrivenSub_ndim.h"
#include "StructEventDriven_ndim.h"
#include "VecMathEventDriven_ndim.h"

int main(int argc, char const *argv[]) {
    MSG = "Equilibrium Running.";
    
    double tic = getTimeStamp();
    char ticStr[32];
    getTimeString(ticStr);
    //=============================================================
    Variable *var = (Variable *)calloc(1, sizeof(Variable));
    readCmdLineVar(var, argc, argv);
    
    Box *box = (Box *)calloc(1, sizeof(Box));
    Particle *particle = (Particle *)calloc(1, sizeof(Particle));
    Update *update = (Update *)calloc(1, sizeof(Update));
    
    readConf(box, particle, update, var);
    moveBarycenter(box, particle, update);
    genGaussianVeloc(box, particle, update, 1.0);
    if (checkOverlap(box, particle, update) == -1) {
        Info("overlapping between particles!");
        return -1;
    }
    
    char strCmd[4096];
    char fname[4096];
    FILE *fout = NULL;
    
    // 0. swap MC moving. The following step 1, 2, 3 will running with this SWAP setting.
    sprintf(strCmd, "--swap");
    addVariable(var, strCmd);
    addSwap(box, particle, update, var);
    //=============================================================
    
    // 1. quick edmd running
    sprintf(strCmd, "--edmd");
    addVariable(var, strCmd);
    SimEDMD *edmd1 = addSimEdmd(box, particle, update);
    do {
        colForward_edmd(box, particle, update, -1.0, 1E3);
        
        if (update->rtype & ThermoReturn) {
        }
        if (update->rtype & JamReturn) {
            Abort("Jam\n");
            break;
        }
        if (update->rtype & ErrorReturn) {
            Abort("Error\n");
            break;
        }
        if (update->rtype & HaltReturn) {
            Abort("Exit\n");
            break;
        }
        if (update->stepCol >= 1E4)
            break;
        
        if (__nCatchSignal__ >= 5) {
            syncAll(box, particle, update);
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            char fname[4096];
            sprintf(fname, "%s/restart_%s.bin", var->cwd, var->sf);
            emergWriteConfHalt(box, particle, update, fname);
            break;
        }
    } while (true);
    syncAll(box, particle, update);
    delSimEdmd(box, particle, update);
    
    if (checkOverlap(box, particle, update) == -1) {
        Info("overlapping between particles!");
        return -1;
    }
    //=============================================================
    
    // 2. compressing to target density
    ExpandEDMD *expd = addExpandEdmd(box, particle, update, var);
    do {
        colForward_expand(box, particle, update, 1E3);
        
        if (update->rtype & ThermoReturn) {
//            double iStep = update->stepCol + (double)update->colCnt / (double)particle->nAtom;
//            double Z = (update->Ztensor[spaceIdx2voigt(0, 0)] + update->Ztensor[spaceIdx2voigt(1, 1)])/2.0;
//            safeFprintf(stdout, "%g %g %g %g\n", iStep, update->Tint, Z, update->volFrac);
        }
        if (update->rtype & JamReturn) {
            Abort("Jam\n");
            break;
        }
        if (update->rtype & ErrorReturn) {
            Abort("Error\n");
            break;
        }
        if ((update->rtype & DensityReturn) || (update->rtype & PressReturn)) {
            break;
        }

        if (__nCatchSignal__ >= 5) {
            syncAll(box, particle, update);
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            char fname[4096];
            sprintf(fname, "%s/restart_%s.bin", var->cwd, var->sf);
            emergWriteConfHalt(box, particle, update, fname);
            break;
        }
    } while (true);
    syncAll(box, particle, update);
    delExpandEdmd(box, particle, update);
    
    if (checkOverlap(box, particle, update) == -1) {
        Info("overlapping between particles!");
        return -1;
    }
    //=============================================================
    
    // 3. quick edmd running
    delVariable(var, "edmd");
    sprintf(strCmd, "--edmd");
    addVariable(var, strCmd);
    SimEDMD *edmd3 = addSimEdmd(box, particle, update);
    do {
        colForward_edmd(box, particle, update, -1.0, 100);
        
        if (update->rtype & ThermoReturn) {
        }
        if (update->rtype & JamReturn) {
            Abort("Jam\n");
            break;
        }
        if (update->rtype & ErrorReturn) {
            Abort("Error\n");
            break;
        }
        if (update->rtype & HaltReturn) {
            Abort("Exit\n");
            break;
        }
        if (update->stepCol >= 1E4)
            break;
        
        if (__nCatchSignal__ >= 5) {
            syncAll(box, particle, update);
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            char fname[4096];
            sprintf(fname, "%s/restart_%s.bin", var->cwd, var->sf);
            emergWriteConfHalt(box, particle, update, fname);
            break;
        }
    } while (true);
    syncAll(box, particle, update);
    delSimEdmd(box, particle, update);
    
    if (checkOverlap(box, particle, update) == -1) {
        Info("overlapping between particles!");
        return -1;
    }
    //=============================================================
    
    // 4-0. swap MC moving
    delSwap(update);
    delVariable(var, "swap");
    sprintf(strCmd, "--swap sisf 7.3 0.1 msd 1.0");
    addVariable(var, strCmd);
    SWAP *swap4 = addSwap(box, particle, update, var);
    
    // 4-1. edmd running to produce equilibrium conf.
    sprintf(fname, "%s/simEquHS_%%%dlf_%s.bin", var->cwd, voigtDIM + 5, var->sf);
    fout = createFileReadWrite(fname);
    
    delVariable(var, "edmd");
    sprintf(strCmd, "--edmd");
    addVariable(var, strCmd);
    SimEDMD *edmd4 = addSimEdmd(box, particle, update);
    do {
        colForward_edmd(box, particle, update, -1.0, 1E3);
        
        if (update->rtype & ThermoReturn) {
            double iStep = update->stepCol + (double)update->colCnt / (double)particle->nAtom;
            double ratio = (double)swap4->naccept / (double)swap4->ntrial;
            safeFwrite(fout, &iStep, sizeof(double), 1);
            safeFwrite(fout, &update->Tint, sizeof(double), 1);
            safeFwrite(fout, update->Ztensor, sizeof(uptriMat), 1);
            safeFwrite(fout, &swap4->selfISF, sizeof(double), 1);
            safeFwrite(fout, &swap4->msd, sizeof(double), 1);
            safeFwrite(fout, &ratio, sizeof(double), 1);
        }
        if (update->rtype & JamReturn) {
            Abort("Jam\n");
            break;
        }
        if (update->rtype & ErrorReturn) {
            Abort("Error\n");
            break;
        }
        if (update->rtype & HaltReturn) {
            Info("Halt due to SWAP.");
            break;
        }
        
        if (__nCatchSignal__ >= 5) {
            syncAll(box, particle, update);
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            char fname[4096];
            sprintf(fname, "%s/restart_%s.bin", var->cwd, var->sf);
            emergWriteConfHalt(box, particle, update, fname);
            break;
        }
    } while (true);
    syncAll(box, particle, update);
    {
        calcKinTensor(box, particle, update);
        calcPressure(box, particle, update);
        calcPoten(box, particle, update);
        
        double iStep = update->stepCol + (double)update->colCnt / (double)particle->nAtom;
        double ratio = (double)swap4->naccept / (double)swap4->ntrial;
        safeFwrite(fout, &iStep, sizeof(double), 1);
        safeFwrite(fout, &update->Tint, sizeof(double), 1);
        safeFwrite(fout, update->Ztensor, sizeof(uptriMat), 1);
        safeFwrite(fout, &swap4->selfISF, sizeof(double), 1);
        safeFwrite(fout, &swap4->msd, sizeof(double), 1);
        safeFwrite(fout, &ratio, sizeof(double), 1);
        
        safeCloseFile(fout);
    }
    delSimEdmd(box, particle, update);
    delSwap(update);

    if (checkOverlap(box, particle, update) == -1) {
        Info("overlapping between particles!");
        return -1;
    }
    //=============================================================
    
    // 5-0. swap MC moving
    delSwap(update);
    delVariable(var, "swap");
    sprintf(strCmd, "--swap");
    addVariable(var, strCmd);
    SWAP *swap5 = addSwap(box, particle, update, var);

    //5-1. NPT running
    sprintf(fname, "%s/simNptHS_%%%dlf_%s.bin", var->cwd, voigtDIM + 4, var->sf);
    fout = createFileReadWrite(fname);
    
    sprintf(strCmd, "--npt 30.0 10.0");
    addVariable(var, strCmd);
    NptIsoEDMD *npt = addNptEdmd(box, particle, update, var);
    do {
        stepForward_npt(box, particle, update, 1E3);
        
        if (update->rtype & ThermoReturn) {
            double iStep = update->stepCol + (double)update->colCnt / (double)particle->nAtom;
            double ratio = (double)swap4->naccept / (double)swap4->ntrial;
            safeFwrite(fout, &iStep, sizeof(double), 1);
            safeFwrite(fout, &update->Tint, sizeof(double), 1);
            safeFwrite(fout, &update->volFrac, sizeof(double), 1);
            safeFwrite(fout, update->Ztensor, sizeof(uptriMat), 1);
            safeFwrite(fout, &ratio, sizeof(double), 1);
        }
        if (update->rtype & JamReturn) {
            Abort("Jam\n");
            break;
        }
        if (update->rtype & ErrorReturn) {
            Abort("Error\n");
            break;
        }
        if (update->rtype & HaltReturn) {
            break;
        }
        if(update->stepCol >= 1E5) break;
        
        if (__nCatchSignal__ >= 5) {
            syncAll(box, particle, update);
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            char fname[4096];
            sprintf(fname, "%s/restart_%s.bin", var->cwd, var->sf);
            emergWriteConfHalt(box, particle, update, fname);
            break;
        }
    } while (true);
    syncAll(box, particle, update);
    {
        double iStep = update->stepCol + (double)update->colCnt / (double)particle->nAtom;
        double ratio = (double)swap4->naccept / (double)swap4->ntrial;
        safeFwrite(fout, &iStep, sizeof(double), 1);
        safeFwrite(fout, &update->Tint, sizeof(double), 1);
        safeFwrite(fout, &update->volFrac, sizeof(double), 1);
        safeFwrite(fout, update->Ztensor, sizeof(uptriMat), 1);
        safeFwrite(fout, &ratio, sizeof(double), 1);
        safeCloseFile(fout);
    }
    delNptEdmd(box, particle, update);
    delSwap(update);
    
    if (checkOverlap(box, particle, update) == -1) {
        Info("overlapping between particles!");
        return -1;
    }
    //=============================================================
    
    writeConf(box, particle, update, var);
    
    //=============================================================
    double toc = getTimeStamp();
    char tocStr[32];
    getTimeString(tocStr);
    safeFprintf(stdout, "Interval: %s (%d) -- %s (%d); [Time: %gs]\n", ticStr, (int)tic, tocStr, (int)toc, toc - tic);
    safeFprintf(logFile, "Interval: %s (%d) -- %s (%d); [Time: %gs]\n", ticStr, (int)tic, tocStr, (int)toc, toc - tic);
    return EXIT_SUCCESS;
}
