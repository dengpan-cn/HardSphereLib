#include "SimEventDrivenSub_ndim.h"

#ifdef __cplusplus
#error "This file must NOT be compiled by C++ compliler!"
#endif

void readCmdLineVar(Variable *var, int argc, char const *argv[]) {
    if (setSignalHandler() < 0)
        Abort("Failed when call setSignalHandler()!");
    
    if (screenOutputReadConf)
        safeFprintf(stderr, "Reading cmdline args ...\n");
    // find --trunc or --append tag
    bool nocite = true;
    for (int iarg = 1; iarg < argc; iarg++) {
        if (!strcmp(argv[iarg], "--trunc")) {
            if (truncFileFlag != 0)
                Abort("--trunc or --append");
            truncFileFlag = 1;
        } else if (!strcmp(argv[iarg], "--append")) {
            if (truncFileFlag != 0)
                Abort("--trunc or --append");
            truncFileFlag = 2;
        } else if (!strcmp(argv[iarg], "--cite")) {
            nocite = false;
            screenOutputReadConf += 10000;
        } else if (!strcmp(argv[iarg], "--dumpSourceFile")) {
#ifdef __Linux__
            dumpSourceFile();
#endif
            exit(EXIT_SUCCESS);
        } else if (!strcmp(argv[iarg], "--help")) {
            logFile = stdout;
            safeFprintf(logFile, "Info at compiling time:\n");
            safeFprintf(logFile, "\tCompile time: %s %s\n", __DATE__, __TIME__);
#ifdef __CompilePath__
            safeFprintf(logFile, "\tCompile path: %s\n", __CompilePath__);
#endif
#ifdef __SourceFileName__
            safeFprintf(logFile, "\tSource Files: %s\n", __SourceFileName__);
#endif
#if defined(__triBox__)
            safeFprintf(logFile, "\tBox: triclinic\n");
#elif defined(__orthBox__)
            safeFprintf(logFile, "\tBox: orthogonal\n");
#endif
            safeFprintf(logFile, "\tDIM: %d\n", DIM);
            safeFprintf(logFile, "\tdumpFileRevNum: %d\n", dumpFileRevNum);
            safeFprintf(logFile, "\tInteraction between particles: Hard Sphere.\n");
            
            char *ptr = getcwd(NULL, 0);
            safeFprintf(logFile, "\tRunning path: %s\n", ptr);
            safeFree(ptr);
            
            exit(EXIT_SUCCESS);
        }
    }
    
    // parse cmd line
    var->maxVar = 8;
    var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
    char *msg = NULL;
    // variable start with "--" and the next character is not "-".
    for (int iarg = 1; iarg < argc;) {
        if (!strcmp(argv[iarg], "--log")) {
            iarg += 1;
            if (iarg < argc && strncmp(argv[iarg], "--", 2)) {
                msg = (char *)argv[iarg];
                iarg += 1;
            }
        } else if (!strcmp(argv[iarg], "--trunc") ||
                   !strcmp(argv[iarg], "--append")) {
            iarg++;
        } else if (!strcmp(argv[iarg], "--cite")) {
            iarg++;
        } else if (!strcmp(argv[iarg], "--cwd")) {
            if (iarg + 1 >= argc)
                Abort("--cwd .");
            int len = (int)strlen(argv[iarg + 1]) + 2;
            var->cwd = (char *)calloc(len, sizeof(char));
            sprintf(var->cwd, "%s", argv[iarg + 1]);
            iarg += 2;
        } else if (!strcmp(argv[iarg], "--sf")) {
            if (iarg + 1 >= argc)
                Abort("--sf id");
            int len = (int)strlen(argv[iarg + 1]) + 2;
            var->sf = (char *)calloc(len, sizeof(char));
            sprintf(var->sf, "%s", argv[iarg + 1]);
            iarg += 2;
        } else if (!strncmp(argv[iarg], "--", 2)) {
            if (strlen(argv[iarg]) == 2) {
                Abort("--varName varArg ...");
            } else if (argv[iarg][2] == '-') {
                Abort("--varName varArg ...");
            }
            
            if (findVariable(var, (char *)argv[iarg]))
                Abort("Repetitive cmd: %s", argv[iarg]);
            if (var->nVar == var->maxVar) {
                var->maxVar += 8;
                var->cmd = (cmdArg *)realloc(var->cmd, var->maxVar * sizeof(cmdArg));
            }
            
            int nsize = (int)strlen(argv[iarg]) + 1;
            var->cmd[var->nVar].cmdArgc = 0;
            var->cmd[var->nVar].cmdArgv = NULL;
            int cmdArgvStart = iarg + 1;
            
            iarg++;
            while (iarg < argc) {
                bool isVar = true;
                if (strncmp(argv[iarg], "--", 2))
                    isVar = false;
                if (isVar) {
                    if (strlen(argv[iarg]) <= 2)
                        isVar = false;
                    else if (argv[iarg][2] == '-')
                        isVar = false;
                }
                if (isVar)
                    break;
                
                var->cmd[var->nVar].cmdArgc++;
                nsize += strlen(argv[iarg]) + 1;
                iarg++;
            }
            if (var->cmd[var->nVar].cmdArgc != 0) {
                var->cmd[var->nVar].cmdArgv =
                (char **)calloc(var->cmd[var->nVar].cmdArgc, sizeof(char *));
            }
            char *ptr = (char *)calloc(nsize + 5, sizeof(char));
            
            var->cmd[var->nVar].cmdType = ptr;
            memcpy(ptr, argv[cmdArgvStart - 1], (strlen(argv[cmdArgvStart - 1]) + 1) * sizeof(char));
            ptr += strlen(argv[cmdArgvStart - 1]) + 1;
            for (int ith = cmdArgvStart; ith < iarg; ith++) {
                var->cmd[var->nVar].cmdArgv[ith - cmdArgvStart] = ptr;
                memcpy(ptr, argv[ith], (strlen(argv[ith]) + 1) * sizeof(char));
                ptr += strlen(argv[ith]) + 1;
            }
            
            var->nVar++;
        } else
            Abort("Unrecognized cmd %s!", argv[iarg]);
    }
    
    if (var->cwd == NULL) {
        var->cwd = (char *)calloc(3, sizeof(char));
        sprintf(var->cwd, ".");
    }
    if (var->sf == NULL) {
        var->sf = (char *)calloc(9, sizeof(char));
        sprintf(var->sf, "default");
    }
    
    // logfile
    if (!nocite) {
        char timeString[32];
        getTimeString(timeString);
        
        char fname[40960] = {0};
        snprintf(fname, 40960, "%s/logFile_%s_%d_%s.dat", var->cwd, timeString, getpid(), var->sf);
        
        logFile = createFileReadWrite(fname);
        safeFprintf(logFile, "Info at compiling time:\n");
        safeFprintf(logFile, "\tCompile time: %s %s\n", __DATE__, __TIME__);
#ifdef __CompilePath__
        safeFprintf(logFile, "\tCompile path: %s\n", __CompilePath__);
#endif
#ifdef __SourceFileName__
        safeFprintf(logFile, "\tSource Files: %s\n", __SourceFileName__);
#endif
#if defined(__triBox__)
        safeFprintf(logFile, "\tBox: triclinic\n");
#elif defined(__orthBox__)
        safeFprintf(logFile, "\tBox: orthogonal\n");
#endif
        safeFprintf(logFile, "\tDIM: %d\n", DIM);
        safeFprintf(logFile, "\tdumpFileRevNum: %d\n", dumpFileRevNum);
        safeFprintf(logFile, "\tInteraction between particles: Hard Sphere.\n");
        
        safeFprintf(logFile, "Info at running time:\n");
        safeFprintf(logFile, "\tProg. Name: %s\n", argv[0]);
        safeFprintf(logFile, "\tRunning time: %s\n", timeString);
        safeFprintf(logFile, "\tCmdLine arguments: ");
        for (int iarg = 1; iarg < argc; iarg++) {
            safeFprintf(logFile, "%s ", argv[iarg]);
        }
        safeFprintf(logFile, "\n");
        
        char *ptr = getcwd(NULL, 0);
        safeFprintf(logFile, "\tRunning path: %s\n", ptr);
        safeFree(ptr);
        
        if (MSG)
            safeFprintf(logFile, "Readme of Program:\n\t%s\n", MSG);
        if (msg)
            safeFprintf(logFile, "Extra info from cmdline: \n\t%s\n", msg);
    }
    
    if (screenOutputReadConf) {
        // screen output
        safeFprintf(stderr, "Prog. Name: \n\t%s\n", argv[0]);
        safeFprintf(stderr, "CmdLine arguments: \n\t");
        for (int iarg = 1; iarg < argc; iarg++) {
            safeFprintf(stderr, "%s ", argv[iarg]);
        }
        safeFprintf(stderr, "\n");
        
        char *ptr = getcwd(NULL, 0);
        safeFprintf(stderr, "Running path: \n\t%s\n", ptr);
        safeFree(ptr);
        if (MSG)
            safeFprintf(stderr, "Readme of Program:\n\t%s\n", MSG);
        if (msg)
            safeFprintf(stderr, "Extra info from cmdline: \n\t%s\n", msg);
    }
}

//===============================================================================
void initConfInfo(Box *box, Particle *particle, Update *update) {
    double vol = 0;
    double minScale = particle->diameterScale[0],
    maxScale = particle->diameterScale[0];
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vol += VolUnitSphere * pow(particle->diameterScale[iatom], DIM);
        double tmp = particle->diameterScale[iatom];
        minScale = (minScale < tmp ? minScale : tmp);
        maxScale = (maxScale > tmp ? maxScale : tmp);
    }
    update->nebrList.minDiameterScale = minScale;
    update->nebrList.maxDiameterScale = maxScale;
    
    update->volFrac = vol * pow(particle->meanDiameter * 0.5, DIM) / box->volume;
    update->dof = DIM * particle->nAtom;
}

void readConf_dump(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--rd");
    if (!cmd || cmd->cmdArgc != 2)
        Abort("--rd dump.bin whichStep");
    int whichStep = (int)atoi(cmd->cmdArgv[1]);
    
    mmapBinFile *binFile = openBinFile(cmd->cmdArgv[0]);
    if (whichStep < 0)
        whichStep = binFile->nStep - 1;
    if (whichStep >= binFile->nStep) {
        Abort("Step %d is out of Range: [-1,%d];", whichStep, binFile->nStep - 1);
    }
    
    readDump(box, particle, update, binFile, whichStep);
    
    closeBinFile(&binFile);
}

void readConf_data(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "--rf");
    if (!cmd || cmd->cmdArgc == 0)
        Abort("read data: --rf lmp.bin");
    FILE *fp = openExistFileReadOnly(cmd->cmdArgv[0]);
    if (particle->pos != NULL)
        Abort("--rf or --rd");
    
    char str[4096];
    fgets(str, 4096, fp);
    if (strcmp(str, "binary") == 0) {
        bool hasMeanDiameter = false, hasDimension = false;
        while (fread(str, sizeof(char), 32, fp) == 32) {
            if (strcmp(str, "dimension") == 0) {
                fread(&box->dim, sizeof(int), 1, fp);
                if (box->dim != DIM)
                    Abort("The file is for d = %d, while the code is for d = %d!",
                          box->dim, DIM);
                hasDimension = true;
            } else if (strcmp(str, "atoms") == 0) {
                fread(&particle->nAtom, sizeof(int), 1, fp);
                if (particle->nAtom <= 0)
                    Abort("Wrong File!");
                
                particle->pos =
                (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
                particle->veloc =
                (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
                particle->timeStamp = (double *)calloc(particle->nAtom, sizeof(double));
                particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
                particle->img = (intVector *)calloc(particle->nAtom, sizeof(intVector));
                particle->type = (int *)calloc(particle->nAtom, sizeof(int));
                particle->diameterScale =
                (double *)calloc(particle->nAtom, sizeof(double));
                particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
                particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
            } else if (strcmp(str, "atom types") == 0) {
                fread(&particle->nAtomType, sizeof(int), 1, fp);
                if (particle->nAtomType <= 0)
                    Abort("Wrong File!");
                
                particle->massPerType =
                (double *)calloc(particle->nAtomType, sizeof(double));
                for (int itype = 0; itype < particle->nAtomType; itype++)
                    particle->massPerType[itype] = 1.0;
            } else if (strcmp(str, "box Hvoigt") == 0) {
                if (!hasDimension) {
                    box->dim = 3;
                    if (box->dim != DIM)
                        Abort("Dimension is not consistent!");
                    
                    fread(&box->boxH[spaceIdx2voigt(0, 0)], sizeof(double), 1, fp); // xx
                    fread(&box->boxH[spaceIdx2voigt(1, 1)], sizeof(double), 1, fp); // yy
                    fread(&box->boxH[spaceIdx2voigt(2, 2)], sizeof(double), 1, fp); // zz
                    fread(&box->boxH[spaceIdx2voigt(1, 2)], sizeof(double), 1, fp); // yz
                    fread(&box->boxH[spaceIdx2voigt(0, 2)], sizeof(double), 1, fp); // xz
                    fread(&box->boxH[spaceIdx2voigt(0, 1)], sizeof(double), 1, fp); // xy
                } else {
                    fread(&box->boxH, sizeof(uptriMat), 1, fp);
                }
            } else if (strcmp(str, "mean diameter") == 0) {
                fread(&particle->meanDiameter, sizeof(double), 1, fp);
                hasMeanDiameter = true;
            } else if (strcmp(str, "Atoms") == 0) {
                for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                    fread(&particle->type[iatom], sizeof(int), 1, fp);
                    fread(&particle->diameterScale[iatom], sizeof(double), 1, fp);
                    fread(&particle->pos[iatom], sizeof(doubleVector), 1, fp);
                    fread(&particle->img[iatom], sizeof(intVector), 1, fp);
                    
                    particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
                    
                    particle->tag2id[iatom] = iatom;
                    particle->id2tag[iatom] = iatom;
                }
            } else if (strcmp(str, "Velocities") == 0) {
                for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                    fread(&particle->veloc[iatom], sizeof(doubleVector), 1, fp);
                }
            } else {
                int sizeByte = 0;
                fread(&sizeByte, sizeof(int), 1, fp);
                void *ptr = calloc(sizeByte, 1);
                fread(ptr, 1, sizeByte, fp);
                if (findToolkit(&update->toolkit, str) >= 0)
                    Abort("Wrong binary file!");
                addToolkit(&update->toolkit, ptr, NULL, str);
            }
        }
        if (!hasMeanDiameter) {
            double meanDiameter = 0.0;
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                meanDiameter += particle->diameterScale[iatom];
            }
            meanDiameter = meanDiameter / particle->nAtom;
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                particle->diameterScale[iatom] /= meanDiameter;
            }
            particle->meanDiameter = meanDiameter;
        }
    } else
        Abort("Wrong File!");
    safeCloseFile(fp);
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    particle->isSync = true;
    update->Edone = update->Zdone = update->Tdone = false;
    update->nebrList.isValid = false;
    update->nebrList.doSort = true;
    update->eventList.isValid = false;
    
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
    adjustImg(box, particle);
}
void readConf(Box *box, Particle *particle, Update *update, Variable *var) {
    int cntR = 0;
    cmdArg *cmd = findVariable(var, "--rf");
    if (cmd) {
        readConf_data(box, particle, update, var);
        cntR++;
    }
    cmd = findVariable(var, "--rd");
    if (cmd && cmd->cmdArgc == 2) {
        if (cntR == 1)
            Abort("--rf or --rd");
        readConf_dump(box, particle, update, var);
        cntR++;
    }
    if (cntR != 1)
        Abort("--rd dump.bin step or --rf conf.bin"); //--rd dump.bin or no "--rd"
    // and "--rf"
    
    //===get global setting===
    cmd = findVariable(var, "skin");
    if (cmd) {
        if (cmd->cmdArgc <= 0)
            Abort("--skin 0.3");
        update->nebrList.skinSet = atof(cmd->cmdArgv[0]);
    } else {
        // estimated from spaces between particles
        double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
        double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
        double lmsset = (lbox - lsph) / (update->nebrList.minDiameterScale * particle->meanDiameter);
        // estimated from box side length
        doubleVector distPlane;
        calcDistBoxPlane(distPlane, box->boxEdge);
        double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
        for (int idim = 0; idim < DIM; idim++) {
            minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
        }
        double maxRcut = update->nebrList.maxDiameterScale * particle->meanDiameter;
        double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
        double maxRskinSet = maxRskin / (update->nebrList.minDiameterScale * particle->meanDiameter);
        
        update->nebrList.skinSet = cpuMin(lmsset, maxRskinSet);
    }
    
    cmd = findVariable(var, "zero");
    if (!cmd) {
        Info("--zero ZMperiod ZCperiod;\n"
             "\tPerforming zero total momentum every ZMperiod step.\n"
             "\tMoving center of mass to origin every ZCperiod step.\n"
             "\tNegative value means Never performing these two action.\n"
             "\tdefault: --zero 1e3 1e4.");
        update->ZMperiod = 1000;
        update->ZCperiod = 10000;
    } else {
        update->ZMperiod = (int)atof(cmd->cmdArgv[0]);
        update->ZCperiod = (int)atof(cmd->cmdArgv[1]);
    }
    
    cmd = findVariable(var, "update");
    if (!cmd) {
        Info("--update updateOutputStep; default --update 1E3");
        update->outputPeriod = 1000;
    } else {
        update->outputPeriod = (int)atof(cmd->cmdArgv[0]);
    }
    //===get global setting===
    update->kTt = 1.0;
    
#ifdef __orthBox__
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            if (fabs(box->boxH[spaceIdx2voigt(idim, jdim)]) >= 5E-16)
                Abort("Not orthogonal Box!");
            box->boxH[spaceIdx2voigt(idim, jdim)] = 0;
        }
    }
#endif
    
    if (screenOutputReadConf) {
        printf("===========System Info==========\n");
        printf("Dimension: %d\n", DIM);
        printf("Number of Particles: %d;\n", particle->nAtom);
        printf("Volume Fraction: %g;\n", update->volFrac);
        printf("Min(diameter): %g;\nMax(diameter): %g;\n",
               update->nebrList.minDiameterScale,
               update->nebrList.maxDiameterScale);
        printf("Edges of Simulation Box: \n");
        for (int iedge = 0; iedge < DIM; iedge++) {
            printf("\t");
            for (int jdim = 0; jdim < DIM; jdim++) {
                printf("%-8.6e\t", box->boxEdge[iedge][jdim] / update->distanceUnits);
            }
            printf("\n");
        }
        printf("Interaction between particles: Hard Sphere.\n");
        printf("===========System Info==========\n");
        
        if (screenOutputReadConf >= 10000)
            screenOutputReadConf -= 10000;
    }
}

void writeConf(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "wf");
    if (!cmd)
        return;
    if (cmd->cmdArgc == 0)
        Abort("--wf output.bin");
    
    syncAll(box, particle, update);
    
    FILE *fbin = createFileReadWrite(cmd->cmdArgv[0]);
    {
        char str[4096];
        // header
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "binary");
        str[31] = '\n';
        fwrite(str, sizeof(char), 32, fbin); // 32 byte;
        // dimension
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "dimension");
        fwrite(str, sizeof(char), 32, fbin); // 32 byte;
        fwrite(&box->dim, sizeof(int), 1, fbin);
        // nAtom atoms
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "atoms");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&particle->nAtom, sizeof(int), 1, fbin);
        // nAtomType atom types
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "atom types");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&particle->nAtomType, sizeof(int), 1, fbin);
        // box Hvoigt
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "box Hvoigt");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&box->boxH, sizeof(uptriMat), 1, fbin);
        // mean diameter
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "mean diameter");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&particle->meanDiameter, sizeof(double), 1, fbin);
        // Atoms Info
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "Atoms");
        fwrite(str, sizeof(char), 32, fbin);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            
            int type = particle->type[idx];
            double diaScale = particle->diameterScale[idx];
            doubleVecPtr posPtr = particle->pos[idx];
            intVecPtr imgPtr = particle->img[idx];
            
            fwrite(&type, sizeof(int), 1, fbin);
            fwrite(&diaScale, sizeof(double), 1, fbin);
            fwrite(posPtr, sizeof(doubleVector), 1, fbin);
            fwrite(imgPtr, sizeof(intVector), 1, fbin);
        }
        // Velocities Info
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "Velocities");
        fwrite(str, sizeof(char), 32, fbin);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            
            doubleVecPtr velocPtr = particle->veloc[idx];
            
            fwrite(velocPtr, sizeof(doubleVector), 1, fbin);
        }
    }
    safeCloseFile(fbin);
}
//write configurations then halt.
int emergWriteConfHalt(Box *box, Particle *particle, Update *update, char *fname){
    if (access(fname, F_OK) == 0) {
        char tsring[32];
        getTimeString(tsring);
        Info("File \"%s\" will be truncated! Operation Time: %s", fname, tsring);
    }
    syncAll(box, particle, update);
    
    FILE *fbin = fopen(fname, "wb+");
    {
        char str[4096];
        // header
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "binary");
        str[31] = '\n';
        fwrite(str, sizeof(char), 32, fbin); // 32 byte;
        // dimension
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "dimension");
        fwrite(str, sizeof(char), 32, fbin); // 32 byte;
        fwrite(&box->dim, sizeof(int), 1, fbin);
        // nAtom atoms
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "atoms");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&particle->nAtom, sizeof(int), 1, fbin);
        // nAtomType atom types
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "atom types");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&particle->nAtomType, sizeof(int), 1, fbin);
        // box Hvoigt
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "box Hvoigt");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&box->boxH, sizeof(uptriMat), 1, fbin);
        // mean diameter
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "mean diameter");
        fwrite(str, sizeof(char), 32, fbin);
        fwrite(&particle->meanDiameter, sizeof(double), 1, fbin);
        // Atoms Info
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "Atoms");
        fwrite(str, sizeof(char), 32, fbin);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            
            int type = particle->type[idx];
            double diaScale = particle->diameterScale[idx];
            doubleVecPtr posPtr = particle->pos[idx];
            intVecPtr imgPtr = particle->img[idx];
            
            fwrite(&type, sizeof(int), 1, fbin);
            fwrite(&diaScale, sizeof(double), 1, fbin);
            fwrite(posPtr, sizeof(doubleVector), 1, fbin);
            fwrite(imgPtr, sizeof(intVector), 1, fbin);
        }
        // Velocities Info
        memset(str, '\0', 4096 * sizeof(char));
        sprintf(str, "Velocities");
        fwrite(str, sizeof(char), 32, fbin);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            
            doubleVecPtr velocPtr = particle->veloc[idx];
            
            fwrite(velocPtr, sizeof(doubleVector), 1, fbin);
        }
        
        fflush(fbin);
    }
    safeCloseFile(fbin);
    
    Abort("Halt.");
    
    return 0;
}

int checkMemSortFlag(Box *box, Particle *particle, Update *update){
    if (!particle->__isSortForbidden) {
        return 1;
    }
    if (particle->nAtom == 0) {
        return 0;
    }
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        if (particle->tag2id[iatom] != iatom) {
            Abort("Abort due to wrong __isSortForbidden.")
        }
    }
    
    return 1;
}
int turnOnMemSort(Box *box, Particle *particle, Update *update){
    checkMemSortFlag(box, particle, update);
    
    particle->__isSortForbidden = false;
    
    return 1;
}
int turnOffMemSort(Box *box, Particle *particle, Update *update){
    particle->__isSortForbidden = true;
    checkMemSortFlag(box, particle, update);
    
    return 1;
}

//===============================================================================
void setBoxPara(Box *box) {
    box->dim = DIM;
    box->volume = 1.0;
    for (int idim = 0; idim < DIM; idim++) {
        if (box->boxH[spaceIdx2voigt(idim, idim)] <= 0.0) {
            Abort("Not rihgt hand basis!");
        }
        box->volume *= box->boxH[spaceIdx2voigt(idim, idim)];
        
        vZeros(box->boxEdge[idim]);
        for (int jdim = 0; jdim <= idim; jdim++) {
            box->boxEdge[idim][jdim] = box->boxH[spaceIdx2voigt(jdim, idim)];
        }
    }
    
    vZeros(box->cornerLo);
    for (int idim = 0; idim < DIM; idim++) {
        vAdd(box->cornerLo, box->cornerLo, box->boxEdge[idim]);
    }
    vScale(box->cornerLo, -0.5, box->cornerLo);
    vScale(box->cornerHi, -1.0, box->cornerLo);
    
    uptriMatZeros(box->invBoxH);
    for (int idim = DIM - 1; idim >= 0; idim--) {
        box->invBoxH[spaceIdx2voigt(idim, idim)] =
        1.0 / box->boxH[spaceIdx2voigt(idim, idim)];
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            double cij = 0.0;
            for (int k = idim + 1; k <= jdim; k++) {
                cij += box->boxH[spaceIdx2voigt(idim, k)] *
                box->invBoxH[spaceIdx2voigt(k, jdim)];
            }
            box->invBoxH[spaceIdx2voigt(idim, jdim)] =
            -cij / box->boxH[spaceIdx2voigt(idim, idim)];
        }
    }
}
void calcBoxLenAngle(Box *box, uptriMat para) {
    for (int idim = 0; idim < DIM; idim++) {
        para[spaceIdx2voigt(idim, idim)] = box->boxH[spaceIdx2voigt(idim, idim)];
        for (int jdim = idim + 1; jdim < DIM; jdim++) {
            double len1 = sNorm(box->boxEdge[idim]);
            double len2 = sNorm(box->boxEdge[jdim]);
            double c12 = sDot(box->boxEdge[idim], box->boxEdge[jdim]);
            para[spaceIdx2voigt(idim, jdim)] = acos(c12 / (len1 * len2)) / PI * 180.0;
        }
    }
}
void calcDistBoxPlane(doubleVector distPlane, doubleVector boxEdge[DIM]) {
    doubleVector edge[DIM];
    
    for (int idim = 0; idim < DIM; idim++) {
        vZeros(edge[idim]);
        edge[idim][idim] = 1.0;
    }
    for (int idim = DIM - 1; idim >= 0; idim--) {
        for (int ith = idim; ith < DIM - 1; ith++) {
            vCpy(edge[ith], boxEdge[ith + 1]);
            memset(edge[ith], '\0', idim * sizeof(double));
        }
        vCpy(edge[DIM - 1], boxEdge[idim]);
        memset(edge[DIM - 1], '\0', idim * sizeof(double));
        
        for (int ith = idim; ith < DIM; ith++) {
            vUnit(edge[ith], edge[ith]);
            for (int jth = ith + 1; jth < DIM; jth++) {
                double dotProd = sDot(edge[ith], edge[jth]);
                vScaleAdd(edge[jth], edge[jth], -dotProd, edge[ith]);
            }
        }
        
        distPlane[idim] = sDot(boxEdge[idim], edge[DIM - 1]);
        distPlane[idim] = fabs(distPlane[idim]);
    }
}
void setUnits(Update *update, double distanceUnits) {
    update->massUnits = 1.0;
    update->energyUnits = 1.0;
    update->distanceUnits = distanceUnits;
    
    update->timeUnits = sqrt(update->massUnits / update->energyUnits) * update->distanceUnits;
    update->forceUnits = update->energyUnits / update->distanceUnits;
    update->velocityUnits = sqrt(update->energyUnits / update->massUnits);
    update->pressureUnits = update->energyUnits / pow(update->distanceUnits, DIM);
    update->volumeUnits = pow(update->distanceUnits, DIM);
}
void reInitSim(Box *box, Particle *particle, Update *update) {
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    update->Edone = update->Zdone = update->Tdone = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.doSort = true;
    
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
}

//===============================================================================
void syncAll(Box *box, Particle *particle, Update *update) {
    if (particle->isSync)
        return;
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        double deltaTime = (update->currentStamp - particle->timeStamp[iatom]);
        vScaleAdd(particle->pos[iatom], particle->pos[iatom], deltaTime,
                  particle->veloc[iatom]);
    }
    if (!particle->isSizeFixed) {
        double sfact = (1.0 + update->rrateSet / update->timeUnits * update->currentStamp);
        particle->meanDiameter *= sfact;
        update->volFrac *= pow(sfact, DIM);
        setUnits(update, particle->meanDiameter);
    }
    memset(particle->timeStamp, '\0', particle->nAtom * sizeof(double));
    
    update->currentStamp = 0.0;
    particle->isSync = true;
    update->eventList.isValid = false;
}
void adjustImg(Box *box, Particle *particle) {
    if (!particle->isSync)
        Abort("Please Sync the particles!");
    
#ifdef __triBox__
    bool isFlip = false;
    doubleVector maxTilt;
    for (int idim = 0; idim < DIM; idim++) {
        maxTilt[idim] = 0.505 * box->boxH[spaceIdx2voigt(idim, idim)];
    }
    doubleVector boxEdge[DIM];
    memcpy(boxEdge, box->boxEdge, DIM * sizeof(doubleVector));
    for (int idim = DIM - 1; idim >= 0; idim--) {
        for (int tilt = idim - 1; tilt >= 0; tilt--) {
            if (boxEdge[idim][tilt] >= maxTilt[tilt]) {
                vSub(boxEdge[idim], boxEdge[idim], boxEdge[tilt]);
                isFlip = true;
            } else if (boxEdge[idim][tilt] < -maxTilt[tilt]) {
                vAdd(boxEdge[idim], boxEdge[idim], boxEdge[tilt]);
                isFlip = true;
            }
        }
    }
    if (isFlip) {
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            unwrapPos(particle->pos[iatom], particle->pos[iatom],
                      particle->img[iatom], box->boxH);
        }
        memset(particle->img, '\0', particle->nAtom * sizeof(intVector));
        
        for (int idim = 0; idim < DIM; idim++) {
            for (int jdim = idim; jdim < DIM; jdim++) {
                box->boxH[spaceIdx2voigt(idim, jdim)] = boxEdge[jdim][idim];
            }
        }
        setBoxPara(box);
    }
#endif
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr posPtr = particle->pos[iatom];
        intVecPtr imgPtr = particle->img[iatom];
        
        // lamda_O = invH * r; lamda_lo = lamda_O + (0.5, 0.5, 0.5)
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, posPtr);
        vShiftAll(lamda, 0.5);
        
        intVector shiftImg;
        vFloor(shiftImg, lamda);
        vAdd(imgPtr, imgPtr, shiftImg);
        
        doubleVector shiftPos;
        MatMulVec(shiftPos, box->boxH, shiftImg);
        vSub(posPtr, posPtr, shiftPos);
    }
}
//===============================================================================

//===============================================================================
void calcKinTensor(Box *box, Particle *particle, Update *update) {
    if (update->Tdone)
        return;
    
    uptriMat totKinTensor;
    uptriMatZeros(totKinTensor);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        double m = particle->mass[iatom];
        for (int idim = 0; idim < DIM; idim++) {
            for (int jdim = idim; jdim < DIM; jdim++) {
                totKinTensor[spaceIdx2voigt(idim, jdim)] +=
                m * particle->veloc[iatom][idim] * particle->veloc[iatom][jdim];
            }
        }
    }
    double totKin = 0;
    for (int idim = 0; idim < DIM; idim++) {
        totKin += totKinTensor[spaceIdx2voigt(idim, idim)];
    }
    
    update->Tint = totKin / update->dof;
    uptriMatCpy(update->Kintensor, totKinTensor);
    update->Tdone = true;
}
void calcPressure(Box *box, Particle *particle, Update *update) {
    if (update->Zdone)
        return;
    if (!update->Tdone) {
        calcKinTensor(box, particle, update);
    }
    
    double deltaTime = update->accTimePeriod; // update->duringTime;
    double N_kB_T = (particle->nAtom) * update->Tint;
    
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            update->Ztensor[spaceIdx2voigt(idim, jdim)] = update->Kintensor[spaceIdx2voigt(idim, jdim)] / N_kB_T + update->sColVirialTensor[spaceIdx2voigt(idim, jdim)] / N_kB_T / deltaTime;
        }
    }
    
    update->Z = 0.0;
    for (int idim = 0; idim < DIM; idim++) {
        update->Z += update->Ztensor[spaceIdx2voigt(idim, idim)];
    }
    update->Z = update->Z / DIM;
    
    update->Zdone = true;
    uptriMatZeros(update->sColVirialTensor);
    update->accTimePeriod = 0;
    
    // printf("%g %g (%d %d)\n",deltaTime,update->Z, update->stepCol,update->colCnt);
}
void calcPoten(Box *box, Particle *particle, Update *update) {
    if (update->Edone)
        return;
    
    update->Epair = 0;
    update->Edone = true;
}

void thermostat(Box *box, Particle *particle, Update *update) {
    syncAll(box, particle, update);
    
    if (!update->Tdone)
        calcKinTensor(box, particle, update);
    
    double sfact = sqrt(update->kTt / update->Tint);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vScale(particle->veloc[iatom], sfact, particle->veloc[iatom]);
    }
    
    update->nextTtstep = update->stepCol + update->Ttperiod;
    update->Tdone = false;
    update->Zdone = false;
    update->eventList.isValid = false;
}
void zeroMomentum(Box *box, Particle *particle, Update *update) {
    if (update->ZMperiod < 0) {
        return;
    }
    syncAll(box, particle, update);
    
    doubleVector sV;
    vZeros(sV);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vScaleAdd(sV, sV, particle->mass[iatom], particle->veloc[iatom]);
    }
    vScale(sV, 1.0 / particle->nAtom, sV);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vScaleAdd(particle->veloc[iatom], particle->veloc[iatom],
                  -1.0 / particle->mass[iatom], sV);
    }
    
    update->nextZMstep = update->stepCol + update->ZMperiod;
    update->Tdone = false;
    update->eventList.isValid = false;
}
void moveBarycenter(Box *box, Particle *particle, Update *update) {
    if (update->ZCperiod < 0) {
        return;
    }
    syncAll(box, particle, update);
    
    doubleVector sCent;
    vZeros(sCent);
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector imgPos;
        unwrapPos(imgPos, particle->pos[iatom], particle->img[iatom], box->boxH);
        vAdd(sCent, imgPos, sCent);
    }
    vScale(sCent, 1.0 / particle->nAtom, sCent);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        vSub(particle->pos[iatom], particle->pos[iatom], sCent);
    }
    
    update->nextZCstep = update->stepCol + update->ZCperiod;
    update->nebrList.isValid = false;
    update->eventList.isValid = false;
}

//=========================
void genGaussianVeloc(Box *box, Particle *particle, Update *update, double targT) {
    syncAll(box, particle, update);
    // generate
    {
        double totKin = 0.0;
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            for (int idim = 0; idim < DIM; idim++) {
                particle->veloc[iatom][idim] = rndStdNorm();
            }
            
            totKin += sNormP2(particle->veloc[iatom]);
        }
        update->Tdone = false;
        update->Zdone = false;
    }
    // zero total momentum
    {
        doubleVector sV;
        vZeros(sV);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            vScaleAdd(sV, sV, particle->mass[iatom], particle->veloc[iatom]);
        }
        vScale(sV, 1.0 / particle->nAtom, sV);
        
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            vScaleAdd(particle->veloc[iatom], particle->veloc[iatom],
                      -1.0 / particle->mass[iatom], sV);
        }
        update->Tdone = false;
        update->Zdone = false;
    }
    // scale velocity
    {
        calcKinTensor(box, particle, update);
        
        double sfact = sqrt(targT / update->Tint);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            vScale(particle->veloc[iatom], sfact, particle->veloc[iatom]);
        }
        
        update->Tdone = false;
        update->Zdone = false;
        update->eventList.isValid = false;
    }
}
void instant_inflate(Box *box, Particle *particle, Update *update, double deltaVF) {
    double volfrac_target = update->volFrac + deltaVF;
    double sfact = pow(volfrac_target / update->volFrac, 1.0 / DIM);
    particle->meanDiameter *= sfact;
    update->volFrac = pow(sfact, DIM) * update->volFrac;
    
    setUnits(update, particle->meanDiameter);
    
    update->Edone = update->Zdone = 0;
    update->eventList.isValid = false;
}

//=================
int checkOverlap(Box *box, Particle *particle, Update *update) {
    if (!particle->isSync)
        Abort("Not Sync!");
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        for (int jatom = iatom + 1; jatom < particle->nAtom; jatom++) {
            doubleVector dRij;
            vSub(dRij, particle->pos[iatom], particle->pos[jatom]);
            PBC(dRij, box);
            double sRc =
            (particle->diameterScale[iatom] + particle->diameterScale[jatom]) *
            particle->meanDiameter * 0.5;
            double rij = sNorm(dRij);
            if (rij <= sRc - 10.0 * DBL_EPSILON) {
                safeFprintf(stderr, "overlap: %d %d: %g %g => %10.8lf%%\n", iatom, jatom, rij, sRc,
                       (rij - sRc) / sRc * 100.0);
                return -1;
            }
        }
    }
    return 0;
}

//===============================================================================
writeDumpFile *getWriteDumpFile(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__writeDumpFile__");
    if (whichTool < 0)
        return (writeDumpFile *)NULL;
    return (writeDumpFile *)update->toolkit.toolkit[whichTool];
}
writeDumpFile *addWriteDumpFile(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "dump");
    if (!cmd) {
        return (writeDumpFile *)NULL;
    } else if (getWriteDumpFile(update) != NULL) {
        Info("repetitive initWriteDumpFile!");
        return getWriteDumpFile(update);
    }
    writeDumpFile *dinfo = (writeDumpFile *)calloc(sizeof(writeDumpFile), 1);
    addToolkit(&update->toolkit, (void *)dinfo, NULL, "__writeDumpFile__");
    
    if (!particle->isSync) {
        syncAll(box, particle, update);
    }
    
    char fname[4096];
    if (cmd->cmdArgc == 1) {
        sprintf(fname, "%s/dump_%s_%s.bin", var->cwd, cmd->cmdArgv[0], var->sf);
    } else {
        sprintf(fname, "%s/dump_%s.bin", var->cwd, var->sf);
    }
    
    // bool isFileExist = true;
    dinfo->fdump = createFileReadWrite(fname);
    // bool isContinueRun = ((truncFileFlag == 2) && isFileExist);
    {
        // write header
        char str[32];
        memset(str, '\0', 32 * sizeof(char));
        sprintf(str, "Revised (HS) Binary File");
        str[31] = '\n';
        fwrite(str, sizeof(char), 32, dinfo->fdump);
        dinfo->revNum = dumpFileRevNum;
        fwrite(&dinfo->revNum, sizeof(int), 1, dinfo->fdump);
        
        fwrite(&box->dim, sizeof(int), 1, dinfo->fdump);
        fwrite(&particle->nAtom, sizeof(int), 1, dinfo->fdump);
        fwrite(&particle->nAtomType, sizeof(int), 1, dinfo->fdump);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            fwrite(&particle->diameterScale[idx], sizeof(double), 1, dinfo->fdump);
        }
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            fwrite(&particle->type[idx], sizeof(int), 1, dinfo->fdump);
        }
        fflush(dinfo->fdump);
    }
    /*
     else {
     // do minor checks
     {
     fseek(dinfo->fdump, 0, SEEK_SET);
     
     char str[32];
     fread(str, sizeof(char), 32, dinfo->fdump);
     if (strcmp(str, "Revised Binary File") == 0) {
     int dumpRevNum = -1;
     fread(&dumpRevNum, sizeof(int), 1, dinfo->fdump);
     if (dumpRevNum != dumpFileRevNum)
     Abort("The dump file is for Rev. %d, code is for %d.", dumpRevNum,
     dumpFileRevNum);
     int dim = -1;
     fread(&dim, sizeof(int), 1, dinfo->fdump);
     if (dim != box->dim)
     Abort("The dumpfile is for d = %d, dim of initial file is %d.", dim,
     box->dim);
     } else {
     fseek(dinfo->fdump, 0, SEEK_SET);
     if (DIM != 3)
     Abort("The dumpfile is for d = 3, the code is for d = %d.", DIM);
     }
     
     int nAtom = -1;
     fread(&nAtom, sizeof(int), 1, dinfo->fdump);
     if (particle->nAtom != nAtom)
     Abort("Not match!");
     
     int nAtomType = -1;
     fread(&nAtomType, sizeof(int), 1, dinfo->fdump);
     if (particle->nAtomType != nAtomType)
     Abort("Not match!");
     }
     fseek(dinfo->fdump, 0, SEEK_END);
     }
     */
    
    return dinfo;
}
int writeDump(Box *box, Particle *particle, Update *update) {
    writeDumpFile *wdf = getWriteDumpFile(update);
    if (wdf == NULL)
        return -1;
    
    if (particle->isSync) {
        safeFwrite(wdf->fdump, &particle->nAtom, sizeof(int), 1);
        safeFwrite(wdf->fdump, &box->boxH, sizeof(uptriMat), 1);
        safeFwrite(wdf->fdump, &particle->meanDiameter, sizeof(double), 1);
        double timestep = update->stepCol + (double)update->colCnt/(double)particle->nAtom;
        safeFwrite(wdf->fdump, &timestep, sizeof(double), 1);
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            doubleVector uxyz;
            unwrapPos(uxyz, particle->pos[idx], particle->img[idx], box->boxH);
            
            safeFwrite(wdf->fdump, uxyz, sizeof(doubleVector), 1);
        }
        
    } else {
        safeFwrite(wdf->fdump, &particle->nAtom, sizeof(int), 1);
        safeFwrite(wdf->fdump, &box->boxH, sizeof(uptriMat), 1);
        if (!particle->isSizeFixed) {
            double sfact = (1.0 + update->rrateSet / update->timeUnits * update->currentStamp);
            double meanDiameter = particle->meanDiameter * sfact;
            safeFwrite(wdf->fdump, &meanDiameter, sizeof(double), 1);
        } else {
            safeFwrite(wdf->fdump, &particle->meanDiameter, sizeof(double), 1);
        }
        double timestep = update->stepCol + (double)update->colCnt/(double)particle->nAtom;
        safeFwrite(wdf->fdump, &timestep, sizeof(double), 1);
        
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            
            double deltaTime = (update->currentStamp - particle->timeStamp[idx]);
            doubleVector iPos;
            vScaleAdd(iPos, particle->pos[idx], deltaTime, particle->veloc[idx]);
            
            doubleVector uxyz;
            unwrapPos(uxyz, iPos, particle->img[idx], box->boxH);
            
            safeFwrite(wdf->fdump, uxyz, sizeof(doubleVector), 1);
        }
    }
    
    return 0;
}
int delWriteDumpFile(Update *update) {
    writeDumpFile *wdf = getWriteDumpFile(update);
    if (!wdf)
        return -1;
    safeCloseFile(wdf->fdump);
    delToolkit(&update->toolkit, "__writeDumpFile__");
    
    return 0;
}

mmapBinFile *openBinFile(char *fname) {
    mmapBinFile *binFile = (mmapBinFile *)calloc(sizeof(mmapBinFile), 1);
    
    binFile->fd = open(fname, O_RDONLY);
    fstat(binFile->fd, &binFile->binStat);
    binFile->dataSection = mmap(NULL, binFile->binStat.st_size, PROT_READ,
                                MAP_PRIVATE, binFile->fd, 0);
    if (binFile->dataSection == MAP_FAILED) {
        Abort("Map Binary file Failed !\n");
        return NULL;
    }
    
    // header
    char *header = (char *)binFile->dataSection;
    void *data = binFile->dataSection;
    int *nElement = (int *)data;
    if (strcmp(header, "Revised Binary File") == 0) {
        int *revNum = (int *)(header + 32);
        binFile->revNum = revNum[0];
        if (binFile->revNum != 1)
            Abort("Wrong Binary File!");
        
        int *dim = (int *)(revNum + 1);
        if (dim[0] != DIM)
            Abort("The dumpfile is for d = %d, while the code is for d = %d!", dim[0], DIM);
        nElement = (int *)(dim + 1);
        
        Info("revNum of %s: %d. (0: SS 3d; 1: SS; 2: HS;)", fname, binFile->revNum);
    } else if (strcmp(header, "Revised (HS) Binary File") == 0) {
        int *revNum = (int *)(header + 32);
        binFile->revNum = revNum[0];
        if (binFile->revNum != 2)
            Abort("Wrong Binary File!");
        
        int *dim = (int *)(revNum + 1);
        if (dim[0] != DIM)
            Abort("The dumpfile is for d = %d, while the code is for d = %d!", dim[0], DIM);
        nElement = (int *)(dim + 1);
        
        Info("revNum of %s: %d. (0: SS 3d; 1: SS; 2: HS;)", fname, binFile->revNum);
    } else {
        binFile->revNum = 0;
        if (DIM != 3)
            Abort("The dumpfile is for d = 3, while the code is for d = %d!", DIM);
        Info("revNum of %s: 0. (0: SS 3d; 1: SS; 2: HS;)", fname);
    }
    
    switch (binFile->revNum) {
    case 0:
        binFile->headerSize = sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
        break;
    case 1:
        binFile->headerSize = 32 * sizeof(char) + 2 * sizeof(int) + sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
        break;
    case 2:
        binFile->headerSize = 32 * sizeof(char) + 2 * sizeof(int) + sizeof(int) + sizeof(int) + nElement[0] * (sizeof(double) + sizeof(int));
        break;
    default:
        Abort("Not code!");
        break;
    }
    
    // We assume the stepSize < INT_MAX.
    switch (binFile->revNum) {
    case 0:
        binFile->stepSize = sizeof(int) + sizeof(double) * 6 + sizeof(double) + nElement[0] * sizeof(doubleVector);
        break;
    case 1:
        binFile->stepSize = sizeof(int) + sizeof(uptriMat) + sizeof(double) + nElement[0] * sizeof(doubleVector);
        break;
    case 2:
        binFile->stepSize = sizeof(int) + sizeof(uptriMat) + sizeof(double) * 2 + nElement[0] * sizeof(doubleVector);
        break;
    default:
        Abort("Not code!");
        break;
    }
    binFile->nStep = (int)((binFile->binStat.st_size - binFile->headerSize) / binFile->stepSize);
    
    return binFile;
}
int readSimInfo(Box *box, Particle *particle, Update *update, mmapBinFile *binFile) {
    return readDump(box, particle, update, binFile, 0);
}
int readDump(Box *box, Particle *particle, Update *update, mmapBinFile *binFile, int whichStep) {
    if (whichStep == -1) {
        whichStep = binFile->nStep - 1;
    }
    if (whichStep >= binFile->nStep || whichStep < 0) {
        Abort("Step %d is out of Range: [0,%d];", whichStep, binFile->nStep - 1);
    }
    
    {
        void *data = (void *)(binFile->dataSection);
        if (binFile->revNum != 0) {
            data = data + 32 * sizeof(char) + 2 * sizeof(int);
        }
        particle->nAtom = *((int *)data);
        data = (void *)((char *)data + sizeof(int));
        particle->nAtomType = *((int *)data);
        data = (void *)((char *)data + sizeof(int));
        if (particle->pos == NULL) {
            particle->pos = (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
            particle->veloc = (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
            particle->timeStamp = (double *)calloc(particle->nAtom, sizeof(double));
            particle->mass = (double *)calloc(particle->nAtom, sizeof(double));
            particle->img = (intVector *)calloc(particle->nAtom, sizeof(intVector));
            particle->type = (int *)calloc(particle->nAtom, sizeof(int));
            particle->diameterScale = (double *)calloc(particle->nAtom, sizeof(double));
            particle->id2tag = (int *)calloc(particle->nAtom, sizeof(int));
            particle->tag2id = (int *)calloc(particle->nAtom, sizeof(int));
            particle->massPerType = (double *)calloc(particle->nAtomType, sizeof(double));
            for (int itype = 0; itype < particle->nAtomType; itype++)
                particle->massPerType[itype] = 1.0;
        }
        memcpy(particle->diameterScale, data, particle->nAtom * sizeof(double));
        data = (void *)((char *)data + sizeof(double) * particle->nAtom);
        memcpy(particle->type, data, particle->nAtom * sizeof(int));
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            particle->mass[iatom] = particle->massPerType[particle->type[iatom]];
        }
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            particle->id2tag[iatom] = iatom;
            particle->tag2id[iatom] = iatom;
        }
        
        //====increase step by step=====
        data = (char *)binFile->dataSection + binFile->headerSize;
        for (int istep = 0; istep < whichStep; istep++) {
            data = ((char *)data + binFile->stepSize);
        }
        data = ((char *)data + sizeof(int));
        
        if (binFile->revNum == 0) {
            double *tmp = (double *)data;
            box->boxH[spaceIdx2voigt(0, 0)] = tmp[0];
            box->boxH[spaceIdx2voigt(1, 1)] = tmp[1];
            box->boxH[spaceIdx2voigt(2, 2)] = tmp[2];
            box->boxH[spaceIdx2voigt(1, 2)] = tmp[3];
            box->boxH[spaceIdx2voigt(0, 2)] = tmp[4];
            box->boxH[spaceIdx2voigt(0, 1)] = tmp[5];
            data = ((char *)data + 6 * sizeof(double));
        } else {
            memcpy(box->boxH, data, sizeof(uptriMat));
            data = ((char *)data + sizeof(uptriMat));
        }
        
        particle->meanDiameter = *((double *)data);
        data = ((char *)data + sizeof(double));
        update->runtimeReal = 0.0;
        if (binFile->revNum == 2) {
            double timestep = *((double *)data);
            update->stepCol = (int)floor(timestep);
            update->colCnt = (int)round((timestep - update->stepCol)*particle->nAtom);
            data = ((char *)data + sizeof(double));
        }
        memcpy(particle->pos, data, particle->nAtom * sizeof(doubleVector));
        
        memset(particle->veloc, '\0', particle->nAtom * sizeof(doubleVector));
        memset(particle->timeStamp, '\0', particle->nAtom * sizeof(double));
        memset(particle->img, '\0', particle->nAtom * sizeof(intVector));
    }
    
    particle->isSync = true;
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    update->Edone = update->Zdone = update->Tdone = false;
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.doSort = true;
    update->nebrList.skinSet = __minSkinSet__;
    
    setBoxPara(box);
    initConfInfo(box, particle, update);
    setUnits(update, particle->meanDiameter);
    reInitSim(box, particle, update);
    adjustImg(box, particle);
    
    return 0;
}
int closeBinFile(mmapBinFile **binFilePtr) {
    if (binFilePtr[0] == NULL)
        return 0;
    munmap(binFilePtr[0]->dataSection, binFilePtr[0]->binStat.st_size);
    close(binFilePtr[0]->fd);
    safeFree(binFilePtr[0]);
    return 0;
}

//===SWAP MC move===
SWAP *getSwap(Update *update) {
    int whichTool = findToolkit(&update->toolkit, "__SWAP_MC__");
    if (whichTool < 0)
        return (SWAP *)NULL;
    return (SWAP *)update->toolkit.toolkit[whichTool];
}
SWAP *addSwap(Box *box, Particle *particle, Update *update, Variable *var) {
    cmdArg *cmd = findVariable(var, "swap");
    if (!cmd || (cmd->cmdArgc != 0 && cmd->cmdArgc != 2 && cmd->cmdArgc != 3 && cmd->cmdArgc != 5))
        Abort(
              "--swap: Never stop running by swap;"
              "\n--swap [sisf q x0] [msd y0]: stop running if [sisf(q;t) <= x0] [msd >= y0];");
    
    // set swap parameter
    SWAP *swap = getSwap(update);
    if (!swap) {
        swap = (SWAP *)calloc(sizeof(SWAP), 1);
        addToolkit(&update->toolkit, (void *)swap, NULL, "__SWAP_MC__");
    } else {
        Info("Changing the parameter of existing SWAP.");
    }
    
    swap->max_fs = 1.5;
    swap->min_msd = -1.0;
    if (cmd->cmdArgc == 0) {
        swap->stopBySwap = 0;
    } else {
        for (int ith = 0; ith < cmd->cmdArgc;) {
            if (strcmp(cmd->cmdArgv[ith], "sisf") == 0) {
                if (ith + 3 > cmd->cmdArgc) {
                    Abort(
                          "--swap: Never stop running by swap;"
                          "\n\t\t--swap [sisf q x0] [msd y0]: stop running if [sisf(q;t) <= x0] [msd >= y0];");
                }
                swap->k_0 = atof(cmd->cmdArgv[ith + 1]) / update->distanceUnits;
                swap->max_fs = atof(cmd->cmdArgv[ith + 2]);
                ith += 3;
                if (swap->k_0 <= 1E-6) {
                    Abort("--sisf 7.3 0.2");
                }
                swap->stopBySwap = 1;
            } else if (strcmp(cmd->cmdArgv[ith], "msd") == 0) {
                if (ith + 2 > cmd->cmdArgc) {
                    Abort(
                          "--swap: Never stop running by swap;"
                          "\n\t\t--swap [sisf q x0] [msd y0]: stop running if [sisf(q;t) <= x0] [msd >= y0];");
                }
                swap->min_msd = atof(cmd->cmdArgv[ith + 1]);
                ith += 2;
                swap->stopBySwap = 1;
            } else
                Abort(
                      "--swap: Never stop running by swap;"
                      "\n\t\t--swap [sisf q x0] [msd y0]: stop running if [sisf(q;t) <= x0] [msd >= y0];");
        }
        
        if (swap->xyz0 == NULL)
            swap->xyz0 = (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
        for (int iatom = 0; iatom < particle->nAtom; iatom++) {
            int idx = particle->tag2id[iatom];
            unwrapPos(swap->xyz0[iatom], particle->pos[idx], particle->img[idx],
                      box->boxH);
        }
        
        swap->ntrial = swap->naccept = 0;
    }
    swap->max_dsig = 0.2;
    swap->swap_frac = 0.8;
    
    swap->ntrial = swap->naccept = 0;
    swap->lastStepCol = update->stepCol;
    swap->lastColCnt = update->colCnt;
    return swap;
}
void calcBinParam4swap(Box *box, Particle *particle, Update *update, SWAP *swap) {
    if (box->isShapeFixed && particle->isSizeFixed && swap->deltaAdjBin)
        return;
    
    doubleVector distPlane;
    calcDistBoxPlane(distPlane, box->boxEdge);
    
    double sysRcs = update->nebrList.maxDiameterScale * particle->meanDiameter;
    double maxRcut = update->nebrList.maxDiameterScale * particle->meanDiameter;
    swap->totBin = 1;
    intVector bstart, bstop, ndb;
    int nAdjBin = 1;
    for (int idim = 0; idim < DIM; idim++) {
        int inbin = (int)floor(distPlane[idim] / sysRcs);
        inbin = cpuMax(inbin, 1);
        double ibinLen = distPlane[idim] / inbin;
        int ixyzNum = (int)ceil(maxRcut / ibinLen);
        int ibstart = -ixyzNum;
        int ibstop = ixyzNum;
        int indb = 2 * ixyzNum + 1;
        if (indb > inbin) {
            inbin = 2 * (inbin / 2) + 1;
            ibinLen = distPlane[idim] / inbin;
            ixyzNum = (inbin / 2);
            ibstart = -ixyzNum;
            ibstop = ixyzNum;
            indb = 2 * ixyzNum + 1;
        }
        
        swap->nbin[idim] = inbin;
        swap->binLen[idim] = ibinLen;
        swap->totBin *= inbin;
        
        bstart[idim] = ibstart;
        bstop[idim] = ibstop;
        ndb[idim] = indb;
        nAdjBin *= ndb[idim];
    }
    
    if (nAdjBin > swap->nAdjBin)
        swap->deltaAdjBin =
        (intVector *)realloc(swap->deltaAdjBin, nAdjBin * sizeof(intVector));
    
    swap->nAdjBin = 0;
    for (int idx = 0; idx < nAdjBin; idx++) {
        intVector delta;
        int In = idx;
        for (int idim = 0; idim < DIM; idim++) {
            delta[idim] = In % ndb[idim] + bstart[idim];
            In = In / ndb[idim];
        }
        vAdd(delta, delta, swap->nbin);
        vCpy(swap->deltaAdjBin[swap->nAdjBin], delta);
        swap->nAdjBin++;
    }
    
    if (swap->totBin > swap->maxBinHead) {
        swap->maxBinHead = swap->totBin;
        swap->binHead = (idPosRadius *)realloc(
                                               swap->binHead, swap->maxBinHead * sizeof(idPosRadius));
    }
    if (swap->binIdx == NULL) {
        swap->binIdx = (intVector *)calloc(particle->nAtom, sizeof(intVector));
        swap->binList = (idPosRadius *)calloc(particle->nAtom, sizeof(idPosRadius));
    }
}
void binParticle4swap(Box *box, Particle *particle, Update *update, SWAP *swap) {
    adjustImg(box, particle);
    calcBinParam4swap(box, particle, update, swap);
    
    for (int ibin = 0; ibin < swap->totBin; ibin++) {
        swap->binHead[ibin].id = -1;
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        idPosRadius data;
        data.id = iatom;
        vCpy(data.pos, particle->pos[iatom]);
        data.radius = -1;
        
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, data.pos);
        vShiftAll(lamda, 0.5);
        int binIdx = 0;
        intVector vBinIdx;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            int cIdx = (int)floor(lamda[idim] * swap->nbin[idim]);
            cIdx = (cIdx < 0 ? swap->nbin[idim] - 1 : cIdx);
            cIdx = (cIdx >= swap->nbin[idim] ? 0 : cIdx);
            binIdx = binIdx * swap->nbin[idim] + cIdx;
            vBinIdx[idim] = cIdx;
        }
        vCpy(swap->binIdx[iatom], vBinIdx);
        
        swap->binList[data.id] = swap->binHead[binIdx];
        swap->binHead[binIdx] = data;
    }
}
ReturnType swapMCmove(Box *box, Particle *particle, Update *update) {
    if (!particle->isSync) {
        Abort("Please sync the particles!");
    }
    
    SWAP *swap = getSwap(update);
    
    long int ntrial = (update->stepCol - swap->lastStepCol) * particle->nAtom + (update->colCnt - swap->lastColCnt);
    ntrial = (int)ceil(ntrial * swap->swap_frac);
    swap->lastStepCol = update->stepCol;
    swap->lastColCnt = update->colCnt;
    double meanRadius = particle->meanDiameter * 0.5;
    // double maxRadius = update->nebrList.maxDiameterScale * meanRadius;
    
    binParticle4swap(box, particle, update, swap);
    for (long int itrial = 0; itrial < ntrial; itrial++) {
        swap->ntrial++;
        
        int iatom = ((int)(rndUniform() * particle->nAtom)) % particle->nAtom;
        int jatom = ((int)(rndUniform() * particle->nAtom)) % particle->nAtom;
        while (jatom == iatom) {
            jatom = ((int)(rndUniform() * particle->nAtom)) % particle->nAtom;
        }
        
        int small = iatom, large = jatom;
        if (particle->diameterScale[small] > particle->diameterScale[large]) {
            small = jatom;
            large = iatom;
        }
        
        double sR = particle->diameterScale[small];
        double bR = particle->diameterScale[large];
        if (bR - sR > swap->max_dsig)
            continue;
        if (sR == bR && particle->type[iatom] == particle->type[jatom])
            continue;
        
        doubleVector sPos;
        vCpy(sPos, particle->pos[small]);
        intVector cIdx;
        vCpy(cIdx, swap->binIdx[small]);
        double newRc = bR * meanRadius;
        
        bool accept = true;
        for (int adj = 0; adj < swap->nAdjBin; adj++) {
            intVector binAdjVec;
            for (int idim = 0; idim < DIM; idim++) {
                binAdjVec[idim] = (cIdx[idim] + swap->deltaAdjBin[adj][idim]) % swap->nbin[idim];
            }
            int adjBidx = 0;
            for (int idim = DIM - 1; idim >= 0; idim--) {
                adjBidx = adjBidx * swap->nbin[idim] + binAdjVec[idim];
            }
            
            for (idPosRadius jdata = swap->binHead[adjBidx]; jdata.id >= 0; jdata = swap->binList[jdata.id]) {
                if (jdata.id == iatom || jdata.id == jatom)
                    continue;
                
                doubleVector dRij;
                vSub(dRij, sPos, jdata.pos);
                PBC(dRij, box);
                double rij = sNorm(dRij);
                double sRc = newRc + particle->diameterScale[jdata.id] * meanRadius;
                sRc -= 10 * DBL_EPSILON;
                // reserve enough space
                
                if (rij < sRc) {
                    accept = false;
                    goto __afterCheck__;
                }
            }
        }
        
    __afterCheck__:
        if (!accept)
            continue;
        
        // switch
        particle->diameterScale[small] = bR;
        particle->diameterScale[large] = sR;
        
        int jtype = particle->type[small];
        particle->type[small] = particle->type[large];
        particle->type[large] = jtype;
        
        double jmass = particle->mass[small];
        particle->mass[small] = particle->mass[large];
        particle->mass[large] = jmass;
        
        swap->naccept++;
    }
    update->nebrList.isValid = false;
    update->eventList.isValid = false;
    if (swap->stopBySwap == 0) {
        update->rtype |= TrivialReturn;
        return TrivialReturn;
    }
    
    double fs = 0.0, msd = 0.0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        int idx = particle->tag2id[iatom];
        doubleVector xyz, vR;
        unwrapPos(xyz, particle->pos[idx], particle->img[idx], box->boxH);
        vSub(vR, xyz, swap->xyz0[iatom]);
        
        double dr2 = sNormP2(vR);
        double dr = sNorm(vR);
        double qr = swap->k_0 * dr;
        if (qr <= 1E-12)
            fs += 1.0;
        else
            fs += sin(qr) / qr;
        msd += dr2;
    }
    swap->selfISF = fs / particle->nAtom;
    swap->msd = msd / particle->nAtom / update->distanceUnits / update->distanceUnits;
    
    if ((swap->selfISF <= swap->max_fs) && (swap->msd >= swap->min_msd)) {
        update->rtype |= HaltReturn;
        return HaltReturn;
    } else {
        update->rtype |= TrivialReturn;
        return TrivialReturn;
    }
}
int delSwap(Update *update) {
    SWAP *swap = getSwap(update);
    if (!swap)
        return -1;
    
    double ratio = (double)swap->naccept / (double)swap->ntrial;
    safeFprintf(stderr, "# of trial move: %lld and %lld is accept (%g);\n", swap->ntrial, swap->naccept, ratio);
    
    safeFree(swap->xyz0);
    safeFree(swap->binHead);
    safeFree(swap->binList);
    safeFree(swap->binIdx);
    safeFree(swap->deltaAdjBin);
    delToolkit(&update->toolkit, "__SWAP_MC__");
    return 0;
}

void generateUnwarpPos(Box *box, Particle *particle) {
    if (particle->pos == NULL)
        return;
    
    if (particle->upos == NULL) {
        particle->upos =
        (doubleVector *)calloc(particle->nAtom, sizeof(doubleVector));
    }
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        unwrapPos(particle->upos[iatom], particle->pos[iatom], particle->img[iatom],
                  box->boxH);
    }
}

#if (DIM == 2 || DIM == 3)
void writeLammpsTopo(Box *box, Particle *particle, Update *update, char *fname) {
    //    cmdArg* cmd = findVariable(var, "topo");
    //    if (!cmd)
    //        return;
    //    if (cmd->cmdArgc != 1)
    //        Abort("./app --topo topo.data");
    
    FILE *fout = createFileReadWrite(fname);
    
    safeFprintf(fout,
                "LAMMPS compatible data file. atom_style: sphere (id type "
                "diameter density x y z ix iy iz).\n\n");
    safeFprintf(fout, "\t%d atoms\n", particle->nAtom);
    safeFprintf(fout, "\t%d atom types\n", particle->nAtomType);
    
    double distanceUnits = update->distanceUnits;
    double velocityUnits = update->velocityUnits;
#if (DIM == 3)
    safeFprintf(fout,
                "\n\t%g %g xlo xhi\n\t%g %g ylo yhi\n\t%g %g "
                "zlo zhi\n\t%g %g %g xy xz yz\n",
                box->cornerLo[0] / distanceUnits,
                box->cornerLo[0] / distanceUnits +
                box->boxH[spaceIdx2voigt(0, 0)] / distanceUnits,
                box->cornerLo[1] / distanceUnits,
                box->cornerLo[1] / distanceUnits +
                box->boxH[spaceIdx2voigt(1, 1)] / distanceUnits,
                box->cornerLo[2] / distanceUnits,
                box->cornerLo[2] / distanceUnits +
                box->boxH[spaceIdx2voigt(2, 2)] / distanceUnits,
                box->boxH[spaceIdx2voigt(0, 1)] / distanceUnits,
                box->boxH[spaceIdx2voigt(0, 2)] / distanceUnits,
                box->boxH[spaceIdx2voigt(1, 2)] / distanceUnits);
#elif (DIM == 2)
    safeFprintf(fout,
                "\n\t%g %g xlo xhi\n\t%g %g ylo yhi\n\t0 0 "
                "zlo zhi\n\t%g 0 0 xy xz yz\n",
                box->cornerLo[0] / distanceUnits,
                box->cornerLo[0] / distanceUnits +
                box->boxH[spaceIdx2voigt(0, 0)] / distanceUnits,
                box->cornerLo[1] / distanceUnits,
                box->cornerLo[1] / distanceUnits +
                box->boxH[spaceIdx2voigt(1, 1)] / distanceUnits,
                box->boxH[spaceIdx2voigt(0, 1)] / distanceUnits);
#endif
    
    safeFprintf(fout, "\nAtoms\n\n");
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        int idx = particle->tag2id[iatom];
        
        int type = particle->type[idx];
        double diameter =
        particle->diameterScale[idx] * particle->meanDiameter / distanceUnits;
        doubleVecPtr posPtr = particle->pos[idx];
        intVecPtr imgPtr = particle->img[idx];
        
#if (DIM == 3)
        safeFprintf(fout, "%d %d %g 1 %g %g %g %d %d %d\n", iatom + 1, type + 1,
                    diameter, posPtr[0] / distanceUnits, posPtr[1] / distanceUnits,
                    posPtr[2] / distanceUnits, imgPtr[0], imgPtr[1], imgPtr[2]);
#elif (DIM == 2)
        safeFprintf(fout, "%d %d %g 1 %g %g 0 %d %d 0\n", iatom + 1, type + 1,
                    diameter, posPtr[0] / distanceUnits, posPtr[1] / distanceUnits,
                    imgPtr[0], imgPtr[1]);
#endif
    }
    
    safeFprintf(fout, "\nVelocities\n\n");
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        int idx = particle->tag2id[iatom];
        doubleVecPtr velocPtr = particle->veloc[idx];
        
#if (DIM == 3)
        safeFprintf(fout, "%d %g %g %g\n", iatom + 1, velocPtr[0] / velocityUnits,
                    velocPtr[1] / velocityUnits, velocPtr[2] / velocityUnits);
#elif (DIM == 2)
        safeFprintf(fout, "%d %g %g %g\n", iatom + 1, velocPtr[0] / velocityUnits,
                    velocPtr[1] / velocityUnits, 0.0);
#endif
    }
    
    safeCloseFile(fout);
}
#endif
