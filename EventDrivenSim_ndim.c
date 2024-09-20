#include "EventDrivenSim_ndim.h"

//=======Event List, Nebr List, do collision=======
void sortParticle(Box* box, Particle* particle, Update* update) {
    if (!update->nebrList.doSort)
        return;
    if (particle->__isSortForbidden) {
        return;
    }
    
    particle->sortFlag++;
    
    adjustImg(box, particle);
    
    NebrList* nebrList = &update->nebrList;
    doubleVector distPlane;
    calcDistBoxPlane(distPlane, box->boxEdge);
    
    double sysRcs = particle->meanDiameter * (1.0 + nebrList->skinSet);
    nebrList->totBin4sort = 1;
    intVector nbin4sort;
    for (int idim = 0; idim < DIM; idim++) {
        nbin4sort[idim] = (int)floor(distPlane[idim] / sysRcs);
        nbin4sort[idim] = cpuMax(nbin4sort[idim], 1);
        nebrList->totBin4sort *= nbin4sort[idim];
    }
    if (nebrList->totBin4sort > nebrList->allocBin4sort) {
        nebrList->binHead4sort = (int*)realloc(nebrList->binHead4sort,
                                               nebrList->totBin4sort * sizeof(int));
        nebrList->allocBin4sort = nebrList->totBin4sort;
    }
    
    for (int ibin = 0; ibin < nebrList->totBin4sort; ibin++) {
        nebrList->binHead4sort[ibin] = -1;
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVecPtr pos = particle->pos[iatom];
        
        intVector cIdx;
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, pos);
        vShiftAll(lamda, 0.5);
        int binIdx = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            cIdx[idim] = (int)floor(lamda[idim] * nbin4sort[idim]);
            cIdx[idim] = (cIdx[idim] < 0 ? nbin4sort[idim] - 1 : cIdx[idim]);
            cIdx[idim] = (cIdx[idim] >= nbin4sort[idim] ? 0 : cIdx[idim]);
            binIdx = binIdx * nbin4sort[idim] + cIdx[idim];
        }
        nebrList->binList4sort[iatom] = nebrList->binHead4sort[binIdx];
        nebrList->binHead4sort[binIdx] = iatom;
    }
    
    for (int idx = 0, noid = 0; idx < nebrList->totBin4sort; idx++) {
        int oid = nebrList->binHead4sort[idx];
        while (oid != -1) {
            nebrList->oid2nid[oid] = noid;
            noid++;
            oid = nebrList->binList4sort[oid];
        }
    }
    
    exchange_doubleVector(particle->pos, (doubleVector*)nebrList->buffer,
                          nebrList->oid2nid, particle->nAtom);
    memcpy(particle->pos, nebrList->buffer,
           particle->nAtom * sizeof(doubleVector));
    
    exchange_doubleVector(particle->veloc, (doubleVector*)nebrList->buffer,
                          nebrList->oid2nid, particle->nAtom);
    memcpy(particle->veloc, nebrList->buffer,
           particle->nAtom * sizeof(doubleVector));
    
    exchange_intVector(particle->img, (intVector*)nebrList->buffer,
                       nebrList->oid2nid, particle->nAtom);
    memcpy(particle->img, nebrList->buffer, particle->nAtom * sizeof(intVector));
    
    exchange_double(particle->mass, (double*)nebrList->buffer, nebrList->oid2nid,
                    particle->nAtom);
    memcpy(particle->mass, nebrList->buffer, particle->nAtom * sizeof(double));
    
    exchange_double(particle->diameterScale, (double*)nebrList->buffer,
                    nebrList->oid2nid, particle->nAtom);
    memcpy(particle->diameterScale, nebrList->buffer,
           particle->nAtom * sizeof(double));
    
    exchange_int(particle->type, (int*)nebrList->buffer, nebrList->oid2nid,
                 particle->nAtom);
    memcpy(particle->type, nebrList->buffer, particle->nAtom * sizeof(int));
    
    exchange_int(particle->id2tag, (int*)nebrList->buffer, nebrList->oid2nid,
                 particle->nAtom);
    memcpy(particle->id2tag, nebrList->buffer, particle->nAtom * sizeof(int));
    
    for (int aid = 0; aid < particle->nAtom; aid++) {
        particle->tag2id[particle->id2tag[aid]] = aid;
    }
    
    nebrList->isValid = false;
    update->eventList.isValid = false;
}
//=======================================
//===build event list===
double QuadraticFormula(double A, double B, double C) {
    // solve Quadratic equation: A * x^2 + 2 * B * x + C = 0;
    double xroot = DBL_INF;
    double det = B * B - A * C;
    
    if (C <= 0.0) {
        if (B < 0.0)
            xroot = 0.0;  // spheres already overlapping and approaching
        else if (A < 0.0) {
            if (det <= 0.0) {
                Abort("Fatal Error!");
            } else {
                xroot = (-B - sqrt(det)) / A;
            }
        }
        return xroot;
    }
    
    if (det > -10.0 * DBL_EPSILON) {
        det = (det < 0.0 ? 0.0 : det);
        if (B < 0.0) {
            xroot = C / (-B + sqrt(det));
        } else if (A < 0.0) {
            xroot = (-B - sqrt(det)) / A;
        }
    }
    
    return xroot;
}

//===do collision===
void doColHeapTop_VR(Box* box, Particle* particle, Update* update) {
    Event data = update->eventList.data[update->eventList.eIdx[1]];
    int iatom = data.atomIdx;
    int jatom = data.ePartner;
    double eTime = data.eTime;
    
    // sumation of diameter at timestamp 0
    double rrate = update->rrateSet / update->timeUnits;
    double sDiameterTime0 = particle->meanDiameter * (particle->diameterScale[iatom] + particle->diameterScale[jatom]);
    double sGrowthRate = sDiameterTime0 * rrate * 0.5;  //(rrate > 0 ? sDiameterTime0 * rrate * 0.5 : 0);
    double sRc = 0.5 * sDiameterTime0 * (1.0 + rrate * eTime);
    
    // get deltaPij and check
    doubleVector dVij;
    vSub(dVij, particle->veloc[iatom], particle->veloc[jatom]);
    
    // update position
    vScaleAdd(particle->pos[iatom], particle->pos[iatom],
              eTime - particle->timeStamp[iatom], particle->veloc[iatom]);
    vScaleAdd(particle->pos[jatom], particle->pos[jatom],
              eTime - particle->timeStamp[jatom], particle->veloc[jatom]);
    
    // calculate unit vector
    doubleVector vRij, Nij;
    vSub(vRij, particle->pos[iatom], particle->pos[jatom]);
    PBC(vRij, box);
    double rij = sNorm(vRij);
    vScale(Nij, 1.0 / rij, vRij);
    
    // get deltaPij and check
    double vc = sDot(dVij, Nij);
    double deltaPij = -vc + sGrowthRate;
    if ((sRc - rij < -10.0 * DBL_EPSILON) || (sRc - rij > 10.0 * DBL_EPSILON)) {
        update->errCnt++;
        if (update->errCnt == particle->nAtom) {
            update->errCnt = 0;
            update->stepErr++;
        }
    }
    
#if 0
    //thermostat at particle level.
    double w1 = sDot(Nij, particle->veloc[iatom]), w2 = sDot(Nij, particle->veloc[jatom]);
    double Eperp = sNormP2(particle->veloc[iatom]) + sNormP2(particle->veloc[jatom]) - w1 * w1 - w2 * w2;
    Eperp = (Eperp < 1E-20 ? 1E-20 : Eperp);
    double sfactPerp = 2.0 * deltaPij * deltaPij - 2.0 * (w2 - w1) * deltaPij;
    sfactPerp = 1.0 - sfactPerp / Eperp;
    sfactPerp = (sfactPerp < 0 ? 0 : sfactPerp);
    sfactPerp = (sfactPerp > 1.0 ? 1.0 : sfactPerp);
    sfactPerp = sqrt(sfactPerp);
    
    vScaleAdd(particle->veloc[iatom], particle->veloc[iatom], -w1, Nij);
    vScale(particle->veloc[iatom], sfactPerp, particle->veloc[iatom]);
    vScaleAdd(particle->veloc[iatom], particle->veloc[iatom], w1, Nij);
    vScaleAdd(particle->veloc[iatom], particle->veloc[iatom], deltaPij, Nij);
    
    vScaleAdd(particle->veloc[jatom], particle->veloc[jatom], -w2, Nij);
    vScale(particle->veloc[jatom], sfactPerp, particle->veloc[jatom]);
    vScaleAdd(particle->veloc[jatom], particle->veloc[jatom], w2, Nij);
    vScaleAdd(particle->veloc[jatom], particle->veloc[jatom], -deltaPij, Nij);
#endif
    
    // get velocity
    vScaleAdd(particle->veloc[iatom], particle->veloc[iatom], deltaPij, Nij);
    vScaleAdd(particle->veloc[jatom], particle->veloc[jatom], -deltaPij, Nij);
    
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            // update->sKinTensorMulTime[spaceIdx2voigt(idim, jdim)] +=
            //     update->sColDeltaKinTensor[spaceIdx2voigt(idim, idim)] *
            //     (update->current_timeStamp - update->lastEventTimeStamp);
            // update->sColDeltaKinTensor[spaceIdx2voigt(idim, jdim)] +=
            //     deltaPij * Nij[idim] * dVij[jdim] +
            //     deltaPij * Nij[jdim] * dVij[idim] +
            //     2.0 * deltaPij * deltaPij * Nij[idim] * Nij[jdim];
            update->sColVirialTensor[spaceIdx2voigt(idim, jdim)] += deltaPij * vRij[idim] * vRij[jdim] / rij;
        }
    }
    particle->timeStamp[iatom] = particle->timeStamp[jatom] = eTime;
    
    update->colCnt++;
    if (update->colCnt == particle->nAtom) {
        update->colCnt = 0;
        update->stepCol++;
    }
    update->nebrList.accColRebuild++;
    update->nebrList.colRebuild++;
}
void doColHeapTop_FR(Box* box, Particle* particle, Update* update) {
    Event data = update->eventList.data[update->eventList.eIdx[1]];
    int iatom = data.atomIdx;
    int jatom = data.ePartner;
    double eTime = data.eTime;
    
    double sRc = 0.5 * (particle->diameterScale[iatom] + particle->diameterScale[jatom]) * particle->meanDiameter;
    
    // get deltaPij and check
    doubleVector dVij;
    vSub(dVij, particle->veloc[iatom], particle->veloc[jatom]);
    
    // update position
    vScaleAdd(particle->pos[iatom], particle->pos[iatom],
              eTime - particle->timeStamp[iatom], particle->veloc[iatom]);
    vScaleAdd(particle->pos[jatom], particle->pos[jatom],
              eTime - particle->timeStamp[jatom], particle->veloc[jatom]);
    
    // calculate unit vector
    doubleVector vRij, Nij;
    vSub(vRij, particle->pos[iatom], particle->pos[jatom]);
    PBC(vRij, box);
    double rij = sNorm(vRij);
    vScale(Nij, 1.0 / rij, vRij);
    
    // get deltaPij and check
    double vc = sDot(dVij, Nij);
    double deltaPij = -vc;
    if ((sRc - rij < -10.0 * DBL_EPSILON) || (sRc - rij > 10.0 * DBL_EPSILON)) {
        update->errCnt++;
        if (update->errCnt == particle->nAtom) {
            update->stepErr++;
            update->errCnt = 0;
        }
    }
    
    // get velocity
    vScaleAdd(particle->veloc[iatom], particle->veloc[iatom], deltaPij, Nij);
    vScaleAdd(particle->veloc[jatom], particle->veloc[jatom], -deltaPij, Nij);
    
    for (int idim = 0; idim < DIM; idim++) {
        for (int jdim = idim; jdim < DIM; jdim++) {
            // update->sKinTensorMulTime[spaceIdx2voigt(idim, jdim)] +=
            //     update->sColDeltaKinTensor[spaceIdx2voigt(idim, idim)] *
            //     (update->current_timeStamp - update->lastEventTimeStamp);
            // update->sColDeltaKinTensor[spaceIdx2voigt(idim, jdim)] +=
            //     deltaPij * Nij[idim] * dVij[jdim] +
            //     deltaPij * Nij[jdim] * dVij[idim] +
            //     2.0 * deltaPij * deltaPij * Nij[idim] * Nij[jdim];
            update->sColVirialTensor[spaceIdx2voigt(idim, jdim)] +=
            deltaPij * vRij[idim] * vRij[jdim] / rij;
        }
    }
    particle->timeStamp[iatom] = particle->timeStamp[jatom] = eTime;
    
    update->colCnt++;
    if (update->colCnt == particle->nAtom) {
        update->stepCol++;
        update->colCnt = 0;
    }
    update->nebrList.accColRebuild++;
    update->nebrList.colRebuild++;
}

//===build event list===
void calcBinParam(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    if (!nebrList->compelInit && box->isShapeFixed && particle->isSizeFixed && nebrList->deltaAdjBin)
        return;
    
    doubleVector distPlane;
    calcDistBoxPlane(distPlane, box->boxEdge);
    
    double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
    for (int idim = 0; idim < DIM; idim++) {
        minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
    }
    double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
    // if {0.5*maxAxByCz >= maxRcut + 2*Rskin} then PBC method is valid before
    // collision
    double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
    nebrList->maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
    if (nebrList->maxRskinSet < __minSkinSet__) {
        double minRskin = __minSkinSet__ * update->nebrList.minDiameterScale * particle->meanDiameter;
        double minLen = (maxRcut + minRskin) * 2.0;
        double minBoxVol = pow(minLen, DIM);
        double maxPhi = update->volFrac * box->volume / minBoxVol;
        
        Abort(
              "The box is too small, the PBC method breaks down! Please use large N, "
              "or rewrite codes! The 0.5*min(box(.,.)) is %g, max(Rcut) is %g, "
              "min(Rskin) is %g. The max(Rcut)+min(Rskin) should be smaller "
              "than 0.5*min(ax,by,cz). The estimated max(volFrac) is %g, currently, "
              "it is %g.",
              0.5 * minAxByCz / update->distanceUnits,
              maxRcut / update->distanceUnits, minRskin / update->distanceUnits,
              maxPhi, update->volFrac);
    }
    
    nebrList->rskin = nebrList->minDiameterScale * particle->meanDiameter * nebrList->skinSet;
    nebrList->rskin = cpuMin(nebrList->rskin, maxRskin);
    double sysRcs = maxRcut + nebrList->rskin;
    
    nebrList->totBin = 1;
    intVector bstart, bstop, ndb;
    int nAdjBin = 1;
    for (int idim = 0; idim < DIM; idim++) {
        int inbin = (int)floor(distPlane[idim] / sysRcs);
        inbin = cpuMax(inbin, 1);
        double ibinLen = distPlane[idim] / inbin;
        int ixyzNum = (int)ceil((maxRcut + nebrList->rskin) / ibinLen);
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
        
        nebrList->nbin[idim] = inbin;
        nebrList->binLen[idim] = ibinLen;
        nebrList->totBin *= inbin;
        
        bstart[idim] = ibstart;
        bstop[idim] = ibstop;
        ndb[idim] = indb;
        nAdjBin *= ndb[idim];
    }
    
    if (nAdjBin > nebrList->nAdjBin)
        nebrList->deltaAdjBin = (intVector*)realloc(nebrList->deltaAdjBin, nAdjBin * sizeof(intVector));
    
    nebrList->nAdjBin = 0;
    for (int idx = 0; idx < nAdjBin; idx++) {
        intVector delta;
        int In = idx;
        for (int idim = 0; idim < DIM; idim++) {
            delta[idim] = (In % ndb[idim]) + bstart[idim];
            In = In / ndb[idim];
        }
        In = vIdx2lIdx(delta, ndb);
        if (In < 0)
            continue;
        
        vAdd(delta, delta, nebrList->nbin);
        vCpy(nebrList->deltaAdjBin[nebrList->nAdjBin], delta);
        nebrList->nAdjBin++;
    }
    
    if (nebrList->totBin > nebrList->maxBinHead) {
        nebrList->maxBinHead = nebrList->totBin;
        nebrList->binHead = (idPosRadius*)realloc(nebrList->binHead, nebrList->maxBinHead * sizeof(idPosRadius));
    }
    
    nebrList->compelInit = false;
}
void binParticle(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    
    calcBinParam(box, particle, update);
    
    for (int ibin = 0; ibin < nebrList->totBin; ibin++) {
        nebrList->binHead[ibin].id = -1;
    }
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        idPosRadius data;
        data.id = iatom;
        vCpy(data.pos, particle->pos[iatom]);
        vCpy(data.veloc, particle->veloc[iatom]);
        data.radius = 0.5 * particle->meanDiameter * particle->diameterScale[iatom];
        
        doubleVector lamda;
        MatMulVec(lamda, box->invBoxH, data.pos);
        vShiftAll(lamda, 0.5);
        int binIdx = 0;
        for (int idim = DIM - 1; idim >= 0; idim--) {
            int cIdx = (int)floor(lamda[idim] * nebrList->nbin[idim]);
            cIdx = (cIdx < 0 ? nebrList->nbin[idim] - 1 : cIdx);
            cIdx = (cIdx >= nebrList->nbin[idim] ? 0 : cIdx);
            binIdx = binIdx * nebrList->nbin[idim] + cIdx;
        }
        
        nebrList->binList[data.id] = nebrList->binHead[binIdx];
        nebrList->binHead[binIdx] = data;
    }
}

//=================
#define __accFact__ (1.0 / 10.0)
void constructEventNebrList_FR(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    MinEventHeap* eventList = &update->eventList;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector iveloc;
        vCpy(iveloc, particle->veloc[iatom]);
        double A = -sDot(iveloc, iveloc);
        double C = nebrList->rskin * 0.5 * nebrList->rskin * 0.5;
        double tcolPredict = QuadraticFormula(A, 0.0, C);
        nebrList->tnebr[iatom] = tcolPredict;
        
        eventList->data[iatom].atomIdx = iatom;
        eventList->data[iatom].ePartner = nebrPartner;
        eventList->data[iatom].eTime = tcolPredict;
        eventList->data[iatom].eventType = NebrListType;
    }
    
    bool overFlow = false;
    do {
        overFlow = false;
        int maxNebrPerAtom = nebrList->maxNebrPerAtom;
        int nColTable = 0;
        memset(nebrList->nNebr, '\0', particle->nAtom * sizeof(int));
        clearMinEventHeap(eventList, particle->nAtom);
        
#if (DIM == 3)
        for (int ibinz = 0; ibinz < nebrList->nbin[2]; ibinz++) {
            for (int ibiny = 0; ibiny < nebrList->nbin[1]; ibiny++) {
                for (int ibinx = 0; ibinx < nebrList->nbin[0]; ibinx++) {
                    int ibin =
                    (ibinz * nebrList->nbin[1] + ibiny) * nebrList->nbin[0] + ibinx;
                    for (idPosRadius idata = nebrList->binHead[ibin]; idata.id >= 0;
                         idata = nebrList->binList[idata.id]) {
                        double iRadiusRskin = idata.radius + nebrList->rskin;
                        Event data = eventList->data[idata.id];
                        
                        for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                            int jbinz =
                            (ibinz + nebrList->deltaAdjBin[adj][2]) % nebrList->nbin[2];
                            int jbiny =
                            (ibiny + nebrList->deltaAdjBin[adj][1]) % nebrList->nbin[1];
                            int jbinx =
                            (ibinx + nebrList->deltaAdjBin[adj][0]) % nebrList->nbin[0];
                            int jbin =
                            (jbinz * nebrList->nbin[1] + jbiny) * nebrList->nbin[0] +
                            jbinx;
                            for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                                 jdata = nebrList->binList[jdata.id]) {
                                if (jbin == ibin && jdata.id <= idata.id)
                                    continue;
                                
                                doubleVector dRij;
                                vSub(dRij, idata.pos, jdata.pos);
                                PBC(dRij, box);
                                double rijP2 = sNormP2(dRij);
                                double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                                if (rijP2 > RcutRskinP2)
                                    continue;
                                
                                if (nebrList->nNebr[idata.id] < nebrList->maxNebrPerAtom &&
                                    nebrList->nNebr[jdata.id] < nebrList->maxNebrPerAtom) {
                                    double sRc = idata.radius + jdata.radius;
                                    doubleVector dVij;
                                    vSub(dVij, idata.veloc, jdata.veloc);
                                    double A = sDot(dVij, dVij);
                                    double B = sDot(dRij, dVij);
                                    double C = sDot(dRij, dRij) - sRc * sRc;
                                    double tcolPredict = QuadraticFormula(A, B, C);
                                    if (nColTable >= nebrList->allocColTable) {
                                        nebrList->allocColTable += particle->nAtom * 8;
                                        nebrList->colPairTable = (double*)realloc(
                                                                                  nebrList->colPairTable,
                                                                                  nebrList->allocColTable * sizeof(double));
                                    }
                                    nebrList->colPairTable[nColTable] = tcolPredict;
                                    
                                    int2 info;
                                    info.first = jdata.id;
                                    info.second = nColTable;
                                    nebrList->list[idata.id * nebrList->maxNebrPerAtom +
                                                   nebrList->nNebr[idata.id]] = info;
                                    info.first = idata.id;
                                    info.second = nColTable;
                                    nebrList->list[jdata.id * nebrList->maxNebrPerAtom +
                                                   nebrList->nNebr[jdata.id]] = info;
                                    
                                    if (tcolPredict <= eventList->data[jdata.id].eTime) {
                                        eventList->data[jdata.id].eTime = tcolPredict;
                                        eventList->data[jdata.id].ePartner = idata.id;
                                        eventList->data[jdata.id].eventType = HardCoreType;
                                    }
                                    if (tcolPredict <= data.eTime) {
                                        data.eTime = tcolPredict;
                                        data.ePartner = jdata.id;
                                        data.eventType = HardCoreType;
                                    }
                                    
                                    nebrList->nNebr[idata.id]++;
                                    nebrList->nNebr[jdata.id]++;
                                    nColTable++;
                                } else {
                                    nebrList->nNebr[idata.id]++;
                                    nebrList->nNebr[jdata.id]++;
                                    maxNebrPerAtom =
                                    cpuMax(nebrList->nNebr[idata.id], maxNebrPerAtom);
                                    maxNebrPerAtom =
                                    cpuMax(nebrList->nNebr[jdata.id], maxNebrPerAtom);
                                    overFlow = true;
                                }
                            }  // end jbin
                            
                        }  // end adjacent bin
                        
                        eventList->data[idata.id] = data;
                    }  // end ibin
                }
            }
        }
#elif (DIM == 2)
        for (int ibiny = 0; ibiny < nebrList->nbin[1]; ibiny++) {
            for (int ibinx = 0; ibinx < nebrList->nbin[0]; ibinx++) {
                int ibin = ibiny * nebrList->nbin[0] + ibinx;
                for (idPosRadius idata = nebrList->binHead[ibin]; idata.id >= 0;
                     idata = nebrList->binList[idata.id]) {
                    double iRadiusRskin = idata.radius + nebrList->rskin;
                    Event data = eventList->data[idata.id];
                    
                    for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                        int jbiny =
                        (ibiny + nebrList->deltaAdjBin[adj][1]) % nebrList->nbin[1];
                        int jbinx =
                        (ibinx + nebrList->deltaAdjBin[adj][0]) % nebrList->nbin[0];
                        int jbin = jbiny * nebrList->nbin[0] + jbinx;
                        for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                             jdata = nebrList->binList[jdata.id]) {
                            if (jbin == ibin && jdata.id <= idata.id)
                                continue;
                            
                            doubleVector dRij;
                            vSub(dRij, idata.pos, jdata.pos);
                            PBC(dRij, box);
                            double rijP2 = sNormP2(dRij);
                            double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                            if (rijP2 > RcutRskinP2)
                                continue;
                            
                            if (nebrList->nNebr[idata.id] < nebrList->maxNebrPerAtom &&
                                nebrList->nNebr[jdata.id] < nebrList->maxNebrPerAtom) {
                                double sRc = idata.radius + jdata.radius;
                                doubleVector dVij;
                                vSub(dVij, idata.veloc, jdata.veloc);
                                double A = sDot(dVij, dVij);
                                double B = sDot(dRij, dVij);
                                double C = sDot(dRij, dRij) - sRc * sRc;
                                double tcolPredict = QuadraticFormula(A, B, C);
                                if (nColTable >= nebrList->allocColTable) {
                                    nebrList->allocColTable += particle->nAtom * 8;
                                    nebrList->colPairTable = (double*)realloc(
                                                                              nebrList->colPairTable,
                                                                              nebrList->allocColTable * sizeof(double));
                                }
                                nebrList->colPairTable[nColTable] = tcolPredict;
                                
                                int2 info;
                                info.first = jdata.id;
                                info.second = nColTable;
                                nebrList->list[idata.id * nebrList->maxNebrPerAtom +
                                               nebrList->nNebr[idata.id]] = info;
                                info.first = idata.id;
                                info.second = nColTable;
                                nebrList->list[jdata.id * nebrList->maxNebrPerAtom +
                                               nebrList->nNebr[jdata.id]] = info;
                                
                                if (tcolPredict <= eventList->data[jdata.id].eTime) {
                                    eventList->data[jdata.id].eTime = tcolPredict;
                                    eventList->data[jdata.id].ePartner = idata.id;
                                    eventList->data[jdata.id].eventType = HardCoreType;
                                }
                                if (tcolPredict <= data.eTime) {
                                    data.eTime = tcolPredict;
                                    data.ePartner = jdata.id;
                                    data.eventType = HardCoreType;
                                }
                                
                                nebrList->nNebr[idata.id]++;
                                nebrList->nNebr[jdata.id]++;
                                nColTable++;
                            } else {
                                nebrList->nNebr[idata.id]++;
                                nebrList->nNebr[jdata.id]++;
                                maxNebrPerAtom =
                                cpuMax(nebrList->nNebr[idata.id], maxNebrPerAtom);
                                maxNebrPerAtom =
                                cpuMax(nebrList->nNebr[jdata.id], maxNebrPerAtom);
                                overFlow = true;
                            }
                        }  // end jbin
                        
                    }  // end adjacent bin
                    
                    eventList->data[idata.id] = data;
                }  // end ibin
            }
        }
#else
        for (int ibin = 0; ibin < nebrList->totBin; ibin++) {
            intVector vibin;
            lidx2vidx(ibin, nebrList->nbin, vibin);
            for (idPosRadius idata = nebrList->binHead[ibin]; idata.id >= 0;
                 idata = nebrList->binList[idata.id]) {
                double iRadiusRskin = idata.radius + nebrList->rskin;
                Event data = eventList->data[idata.id];
                
                for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                    intVector vjbin;
                    vAdd(vjbin, vibin, nebrList->deltaAdjBin[adj]);
                    for (int idim = 0; idim < DIM; idim++) {
                        vjbin[idim] = vjbin[idim] % nebrList->nbin[idim];
                    }
                    int jbin = vIdx2lIdx(vjbin, nebrList->nbin);
                    
                    for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                         jdata = nebrList->binList[jdata.id]) {
                        if (jbin == ibin && jdata.id <= idata.id)
                            continue;
                        
                        doubleVector dRij;
                        vSub(dRij, idata.pos, jdata.pos);
                        PBC(dRij, box);
                        double rijP2 = sNormP2(dRij);
                        double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                        if (rijP2 > RcutRskinP2)
                            continue;
                        
                        if (nebrList->nNebr[idata.id] < nebrList->maxNebrPerAtom &&
                            nebrList->nNebr[jdata.id] < nebrList->maxNebrPerAtom) {
                            double sRc = idata.radius + jdata.radius;
                            doubleVector dVij;
                            vSub(dVij, idata.veloc, jdata.veloc);
                            double A = sDot(dVij, dVij);
                            double B = sDot(dRij, dVij);
                            double C = sDot(dRij, dRij) - sRc * sRc;
                            double tcolPredict = QuadraticFormula(A, B, C);
                            if (nColTable >= nebrList->allocColTable) {
                                nebrList->allocColTable += particle->nAtom * 8;
                                nebrList->colPairTable =
                                (double*)realloc(nebrList->colPairTable,
                                                 nebrList->allocColTable * sizeof(double));
                            }
                            nebrList->colPairTable[nColTable] = tcolPredict;
                            
                            int2 info;
                            info.first = jdata.id;
                            info.second = nColTable;
                            nebrList->list[idata.id * nebrList->maxNebrPerAtom + nebrList->nNebr[idata.id]] = info;
                            info.first = idata.id;
                            info.second = nColTable;
                            nebrList->list[jdata.id * nebrList->maxNebrPerAtom + nebrList->nNebr[jdata.id]] = info;
                            
                            if (tcolPredict <= eventList->data[jdata.id].eTime) {
                                eventList->data[jdata.id].eTime = tcolPredict;
                                eventList->data[jdata.id].ePartner = idata.id;
                                eventList->data[jdata.id].eventType = HardCoreType;
                            }
                            if (tcolPredict <= data.eTime) {
                                data.eTime = tcolPredict;
                                data.ePartner = jdata.id;
                                data.eventType = HardCoreType;
                            }
                            
                            nebrList->nNebr[idata.id]++;
                            nebrList->nNebr[jdata.id]++;
                            nColTable++;
                        } else {
                            nebrList->nNebr[idata.id]++;
                            nebrList->nNebr[jdata.id]++;
                            maxNebrPerAtom = cpuMax(nebrList->nNebr[idata.id], maxNebrPerAtom);
                            maxNebrPerAtom = cpuMax(nebrList->nNebr[jdata.id], maxNebrPerAtom);
                            overFlow = true;
                        }
                    }  // end jbin
                    
                }  // end adjacent bin
                
                eventList->data[idata.id] = data;
            }  // end ibin
        }
#endif
        
        if (overFlow) {
            nebrList->list = (int2*)realloc(nebrList->list, maxNebrPerAtom * particle->nAtom * sizeof(int2));
            nebrList->maxNebrPerAtom = maxNebrPerAtom;
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                eventList->data[iatom].atomIdx = iatom;
                eventList->data[iatom].ePartner = nebrPartner;
                eventList->data[iatom].eTime = nebrList->tnebr[iatom];
                eventList->data[iatom].eventType = NebrListType;
            }
        }
        
    } while (overFlow);
}
void buildEventList_FR(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    if (nebrList->xyzHold == NULL) {
        nebrList->xyzHold =
        (doubleVector*)calloc(particle->nAtom, sizeof(doubleVector));
        nebrList->tnebr = (double*)calloc(particle->nAtom, sizeof(double));
        nebrList->nNebr = (int*)calloc(particle->nAtom, sizeof(int));
        nebrList->binList =
        (idPosRadius*)calloc(particle->nAtom, sizeof(idPosRadius));
        
        nebrList->binList4sort = (int*)calloc(particle->nAtom, sizeof(int));
        nebrList->oid2nid = (int*)calloc(particle->nAtom, sizeof(int));
        nebrList->buffer = calloc(particle->nAtom, sizeof(doubleVector));
    }
    
    if (update->stepCol >= update->nextTtstep)
        thermostat(box, particle, update);
    if (update->stepCol >= update->nextZCstep)
        moveBarycenter(box, particle, update);
    if (update->stepCol >= update->nextZMstep)
        zeroMomentum(box, particle, update);
    
    MinEventHeap* eventList = &update->eventList;
    if (eventList->isValid)
        return;
    clearMinEventHeap(eventList, particle->nAtom);
    
    syncAll(box, particle, update);
    
    if (getSwap(update)) {
        swapMCmove(box, particle, update);
    }
    
    // adjust parameter
    if (nebrList->nBuild == 10 * AveNrebuild) {
        nebrList->accColRebuild = 0;
    }
    if (nebrList->nBuild > 10 * AveNrebuild && nebrList->nBuild % AveNrebuild == 0) {
        if (nebrList->accColRebuild > 0.1 * AveNrebuild * colPerAtomBetweenRebuild * (double)particle->nAtom)
            nebrList->skinSet = nebrList->skinSet * (AveNrebuild * colPerAtomBetweenRebuild * (double)particle->nAtom / (double)nebrList->accColRebuild);
        else {
            nebrList->skinSet *= 1.2;
        }
        
        // estimated from spaces between particles
        double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
        double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
        double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
        // estimated from box side length
        doubleVector distPlane;
        calcDistBoxPlane(distPlane, box->boxEdge);
        double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
        for (int idim = 0; idim < DIM; idim++) {
            minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
        }
        double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
        double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
        double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
        
        double lmaxset = cpuMin(lmsset, maxRskinSet);
        nebrList->skinSet = cpuMin(nebrList->skinSet, lmaxset);
        
        nebrList->compelInit = true;
        nebrList->accColRebuild = 0;
    }
    if ((nebrList->nBuild % 100 == 0)) {
        sortParticle(box, particle, update);
    }
    
    adjustImg(box, particle);
    binParticle(box, particle, update);
    
    memcpy(nebrList->xyzHold, particle->pos,
           particle->nAtom * sizeof(doubleVector));
    
    constructEventNebrList_FR(box, particle, update);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        Event data = eventList->data[iatom];
        replaceMinEventHeap(eventList, &data);
    }
    
    eventList->isValid = true;
    nebrList->isValid = true;
    nebrList->nBuild++;
    nebrList->colRebuild = 0;
}
void updateColPair_FR(Box* box, Particle* particle, Update* update, int eIatom, int ePartner) {
    NebrList* nebrList = &update->nebrList;
    MinEventHeap* eventList = &update->eventList;
    double meanRadius = 0.5 * particle->meanDiameter;
    
    {
        int iatom = eIatom;
        Event data;
        double iRc = particle->diameterScale[iatom] * meanRadius;
        double istamp = particle->timeStamp[iatom];
        doubleVector ipos, iveloc;
        vCpy(ipos, particle->pos[iatom]);
        vCpy(iveloc, particle->veloc[iatom]);
        // Qij
        {
            doubleVector Rt0;
            vSub(Rt0, ipos, nebrList->xyzHold[iatom]);
            doubleVector dVel;
            vCpy(dVel, iveloc);
            double rst = nebrList->rskin * 0.5;
            
            double A = -sDot(dVel, dVel);
            double B = -sDot(Rt0, dVel);
            double C = -sDot(Rt0, Rt0) + rst * rst;
            
            data.eTime = istamp + (C <= 0.0 ? 0.0 : QuadraticFormula(A, B, C));
            data.ePartner = nebrPartner;
            data.atomIdx = iatom;
            data.eventType = NebrListType;
            
            nebrList->tnebr[iatom] = data.eTime;
        }
        
        // Pij
        for (int cnt = 0; cnt < nebrList->nNebr[iatom]; cnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + cnt];
            int jatom = info.first, ijTable = info.second;
            
            doubleVector dRij;
            double maxtime = 0;
            double sRc = iRc + particle->diameterScale[jatom] * meanRadius;
            if (istamp > particle->timeStamp[jatom]) {
                maxtime = istamp;
                double deltaTime = maxtime - particle->timeStamp[jatom];
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, -deltaTime, particle->veloc[jatom]);
            } else {
                maxtime = particle->timeStamp[jatom];
                double deltaTime = maxtime - istamp;
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, deltaTime, iveloc);
            }
            PBC(dRij, box);
            
            doubleVector dVij;
            vSub(dVij, iveloc, particle->veloc[jatom]);
            
            double A = sDot(dVij, dVij);
            double B = sDot(dRij, dVij);
            double C = sDot(dRij, dRij) - sRc * sRc;
            
            double tcolPredict = maxtime + QuadraticFormula(A, B, C);
            nebrList->colPairTable[ijTable] = tcolPredict;
            if (tcolPredict <= data.eTime) {
                data.eTime = tcolPredict;
                data.ePartner = jatom;
                data.eventType = HardCoreType;
            }
            
            if (eventList->data[jatom].ePartner != iatom &&
                tcolPredict <= eventList->data[jatom].eTime) {
                Event jdata = eventList->data[jatom];
                jdata.ePartner = iatom;
                jdata.eTime = tcolPredict;
                jdata.eventType = HardCoreType;
                eventList->data[jatom] = jdata;
                
                replaceMinEventHeap(eventList, &jdata);
            }
        }
        eventList->data[iatom] = data;
        replaceMinEventHeap(eventList, &data);
        
        // nebr
        for (int icnt = 0; icnt < nebrList->nNebr[iatom]; icnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + icnt];
            int jatom = info.first;  // ijTable = info.second;
            if (jatom == ePartner)
                continue;
            Event jdata = eventList->data[jatom];
            if (eventList->data[jatom].ePartner != iatom)
                continue;
            
            jdata.eTime = nebrList->tnebr[jatom];
            jdata.ePartner = nebrPartner;
            jdata.eventType = NebrListType;
            for (int jcnt = 0; jcnt < nebrList->nNebr[jatom]; jcnt++) {
                int2 jinfo = nebrList->list[jatom * nebrList->maxNebrPerAtom + jcnt];
                int jjatom = jinfo.first, jjTable = jinfo.second;
                if (nebrList->colPairTable[jjTable] <= jdata.eTime) {
                    jdata.eTime = nebrList->colPairTable[jjTable];
                    jdata.ePartner = jjatom;
                    jdata.eventType = HardCoreType;
                }
            }
            eventList->data[jatom] = jdata;
            replaceMinEventHeap(eventList, &jdata);
        }
    }
    
    {
        int iatom = ePartner;
        Event data;
        double iRc = particle->diameterScale[iatom] * meanRadius;
        double istamp = particle->timeStamp[iatom];
        doubleVector ipos, iveloc;
        vCpy(ipos, particle->pos[iatom]);
        vCpy(iveloc, particle->veloc[iatom]);
        // Qij
        {
            doubleVector Rt0;
            vSub(Rt0, ipos, nebrList->xyzHold[iatom]);
            doubleVector dVel;
            vCpy(dVel, iveloc);
            double rst = nebrList->rskin * 0.5;
            
            double A = -sDot(dVel, dVel);
            double B = -sDot(Rt0, dVel);
            double C = -sDot(Rt0, Rt0) + rst * rst;
            
            data.eTime = istamp + (C <= 0.0 ? 0.0 : QuadraticFormula(A, B, C));
            data.ePartner = nebrPartner;
            data.atomIdx = iatom;
            data.eventType = NebrListType;
            
            nebrList->tnebr[iatom] = data.eTime;
        }
        
        // Pij
        for (int cnt = 0; cnt < nebrList->nNebr[iatom]; cnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + cnt];
            int jatom = info.first, ijTable = info.second;
            
            doubleVector dRij;
            double maxtime = 0;
            double sRc = iRc + particle->diameterScale[jatom] * meanRadius;
            if (istamp > particle->timeStamp[jatom]) {
                maxtime = istamp;
                double deltaTime = maxtime - particle->timeStamp[jatom];
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, -deltaTime, particle->veloc[jatom]);
            } else {
                maxtime = particle->timeStamp[jatom];
                double deltaTime = maxtime - istamp;
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, deltaTime, iveloc);
            }
            PBC(dRij, box);
            
            doubleVector dVij;
            vSub(dVij, iveloc, particle->veloc[jatom]);
            
            double A = sDot(dVij, dVij);
            double B = sDot(dRij, dVij);
            double C = sDot(dRij, dRij) - sRc * sRc;
            
            double tcolPredict = maxtime + QuadraticFormula(A, B, C);
            nebrList->colPairTable[ijTable] = tcolPredict;
            if (tcolPredict <= data.eTime) {
                data.eTime = tcolPredict;
                data.ePartner = jatom;
                data.eventType = HardCoreType;
            }
            
            if (eventList->data[jatom].ePartner != iatom &&
                tcolPredict <= eventList->data[jatom].eTime) {
                Event jdata = eventList->data[jatom];
                jdata.ePartner = iatom;
                jdata.eTime = tcolPredict;
                jdata.eventType = HardCoreType;
                eventList->data[jatom] = jdata;
                
                replaceMinEventHeap(eventList, &jdata);
            }
        }
        eventList->data[iatom] = data;
        replaceMinEventHeap(eventList, &data);
        
        // nebr
        for (int icnt = 0; icnt < nebrList->nNebr[iatom]; icnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + icnt];
            int jatom = info.first;  // ijTable = info.second;
            if (jatom == eIatom)
                continue;
            Event jdata = eventList->data[jatom];
            if (eventList->data[jatom].ePartner != iatom)
                continue;
            
            jdata.eTime = nebrList->tnebr[jatom];
            jdata.ePartner = nebrPartner;
            jdata.eventType = NebrListType;
            for (int jcnt = 0; jcnt < nebrList->nNebr[jatom]; jcnt++) {
                int2 jinfo = nebrList->list[jatom * nebrList->maxNebrPerAtom + jcnt];
                int jjatom = jinfo.first, jjTable = jinfo.second;
                if (nebrList->colPairTable[jjTable] <= jdata.eTime) {
                    jdata.eTime = nebrList->colPairTable[jjTable];
                    jdata.ePartner = jjatom;
                    jdata.eventType = HardCoreType;
                }
            }
            eventList->data[jatom] = jdata;
            replaceMinEventHeap(eventList, &jdata);
        }
    }
}

//=================
void constructEventNebrList_VR(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    MinEventHeap* eventList = &update->eventList;
    double meanRadiusTime0 = 0.5 * particle->meanDiameter;
    double rrate = update->rrateSet / update->timeUnits;
    double rst = nebrList->rskin * 0.5;
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        doubleVector iveloc;
        vCpy(iveloc, particle->veloc[iatom]);
        // to avoid calculutions of t_PBC;
        double iRc = particle->diameterScale[iatom] * meanRadiusTime0;
        
        double drdt = iRc * rrate;
        double A = -sDot(iveloc, iveloc) + drdt * drdt;
        double B = -rst * drdt;
        double C = rst * rst;
        double tcolPredict = QuadraticFormula(A, B, C);
        nebrList->tnebr[iatom] = tcolPredict;
        
        eventList->data[iatom].atomIdx = iatom;
        eventList->data[iatom].eTime = tcolPredict;
        eventList->data[iatom].ePartner = nebrPartner;
        eventList->data[iatom].eventType = NebrListType;
    }
    
    bool overFlow = false;
    do {
        overFlow = false;
        int maxNebrPerAtom = nebrList->maxNebrPerAtom;
        int nColTable = 0;
        memset(nebrList->nNebr, '\0', particle->nAtom * sizeof(int));
        clearMinEventHeap(eventList, particle->nAtom);
        
#if (DIM == 3)
        for (int ibinz = 0; ibinz < nebrList->nbin[2]; ibinz++) {
            for (int ibiny = 0; ibiny < nebrList->nbin[1]; ibiny++) {
                for (int ibinx = 0; ibinx < nebrList->nbin[0]; ibinx++) {
                    int ibin =
                    (ibinz * nebrList->nbin[1] + ibiny) * nebrList->nbin[0] + ibinx;
                    for (idPosRadius idata = nebrList->binHead[ibin]; idata.id >= 0;
                         idata = nebrList->binList[idata.id]) {
                        double iRadiusRskin = idata.radius + nebrList->rskin;
                        Event data = eventList->data[idata.id];
                        
                        for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                            int jbinz =
                            (ibinz + nebrList->deltaAdjBin[adj][2]) % nebrList->nbin[2];
                            int jbiny =
                            (ibiny + nebrList->deltaAdjBin[adj][1]) % nebrList->nbin[1];
                            int jbinx =
                            (ibinx + nebrList->deltaAdjBin[adj][0]) % nebrList->nbin[0];
                            int jbin =
                            (jbinz * nebrList->nbin[1] + jbiny) * nebrList->nbin[0] +
                            jbinx;
                            for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                                 jdata = nebrList->binList[jdata.id]) {
                                if (jbin == ibin && jdata.id <= idata.id)
                                    continue;
                                
                                //============
                                doubleVector dRij;
                                vSub(dRij, idata.pos, jdata.pos);
                                PBC(dRij, box);
                                double rijP2 = sNormP2(dRij);
                                double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                                if (rijP2 > RcutRskinP2)
                                    continue;
                                
                                if (nebrList->nNebr[idata.id] < nebrList->maxNebrPerAtom &&
                                    nebrList->nNebr[jdata.id] < nebrList->maxNebrPerAtom) {
                                    double cij = idata.radius + jdata.radius;
                                    
                                    doubleVector dVij;
                                    vSub(dVij, idata.veloc, jdata.veloc);
                                    double A = sDot(dVij, dVij) - rrate * rrate * cij * cij;
                                    double B = sDot(dRij, dVij) - rrate * cij * cij;
                                    double C = sDot(dRij, dRij) - cij * cij;
                                    double tcolPredict = QuadraticFormula(A, B, C);
                                    if (nColTable >= nebrList->allocColTable) {
                                        nebrList->allocColTable += particle->nAtom * 8;
                                        nebrList->colPairTable = (double*)realloc(
                                                                                  nebrList->colPairTable,
                                                                                  nebrList->allocColTable * sizeof(double));
                                    }
                                    nebrList->colPairTable[nColTable] = tcolPredict;
                                    
                                    int2 info;
                                    info.first = jdata.id;
                                    info.second = nColTable;
                                    nebrList->list[idata.id * nebrList->maxNebrPerAtom +
                                                   nebrList->nNebr[idata.id]] = info;
                                    info.first = idata.id;
                                    info.second = nColTable;
                                    nebrList->list[jdata.id * nebrList->maxNebrPerAtom +
                                                   nebrList->nNebr[jdata.id]] = info;
                                    
                                    if (tcolPredict <= eventList->data[jdata.id].eTime) {
                                        eventList->data[jdata.id].eTime = tcolPredict;
                                        eventList->data[jdata.id].ePartner = idata.id;
                                        eventList->data[jdata.id].eventType = HardCoreType;
                                    }
                                    if (tcolPredict <= data.eTime) {
                                        data.eTime = tcolPredict;
                                        data.ePartner = jdata.id;
                                        data.eventType = HardCoreType;
                                    }
                                    
                                    nebrList->nNebr[idata.id]++;
                                    nebrList->nNebr[jdata.id]++;
                                    nColTable++;
                                } else {
                                    nebrList->nNebr[idata.id]++;
                                    nebrList->nNebr[jdata.id]++;
                                    maxNebrPerAtom =
                                    cpuMax(nebrList->nNebr[idata.id], maxNebrPerAtom);
                                    maxNebrPerAtom =
                                    cpuMax(nebrList->nNebr[jdata.id], maxNebrPerAtom);
                                    overFlow = true;
                                }
                            }  // end jbin
                            
                        }  // end adjacent bin
                        
                        eventList->data[idata.id] = data;
                    }  // end ibin
                }
            }
        }
#elif (DIM == 2)
        for (int ibiny = 0; ibiny < nebrList->nbin[1]; ibiny++) {
            for (int ibinx = 0; ibinx < nebrList->nbin[0]; ibinx++) {
                int ibin = ibiny * nebrList->nbin[0] + ibinx;
                for (idPosRadius idata = nebrList->binHead[ibin]; idata.id >= 0;
                     idata = nebrList->binList[idata.id]) {
                    double iRadiusRskin = idata.radius + nebrList->rskin;
                    Event data = eventList->data[idata.id];
                    
                    for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                        int jbiny =
                        (ibiny + nebrList->deltaAdjBin[adj][1]) % nebrList->nbin[1];
                        int jbinx =
                        (ibinx + nebrList->deltaAdjBin[adj][0]) % nebrList->nbin[0];
                        int jbin = jbiny * nebrList->nbin[0] + jbinx;
                        for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                             jdata = nebrList->binList[jdata.id]) {
                            if (jbin == ibin && jdata.id <= idata.id)
                                continue;
                            
                            //============
                            doubleVector dRij;
                            vSub(dRij, idata.pos, jdata.pos);
                            PBC(dRij, box);
                            double rijP2 = sNormP2(dRij);
                            double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                            if (rijP2 > RcutRskinP2)
                                continue;
                            
                            if (nebrList->nNebr[idata.id] < nebrList->maxNebrPerAtom &&
                                nebrList->nNebr[jdata.id] < nebrList->maxNebrPerAtom) {
                                double cij = idata.radius + jdata.radius;
                                
                                doubleVector dVij;
                                vSub(dVij, idata.veloc, jdata.veloc);
                                double A = sDot(dVij, dVij) - rrate * rrate * cij * cij;
                                double B = sDot(dRij, dVij) - rrate * cij * cij;
                                double C = sDot(dRij, dRij) - cij * cij;
                                double tcolPredict = QuadraticFormula(A, B, C);
                                if (nColTable >= nebrList->allocColTable) {
                                    nebrList->allocColTable += particle->nAtom * 8;
                                    nebrList->colPairTable = (double*)realloc(
                                                                              nebrList->colPairTable,
                                                                              nebrList->allocColTable * sizeof(double));
                                }
                                nebrList->colPairTable[nColTable] = tcolPredict;
                                
                                int2 info;
                                info.first = jdata.id;
                                info.second = nColTable;
                                nebrList->list[idata.id * nebrList->maxNebrPerAtom +
                                               nebrList->nNebr[idata.id]] = info;
                                info.first = idata.id;
                                info.second = nColTable;
                                nebrList->list[jdata.id * nebrList->maxNebrPerAtom +
                                               nebrList->nNebr[jdata.id]] = info;
                                
                                if (tcolPredict <= eventList->data[jdata.id].eTime) {
                                    eventList->data[jdata.id].eTime = tcolPredict;
                                    eventList->data[jdata.id].ePartner = idata.id;
                                    eventList->data[jdata.id].eventType = HardCoreType;
                                }
                                if (tcolPredict <= data.eTime) {
                                    data.eTime = tcolPredict;
                                    data.ePartner = jdata.id;
                                    data.eventType = HardCoreType;
                                }
                                
                                nebrList->nNebr[idata.id]++;
                                nebrList->nNebr[jdata.id]++;
                                nColTable++;
                            } else {
                                nebrList->nNebr[idata.id]++;
                                nebrList->nNebr[jdata.id]++;
                                maxNebrPerAtom =
                                cpuMax(nebrList->nNebr[idata.id], maxNebrPerAtom);
                                maxNebrPerAtom =
                                cpuMax(nebrList->nNebr[jdata.id], maxNebrPerAtom);
                                overFlow = true;
                            }
                        }  // end jbin
                        
                    }  // end adjacent bin
                    
                    eventList->data[idata.id] = data;
                }  // end ibin
            }
        }
#else
        for (int ibin = 0; ibin < nebrList->totBin; ibin++) {
            intVector vibin;
            lidx2vidx(ibin, nebrList->nbin, vibin);
            for (idPosRadius idata = nebrList->binHead[ibin]; idata.id >= 0;
                 idata = nebrList->binList[idata.id]) {
                double iRadiusRskin = idata.radius + nebrList->rskin;
                Event data = eventList->data[idata.id];
                
                for (int adj = 0; adj < nebrList->nAdjBin; adj++) {
                    intVector vjbin;
                    vAdd(vjbin, vibin, nebrList->deltaAdjBin[adj]);
                    for (int idim = 0; idim < DIM; idim++) {
                        vjbin[idim] = vjbin[idim] % nebrList->nbin[idim];
                    }
                    int jbin = vIdx2lIdx(vjbin, nebrList->nbin);
                    
                    for (idPosRadius jdata = nebrList->binHead[jbin]; jdata.id >= 0;
                         jdata = nebrList->binList[jdata.id]) {
                        if (jbin == ibin && jdata.id <= idata.id)
                            continue;
                        
                        doubleVector dRij;
                        vSub(dRij, idata.pos, jdata.pos);
                        PBC(dRij, box);
                        double rijP2 = sNormP2(dRij);
                        double RcutRskinP2 = pow(iRadiusRskin + jdata.radius, 2);
                        if (rijP2 > RcutRskinP2)
                            continue;
                        
                        if (nebrList->nNebr[idata.id] < nebrList->maxNebrPerAtom &&
                            nebrList->nNebr[jdata.id] < nebrList->maxNebrPerAtom) {
                            double cij = idata.radius + jdata.radius;
                            
                            doubleVector dVij;
                            vSub(dVij, idata.veloc, jdata.veloc);
                            double A = sDot(dVij, dVij) - rrate * rrate * cij * cij;
                            double B = sDot(dRij, dVij) - rrate * cij * cij;
                            double C = sDot(dRij, dRij) - cij * cij;
                            double tcolPredict = QuadraticFormula(A, B, C);
                            if (nColTable >= nebrList->allocColTable) {
                                nebrList->allocColTable += particle->nAtom * 8;
                                nebrList->colPairTable = (double*)realloc(nebrList->colPairTable, nebrList->allocColTable * sizeof(double));
                            }
                            nebrList->colPairTable[nColTable] = tcolPredict;
                            
                            int2 info;
                            info.first = jdata.id;
                            info.second = nColTable;
                            nebrList->list[idata.id * nebrList->maxNebrPerAtom + nebrList->nNebr[idata.id]] = info;
                            info.first = idata.id;
                            info.second = nColTable;
                            nebrList->list[jdata.id * nebrList->maxNebrPerAtom + nebrList->nNebr[jdata.id]] = info;
                            
                            if (tcolPredict <= eventList->data[jdata.id].eTime) {
                                eventList->data[jdata.id].eTime = tcolPredict;
                                eventList->data[jdata.id].ePartner = idata.id;
                                eventList->data[jdata.id].eventType = HardCoreType;
                            }
                            if (tcolPredict <= data.eTime) {
                                data.eTime = tcolPredict;
                                data.ePartner = jdata.id;
                                data.eventType = HardCoreType;
                            }
                            
                            nebrList->nNebr[idata.id]++;
                            nebrList->nNebr[jdata.id]++;
                            nColTable++;
                        } else {
                            nebrList->nNebr[idata.id]++;
                            nebrList->nNebr[jdata.id]++;
                            maxNebrPerAtom =
                            cpuMax(nebrList->nNebr[idata.id], maxNebrPerAtom);
                            maxNebrPerAtom =
                            cpuMax(nebrList->nNebr[jdata.id], maxNebrPerAtom);
                            overFlow = true;
                        }
                    }  // end jbin
                    
                }  // end adjacent bin
                
                eventList->data[idata.id] = data;
            }  // end ibin
        }
#endif
        
        if (overFlow) {
            nebrList->list = (int2*)realloc(nebrList->list, maxNebrPerAtom * particle->nAtom * sizeof(int2));
            nebrList->maxNebrPerAtom = maxNebrPerAtom;
            for (int iatom = 0; iatom < particle->nAtom; iatom++) {
                eventList->data[iatom].atomIdx = iatom;
                eventList->data[iatom].ePartner = nebrPartner;
                eventList->data[iatom].eTime = nebrList->tnebr[iatom];
                eventList->data[iatom].eventType = NebrListType;
            }
        }
        
    } while (overFlow);
}
void buildEventList_VR(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    if (nebrList->xyzHold == NULL) {
        nebrList->xyzHold =
        (doubleVector*)calloc(particle->nAtom, sizeof(doubleVector));
        nebrList->tnebr = (double*)calloc(particle->nAtom, sizeof(double));
        nebrList->nNebr = (int*)calloc(particle->nAtom, sizeof(int));
        nebrList->binList =
        (idPosRadius*)calloc(particle->nAtom, sizeof(idPosRadius));
        
        nebrList->binList4sort = (int*)calloc(particle->nAtom, sizeof(int));
        nebrList->oid2nid = (int*)calloc(particle->nAtom, sizeof(int));
        nebrList->buffer = calloc(particle->nAtom, sizeof(doubleVector));
    }
    
    if (update->stepCol >= update->nextTtstep)
        thermostat(box, particle, update);
    if (update->stepCol >= update->nextZCstep)
        moveBarycenter(box, particle, update);
    if (update->stepCol >= update->nextZMstep)
        zeroMomentum(box, particle, update);
    
    MinEventHeap* eventList = &update->eventList;
    if (eventList->isValid)
        return;
    clearMinEventHeap(eventList, particle->nAtom);
    
    syncAll(box, particle, update);
    if (getSwap(update)) {
        swapMCmove(box, particle, update);
    }
    
    if ((nebrList->nBuild % 100 == 0)) {
        sortParticle(box, particle, update);
    }
    
    // printf("%ld %g\n",nebrList->nBuild, nebrList->skinSet);
    
    adjustImg(box, particle);
    binParticle(box, particle, update);
    
    memcpy(nebrList->xyzHold, particle->pos,
           particle->nAtom * sizeof(doubleVector));
    
    constructEventNebrList_VR(box, particle, update);
    
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        Event data = eventList->data[iatom];
        replaceMinEventHeap(eventList, &data);
    }
    
    eventList->isValid = true;
    nebrList->isValid = true;
    nebrList->nBuild++;
    nebrList->colRebuild = 0;
}
void updateColPair_VR(Box* box, Particle* particle, Update* update, int eIatom, int ePartner) {
    NebrList* nebrList = &update->nebrList;
    MinEventHeap* eventList = &update->eventList;
    double meanRadiusTime0 = 0.5 * particle->meanDiameter;
    double rrate = update->rrateSet / update->timeUnits;
    
    {
        int iatom = eIatom;
        Event data;
        // to avoid calculutions of t_PBC;
        double iRc = particle->diameterScale[iatom] * meanRadiusTime0;
        double istamp = particle->timeStamp[iatom];
        doubleVector ipos, iveloc;
        vCpy(ipos, particle->pos[iatom]);
        vCpy(iveloc, particle->veloc[iatom]);
        
        // Qij
        {
            doubleVector Rt0;
            vSub(Rt0, ipos, nebrList->xyzHold[iatom]);
            doubleVector dVel;
            vCpy(dVel, iveloc);
            
            double rst = nebrList->rskin * 0.5 - iRc * rrate * istamp;
            double drdt = iRc * rrate;
            double A = sDot(dVel, dVel);
            A = -(A - drdt * drdt);
            double B = sDot(Rt0, dVel);
            B = -(B + rst * drdt);
            double C = sDot(Rt0, Rt0);
            C = -(C - rst * rst);
            
            data.eTime = istamp + (C <= 0.0 ? 0.0 : QuadraticFormula(A, B, C));
            data.ePartner = nebrPartner;
            data.atomIdx = iatom;
            data.eventType = NebrListType;
            
            nebrList->tnebr[iatom] = data.eTime;
        }
        
        // Pij
        for (int cnt = 0; cnt < nebrList->nNebr[iatom]; cnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + cnt];
            int jatom = info.first, ijTable = info.second;
            
            doubleVector dRij;
            double maxtime = 0;
            double cij = iRc + particle->diameterScale[jatom] * meanRadiusTime0;
            if (istamp > particle->timeStamp[jatom]) {
                maxtime = istamp;
                double deltaTime = maxtime - particle->timeStamp[jatom];
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, -deltaTime, particle->veloc[jatom]);
            } else {
                maxtime = particle->timeStamp[jatom];
                double deltaTime = maxtime - istamp;
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, deltaTime, iveloc);
            }
            double atime0 = 1.0 + rrate * maxtime;
            
            PBC(dRij, box);
            
            doubleVector dVij;
            vSub(dVij, iveloc, particle->veloc[jatom]);
            
            double A = sDot(dVij, dVij);
            A = A - rrate * rrate * cij * cij;
            double B = sDot(dRij, dVij);
            B = B - rrate * atime0 * cij * cij;
            double C = sDot(dRij, dRij);
            C = C - cij * cij * atime0 * atime0;
            
            double tcolPredict = maxtime + QuadraticFormula(A, B, C);
            nebrList->colPairTable[ijTable] = tcolPredict;
            if (tcolPredict <= data.eTime) {
                data.eTime = tcolPredict;
                data.ePartner = jatom;
                data.eventType = HardCoreType;
            }
            
            if (eventList->data[jatom].ePartner != iatom &&
                tcolPredict <= eventList->data[jatom].eTime) {
                Event jdata = eventList->data[jatom];
                jdata.ePartner = iatom;
                jdata.eTime = tcolPredict;
                jdata.eventType = HardCoreType;
                eventList->data[jatom] = jdata;
                
                replaceMinEventHeap(eventList, &jdata);
            }
        }
        eventList->data[iatom] = data;
        replaceMinEventHeap(eventList, &data);
        
        // nebr
        for (int icnt = 0; icnt < nebrList->nNebr[iatom]; icnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + icnt];
            int jatom = info.first;  // ijTable = info.second;
            if (jatom == ePartner)
                continue;
            Event jdata = eventList->data[jatom];
            if (eventList->data[jatom].ePartner != iatom)
                continue;
            
            jdata.eTime = nebrList->tnebr[jatom];
            jdata.ePartner = nebrPartner;
            jdata.eventType = NebrListType;
            for (int jcnt = 0; jcnt < nebrList->nNebr[jatom]; jcnt++) {
                int2 jinfo = nebrList->list[jatom * nebrList->maxNebrPerAtom + jcnt];
                int jjatom = jinfo.first, jjTable = jinfo.second;
                if (nebrList->colPairTable[jjTable] <= jdata.eTime) {
                    jdata.eTime = nebrList->colPairTable[jjTable];
                    jdata.ePartner = jjatom;
                    jdata.eventType = HardCoreType;
                }
            }
            eventList->data[jatom] = jdata;
            replaceMinEventHeap(eventList, &jdata);
        }
    }
    
    {
        int iatom = ePartner;
        Event data;
        // to avoid calculutions of t_PBC;
        double iRc = particle->diameterScale[iatom] * meanRadiusTime0;
        double istamp = particle->timeStamp[iatom];
        doubleVector ipos, iveloc;
        vCpy(ipos, particle->pos[iatom]);
        vCpy(iveloc, particle->veloc[iatom]);
        
        // Qij
        {
            doubleVector Rt0;
            vSub(Rt0, ipos, nebrList->xyzHold[iatom]);
            doubleVector dVel;
            vCpy(dVel, iveloc);
            
            double rst = nebrList->rskin * 0.5 - iRc * rrate * istamp;
            double drdt = iRc * rrate;
            double A = sDot(dVel, dVel);
            A = -(A - drdt * drdt);
            double B = sDot(Rt0, dVel);
            B = -(B + rst * drdt);
            double C = sDot(Rt0, Rt0);
            C = -(C - rst * rst);
            
            data.eTime = istamp + (C <= 0.0 ? 0.0 : QuadraticFormula(A, B, C));
            data.ePartner = nebrPartner;
            data.atomIdx = iatom;
            data.eventType = NebrListType;
            
            nebrList->tnebr[iatom] = data.eTime;
        }
        
        // Pij
        for (int cnt = 0; cnt < nebrList->nNebr[iatom]; cnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + cnt];
            int jatom = info.first, ijTable = info.second;
            
            doubleVector dRij;
            double maxtime = 0;
            double cij = iRc + particle->diameterScale[jatom] * meanRadiusTime0;
            if (istamp > particle->timeStamp[jatom]) {
                maxtime = istamp;
                double deltaTime = maxtime - particle->timeStamp[jatom];
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, -deltaTime, particle->veloc[jatom]);
            } else {
                maxtime = particle->timeStamp[jatom];
                double deltaTime = maxtime - istamp;
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, deltaTime, iveloc);
            }
            double atime0 = 1.0 + rrate * maxtime;
            
            PBC(dRij, box);
            
            doubleVector dVij;
            vSub(dVij, iveloc, particle->veloc[jatom]);
            
            double A = sDot(dVij, dVij);
            A = A - rrate * rrate * cij * cij;
            double B = sDot(dRij, dVij);
            B = B - rrate * atime0 * cij * cij;
            double C = sDot(dRij, dRij);
            C = C - cij * cij * atime0 * atime0;
            
            double tcolPredict = maxtime + QuadraticFormula(A, B, C);
            nebrList->colPairTable[ijTable] = tcolPredict;
            if (tcolPredict <= data.eTime) {
                data.eTime = tcolPredict;
                data.ePartner = jatom;
                data.eventType = HardCoreType;
            }
            
            if (eventList->data[jatom].ePartner != iatom &&
                tcolPredict <= eventList->data[jatom].eTime) {
                Event jdata = eventList->data[jatom];
                jdata.ePartner = iatom;
                jdata.eTime = tcolPredict;
                jdata.eventType = HardCoreType;
                eventList->data[jatom] = jdata;
                
                replaceMinEventHeap(eventList, &jdata);
            }
        }
        eventList->data[iatom] = data;
        replaceMinEventHeap(eventList, &data);
        
        // nebr
        for (int icnt = 0; icnt < nebrList->nNebr[iatom]; icnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + icnt];
            int jatom = info.first;  // ijTable = info.second;
            if (jatom == eIatom)
                continue;
            Event jdata = eventList->data[jatom];
            if (eventList->data[jatom].ePartner != iatom)
                continue;
            
            jdata.eTime = nebrList->tnebr[jatom];
            jdata.ePartner = nebrPartner;
            jdata.eventType = NebrListType;
            for (int jcnt = 0; jcnt < nebrList->nNebr[jatom]; jcnt++) {
                int2 jinfo = nebrList->list[jatom * nebrList->maxNebrPerAtom + jcnt];
                int jjatom = jinfo.first, jjTable = jinfo.second;
                if (nebrList->colPairTable[jjTable] <= jdata.eTime) {
                    jdata.eTime = nebrList->colPairTable[jjTable];
                    jdata.ePartner = jjatom;
                    jdata.eventType = HardCoreType;
                }
            }
            eventList->data[jatom] = jdata;
            replaceMinEventHeap(eventList, &jdata);
        }
    }
}

int checkEventList_FR(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    MinEventHeap* eventList = &update->eventList;
    int errcode = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        Event data;
        double meanRadius = 0.5 * particle->meanDiameter;
        
        double iRc = particle->diameterScale[iatom] * meanRadius;
        double istamp = particle->timeStamp[iatom];
        doubleVector ipos, iveloc;
        vCpy(ipos, particle->pos[iatom]);
        vCpy(iveloc, particle->veloc[iatom]);
        // Qij
        {
            doubleVector Rt0;
            vSub(Rt0, ipos, nebrList->xyzHold[iatom]);
            doubleVector dVel;
            vCpy(dVel, iveloc);
            double rst = nebrList->rskin * 0.5;
            
            double A = -sDot(dVel, dVel);
            double B = -sDot(Rt0, dVel);
            double C = -sDot(Rt0, Rt0) + rst * rst;
            
            data.eTime = istamp + (C <= 0.0 ? 0.0 : QuadraticFormula(A, B, C));
            data.ePartner = nebrPartner;
            data.atomIdx = iatom;
            data.eventType = NebrListType;
            
            if (data.eTime != nebrList->tnebr[iatom]) {
                errcode = -1;
                safeFprintf(stderr, "Q: %d %g %g %g\n", iatom, data.eTime,
                            nebrList->tnebr[iatom],
                            fabs(nebrList->tnebr[iatom] - data.eTime));
            }
        }
        
        // Pij
        for (int cnt = 0; cnt < nebrList->nNebr[iatom]; cnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + cnt];
            int jatom = info.first, ijTable = info.second;
            
            doubleVector dRij;
            double maxtime = 0;
            double sRc = iRc + particle->diameterScale[jatom] * meanRadius;
            if (istamp > particle->timeStamp[jatom]) {
                maxtime = istamp;
                double deltaTime = maxtime - particle->timeStamp[jatom];
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, -deltaTime, particle->veloc[jatom]);
            } else {
                maxtime = particle->timeStamp[jatom];
                double deltaTime = maxtime - istamp;
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, deltaTime, iveloc);
            }
            PBC(dRij, box);
            
            doubleVector dVij;
            vSub(dVij, iveloc, particle->veloc[jatom]);
            
            double A = sDot(dVij, dVij);
            double B = sDot(dRij, dVij);
            double C = sDot(dRij, dRij) - sRc * sRc;
            
            double tcolPredict = maxtime + QuadraticFormula(A, B, C);
            if (tcolPredict < 1e3) {
                double eTime = tcolPredict;
                
                double sRc =
                0.5 *
                (particle->diameterScale[iatom] + particle->diameterScale[jatom]) *
                particle->meanDiameter;
                
                // get deltaPij and check
                doubleVector dVij;
                vSub(dVij, particle->veloc[iatom], particle->veloc[jatom]);
                
                // update position
                doubleVector iPos, jPos;
                vScaleAdd(iPos, particle->pos[iatom],
                          eTime - particle->timeStamp[iatom], particle->veloc[iatom]);
                vScaleAdd(jPos, particle->pos[jatom],
                          eTime - particle->timeStamp[jatom], particle->veloc[jatom]);
                
                // calculate unit vector
                doubleVector vRij, Nij;
                vSub(vRij, iPos, jPos);
                PBC(vRij, box);
                double rij = sNorm(vRij);
                vScale(Nij, 1.0 / rij, vRij);
                
                // get deltaPij and check
                
                // safeFprintf(stderr, "%g <=> %g \n", rij - sRc, sRc * TolRatio);
                
                if ((sRc - rij < -DBL_EPSILON) || (sRc - rij > DBL_EPSILON)) {
                    if (fabs(rij - sRc) >= DBL_EPSILON) {
                        safeFprintf(
                                    stderr,
                                    "checkEL: %g: %d %d: %.15g - %.15g = %.15g => %10.8lf%%\n",
                                    eTime, iatom, jatom, rij, sRc, rij - sRc,
                                    (rij - sRc) / sRc * 100.0);
                    }
                    
                    if (fabs(rij - sRc) >= 100 * DBL_EPSILON) {
                        Abort("%g: %d %d: %10.8lf%%\n", eTime, iatom, jatom,
                              (rij - sRc) / sRc * 100.0);
                    }
                }
            }
            
            if (nebrList->colPairTable[ijTable] != tcolPredict) {
                errcode = -1;
                safeFprintf(stderr, "P: %d %d %g %g %g\n", iatom, jatom, tcolPredict,
                            nebrList->colPairTable[ijTable],
                            nebrList->colPairTable[ijTable] - tcolPredict);
            }
            
            if (tcolPredict <= data.eTime) {
                data.eTime = tcolPredict;
                data.ePartner = jatom;
                data.eventType = HardCoreType;
            }
        }
        
        if (data.eTime != eventList->data[iatom].eTime ||
            data.ePartner != eventList->data[iatom].ePartner) {
            errcode = -1;
            safeFprintf(stderr, "EL: %d %d %g (%d) %s %g (%d) = %g\n", iatom,
                        data.ePartner, data.eTime, data.ePartner,
                        (data.eTime < eventList->data[iatom].eTime ? "<" : ">"),
                        eventList->data[iatom].eTime,
                        eventList->data[iatom].ePartner, fabs(data.eTime - eventList->data[iatom].eTime));
        }
    }
    
    return errcode;
}
int checkEventList_VR(Box* box, Particle* particle, Update* update) {
    NebrList* nebrList = &update->nebrList;
    MinEventHeap* eventList = &update->eventList;
    double meanRadiusTime0 = 0.5 * particle->meanDiameter;
    double rrate = update->rrateSet / update->timeUnits;
    int errcode = 0;
    for (int iatom = 0; iatom < particle->nAtom; iatom++) {
        Event data;
        double meanRadius = 0.5 * particle->meanDiameter;
        
        double iRc = particle->diameterScale[iatom] * meanRadius;
        double istamp = particle->timeStamp[iatom];
        doubleVector ipos, iveloc;
        vCpy(ipos, particle->pos[iatom]);
        vCpy(iveloc, particle->veloc[iatom]);
        {  // Qij
            doubleVector Rt0;
            vSub(Rt0, ipos, nebrList->xyzHold[iatom]);
            doubleVector dVel;
            vCpy(dVel, iveloc);
            
            double rst = nebrList->rskin * 0.5 - iRc * rrate * istamp;
            double drdt = iRc * rrate;
            double A = sDot(dVel, dVel);
            A = -(A - drdt * drdt);
            double B = sDot(Rt0, dVel);
            B = -(B + rst * drdt);
            double C = sDot(Rt0, Rt0);
            C = -(C - rst * rst);
            
            data.eTime = istamp + (C <= 0.0 ? 0.0 : QuadraticFormula(A, B, C));
            data.ePartner = nebrPartner;
            data.atomIdx = iatom;
            data.eventType = NebrListType;
            
            if (data.eTime != nebrList->tnebr[iatom]) {
                errcode = -1;
                safeFprintf(stderr, "Q: %d %g %g %g\n", iatom, data.eTime,
                            nebrList->tnebr[iatom],
                            fabs(nebrList->tnebr[iatom] - data.eTime));
            }
        }
        
        // Pij
        for (int cnt = 0; cnt < nebrList->nNebr[iatom]; cnt++) {
            int2 info = nebrList->list[iatom * nebrList->maxNebrPerAtom + cnt];
            int jatom = info.first, ijTable = info.second;
            
            doubleVector dRij;
            double maxtime = 0;
            double cij = iRc + particle->diameterScale[jatom] * meanRadiusTime0;
            if (istamp > particle->timeStamp[jatom]) {
                maxtime = istamp;
                double deltaTime = maxtime - particle->timeStamp[jatom];
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, -deltaTime, particle->veloc[jatom]);
            } else {
                maxtime = particle->timeStamp[jatom];
                double deltaTime = maxtime - istamp;
                vSub(dRij, ipos, particle->pos[jatom]);
                vScaleAdd(dRij, dRij, deltaTime, iveloc);
            }
            double atime0 = 1.0 + rrate * maxtime;
            
            PBC(dRij, box);
            
            doubleVector dVij;
            vSub(dVij, iveloc, particle->veloc[jatom]);
            
            double A = sDot(dVij, dVij);
            A = A - rrate * rrate * cij * cij;
            double B = sDot(dRij, dVij);
            B = B - rrate * atime0 * cij * cij;
            double C = sDot(dRij, dRij);
            C = C - cij * cij * atime0 * atime0;
            
            double tcolPredict = maxtime + QuadraticFormula(A, B, C);
            if (nebrList->colPairTable[ijTable] != tcolPredict) {
                errcode = -1;
                safeFprintf(stderr, "P: %d %d %g %g %g\n", iatom, jatom, tcolPredict,
                            nebrList->colPairTable[ijTable],
                            nebrList->colPairTable[ijTable] - tcolPredict);
            }
            
            if (tcolPredict <= data.eTime) {
                data.eTime = tcolPredict;
                data.ePartner = jatom;
                data.eventType = HardCoreType;
            }
        }
        
        if (data.eTime != eventList->data[iatom].eTime ||
            data.ePartner != eventList->data[iatom].ePartner) {
            errcode = -1;
            safeFprintf(stderr, "EL: %d %d %g (%d) %s %g (%d) = %g\n", iatom,
                        data.ePartner, data.eTime, data.ePartner,
                        (data.eTime < eventList->data[iatom].eTime ? "<" : ">"),
                        eventList->data[iatom].eTime, eventList->data[iatom].ePartner,
                        fabs(data.eTime - eventList->data[iatom].eTime));
        }
    }
    
    return errcode;
}

//===Event Driven Molecular Dynamics===
SimEDMD* getSimEdmd(Box* box, Particle* particle, Update* update) {
    int whichTool = findToolkit(&update->toolkit, "__SimRun_EDMD__");
    if (whichTool < 0)
        return NULL;
    return (SimEDMD*)update->toolkit.toolkit[whichTool];
}
SimEDMD* addSimEdmd(Box* box, Particle* particle, Update* update) {
    if (getSimEdmd(box, particle, update)) {
        Abort("repetitive addSimEdmd()!");
    }
    syncAll(box, particle, update);
    
    SimEDMD* edmd = (SimEDMD*)calloc(1, sizeof(SimEDMD));
    addToolkit(&update->toolkit, (void*)edmd, NULL, "__SimRun_EDMD__");
    
    if (!particle->isSync)
        Abort("Fatal Error!");
    
    box->isShapeFixed = true;
    particle->isSizeFixed = true;
    
    // set paramter
    update->rrateSet = 0;
    
    update->nebrList.isValid = false;
    update->eventList.isValid = false;
    
    update->kTt = 1.0;
    update->Ttperiod = 0;  // every step
    update->nextTtstep = 0;
    
    update->colCnt = update->errCnt = 0;
    update->stepCol = update->stepErr = 0;
    
    update->nextZMstep = 0;
    update->nextZCstep = 0;
    update->nextOutputStep = 0;
    
    update->currentStamp = 0;
    update->runtimeReal = 0;
    
    update->nebrList.nBuild = 0;
    update->nebrList.accColRebuild = 0;
    
    return 0;
}
ReturnType colForward_edmd(Box* box, Particle* particle, Update* update, double runTimeReal, int maxStep) {
    // return if the following conditions are satisfied.
    // 1. number of collisions (in units of step) is equal to maxStep or stepCol >= 10000.
    // 2. the running time is equal to runTimeReal.
    SimEDMD* edmd = getSimEdmd(box, particle, update);
    if (!edmd)
        Abort("Call addSimEdmd().");
    
    update->rtype = TrivialReturn;
    
    MinEventHeap* eventList = &update->eventList;
    int colCnt = 0, stepCol = 0;
    
    if (!eventList->isValid) {
        buildEventList_FR(box, particle, update);
        update->Tdone = false;
        update->Zdone = false;
        update->Edone = false;
        particle->isSync = false;
        if (update->rtype & HaltReturn) {
            return update->rtype;
        }
    }
    
    int imaxstep = (maxStep < 0 ? INT_MAX : maxStep);
    double irunTimeReal = (runTimeReal < 0 ? DBL_INF : runTimeReal);
    double accRunTimeReal = 0.0;
    
    while (true) {
        Event data = eventList->data[eventList->eIdx[1]];
        int atomIdx = data.atomIdx;
        int partnerTag = data.ePartner;
        
        double dtime = data.eTime - update->currentStamp;
        double dtimeReal = dtime / update->timeUnits;
        if (accRunTimeReal + dtimeReal >= irunTimeReal) {
            dtime = (irunTimeReal - accRunTimeReal) * update->timeUnits;
            update->runtimeReal += (irunTimeReal - accRunTimeReal);
            update->currentStamp = data.eTime;
            update->accTimePeriod += dtime;
            
            break;
        } else {
            update->runtimeReal += dtimeReal;
            update->currentStamp = data.eTime;
            update->accTimePeriod += dtime;
            accRunTimeReal += dtimeReal;
        }
        
        if (partnerTag == nebrPartner) {
            bool breakFlag = false;
            eventList->isValid = false;
            if (update->stepCol >= update->nextOutputStep) {
                calcKinTensor(box, particle, update);
                calcPressure(box, particle, update);
                calcPoten(box, particle, update);
                update->nextOutputStep = update->stepCol + update->outputPeriod;
                if (update->Z >= JammingZmin) {
                    syncAll(box, particle, update);
                    update->rtype |= JamReturn;
                    breakFlag = true;
                }
                update->rtype |= ThermoReturn;
                breakFlag = true;
            }
            
            buildEventList_FR(box, particle, update);
            update->Tdone = false;
            update->Zdone = false;
            update->Edone = false;
            particle->isSync = false;
            if (update->rtype & HaltReturn) {
                breakFlag = true;
            }
            
            if (breakFlag)
                break;
        } else {
            doColHeapTop_FR(box, particle, update);
            updateColPair_FR(box, particle, update, atomIdx, partnerTag);
            colCnt++;
            if (colCnt == particle->nAtom) {
                colCnt = 0;
                stepCol++;
            }
        }
        
        if (stepCol >= imaxstep || stepCol >= 10000) {
            break;
        }
    }
    
    edmd->rtimeStepReal = accRunTimeReal;
    return update->rtype;
}
int delSimEdmd(Box* box, Particle* particle, Update* update) {
    syncAll(box, particle, update);
    if (!getSimEdmd(box, particle, update))
        return -1;
    calcKinTensor(box, particle, update);
    calcPoten(box, particle, update);
    calcPressure(box, particle, update);
    
    safeFprintf(stderr, "ColCnt: %lld (step: %d); InaccurateCnt: %lld (%g %%)\n",
                (long long int)(update->colCnt + update->stepCol * particle->nAtom),
                update->stepCol, (long long int)(update->errCnt + update->stepErr * particle->nAtom),
                (double)(update->errCnt + update->stepErr * particle->nAtom) / (double)(update->colCnt + update->stepCol * particle->nAtom) * 100.0);
    
    update->errCnt = update->colCnt = update->stepCol = update->stepErr = 0;
    update->runtimeReal = 0;
    update->nextTtstep = update->nextZCstep = update->nextZMstep = 0;
    delToolkit(&update->toolkit, "__SimRun_EDMD__");
    
    return 0;
}

//===Event Driven Molecular Dynamics Expand diameter===
ExpandEDMD* addExpandEdmd(Box* box, Particle* particle, Update* update, Variable* var) {
    if (getExpandEdmd(box, particle, update))
        Abort("repetitive addExpandEdmd()!");
    syncAll(box, particle, update);
    
    cmdArg* cmd = findVariable(var, "expand");
    if (!cmd || cmd->cmdArgc != 3)
        Abort("--expand rrate(\\dot{D}/D = rrate) volume/pressure 0.63/1E4");
    ExpandEDMD* expd = (ExpandEDMD*)calloc(1, sizeof(ExpandEDMD));
    addToolkit(&update->toolkit, (void*)expd, NULL, "__SimRun_ExpandEDMD__");
    
    //====expand
    box->isShapeFixed = true;
    particle->isSizeFixed = false;
    
    expd->rrateSet = (double)atof(cmd->cmdArgv[0]);
    update->rrateSet = expd->rrateSet;
    if (fabs(expd->rrateSet) <= __rrateThreshold__) {
        Abort("rrate is too small!");
    }
    if (strcmp(cmd->cmdArgv[1], "volume") == 0) {
        expd->isVolCtrl = true;
        expd->targetVF = (double)atof(cmd->cmdArgv[2]);
        if (expd->targetVF < 1E-3 || expd->targetVF > 1.0)
            Abort("targert packing fraction is invalid.");
        if ((expd->targetVF > update->volFrac && expd->rrateSet <= 0) ||
            (expd->targetVF < update->volFrac && expd->rrateSet >= 0))
            Abort("--expand rrate(\\dot{D}/D = rrate) volume/pressure 0.63/1E4");
    } else if (strcmp(cmd->cmdArgv[1], "pressure") == 0) {
        expd->isVolCtrl = false;
        expd->targetZ = (double)atof(cmd->cmdArgv[2]);
        if (expd->targetZ > JammingZmin) {
            Abort("target Z is larger than jamming criteria (%g).", JammingZmin);
        }
        if (expd->targetZ < 1.0)
            Abort("Target Z is invalid!");
        // expd->targetVF = 1.0;
    } else
        Abort("--expand rrate(\\dot{D}/D = rrate) volume/pressure 0.63/1E4");
    
    // set paramter
    update->nebrList.isValid = false;
    update->eventList.isValid = false;
    
    update->kTt = 1.0;
    update->Ttperiod = 0;  // every step
    update->nextTtstep = 0;
    
    update->colCnt = update->errCnt = 0;
    update->stepCol = update->stepErr = 0;
    
    update->nextZMstep = 0;
    update->nextZCstep = 0;
    update->nextOutputStep = 0;
    
    update->currentStamp = 0;
    update->runtimeReal = 0;
    
    update->nebrList.nBuild = 0;
    update->nebrList.accColRebuild = 0;
    
    return expd;
}
ExpandEDMD* getExpandEdmd(Box* box, Particle* particle, Update* update) {
    int whichTool = findToolkit(&update->toolkit, "__SimRun_ExpandEDMD__");
    if (whichTool < 0)
        return NULL;
    return (ExpandEDMD*)update->toolkit.toolkit[whichTool];
}
ReturnType colForward_expand(Box* box, Particle* particle, Update* update, int maxStep) {
    ExpandEDMD* expd = getExpandEdmd(box, particle, update);
    if (expd == NULL)
        Abort("Call addExpandEdmd().");
    
    update->rtype = TrivialReturn;
    
    bool isVolCtrl = expd->isVolCtrl;
    double tVF = expd->targetVF, tZ = expd->targetZ;
    double tmax = DBL_INF;
    if (isVolCtrl) {
        tmax = (pow(tVF / update->volFrac, 1.0 / DIM) - 1.0) / (update->rrateSet / update->timeUnits);
    }
    int imaxstep = (maxStep < 0 ? INT_MAX : maxStep);
    
    MinEventHeap* eventList = &update->eventList;
    int colCnt = 0, stepCol = 0;
    if (!eventList->isValid) {
        buildEventList_VR(box, particle, update);
        
        if (isVolCtrl) {
            tmax = (pow(tVF / update->volFrac, 1.0 / DIM) - 1.0) / (update->rrateSet / update->timeUnits);
        }
        
        update->Tdone = false;
        update->Zdone = false;
        update->Edone = false;
        particle->isSync = false;
        if (update->rtype & HaltReturn) {
            return update->rtype;
        }
    }
    
    double accRunTimeReal = 0.0;
    while (true) {
        Event data = eventList->data[eventList->eIdx[1]];
        if (data.ePartner != nebrPartner && update->nebrList.colRebuild >= 1.2 * particle->nAtom) {
            // modify Event Structure to rebuild NebrList;
            data.ePartner = nebrPartner;
            data.eTime = update->currentStamp + (data.eTime - update->currentStamp) * 0.5;
            
            //====estimate rskin====
            double skinSet = update->nebrList.skinSet * 0.95;
            // estimated from spaces between particles
            NebrList* nebrList = &update->nebrList;
            double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
            double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
            double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
            // estimated from box side length
            doubleVector distPlane;
            calcDistBoxPlane(distPlane, box->boxEdge);
            double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
            for (int idim = 0; idim < DIM; idim++) {
                minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
            }
            double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
            double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
            double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
            double lmaxset = cpuMin(lmsset, maxRskinSet);
            skinSet = cpuMin(skinSet, lmaxset);
            nebrList->skinSet = skinSet;
            nebrList->compelInit = true;
        }
        
        double dtime = data.eTime - update->currentStamp;
        double dtimeReal = 0.0;
        if (data.eTime >= tmax) {  // reaching target density
            dtime = tmax - update->currentStamp;
            dtimeReal = log(1.0 + update->rrateSet / update->timeUnits * dtime) / (update->rrateSet);
            
            // update->duringTime += dtime;
            update->runtimeReal += dtimeReal;
            update->currentStamp = data.eTime;
            update->accTimePeriod += dtime;
            accRunTimeReal += dtimeReal;
            
            syncAll(box, particle, update);
            eventList->isValid = false;
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            update->rtype |= DensityReturn;
            break;
        }
        if (1.0 + update->rrateSet / update->timeUnits * dtime < 1E-10) {  // reaching zero density
            dtime = (1E-10 - 1.0) / (update->rrateSet / update->timeUnits);
            dtimeReal = log(1E-10) / (update->rrateSet);
            
            // update->duringTime += dtime;
            update->runtimeReal += dtimeReal;
            update->currentStamp = data.eTime;
            update->accTimePeriod += dtime;
            accRunTimeReal += dtimeReal;
            
            syncAll(box, particle, update);
            eventList->isValid = false;
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            update->rtype |= DensityReturn;
            break;
        }
        dtimeReal = log(1.0 + update->rrateSet / update->timeUnits * dtime) / (update->rrateSet);
        
        // update->duringTime += dtime;
        update->accTimePeriod += dtime;
        update->runtimeReal += dtimeReal;
        update->currentStamp = data.eTime;
        accRunTimeReal += dtimeReal;
        
        int atomIdx = data.atomIdx;
        int partnerTag = data.ePartner;
        if (partnerTag == nebrPartner) {
            if (update->nebrList.colRebuild <= 0.8 * particle->nAtom) {
                double skinSet = update->nebrList.skinSet * 1.05;
                // estimated from spaces between particles
                NebrList* nebrList = &update->nebrList;
                double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
                double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
                double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
                // estimated from box side length
                doubleVector distPlane;
                calcDistBoxPlane(distPlane, box->boxEdge);
                double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
                for (int idim = 0; idim < DIM; idim++) {
                    minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
                }
                double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
                double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
                double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
                double lmaxset = cpuMin(lmsset, maxRskinSet);
                skinSet = cpuMin(skinSet, lmaxset);
                nebrList->skinSet = skinSet;
                nebrList->compelInit = true;
            }
            
            eventList->isValid = false;
            calcKinTensor(box, particle, update);
            calcPressure(box, particle, update);
            calcPoten(box, particle, update);
            
            bool breakFlag = false;
            if (isVolCtrl && update->Z >= JammingZmin) {
                syncAll(box, particle, update);
                update->rtype |= JamReturn;
                breakFlag = true;
            }
            if (!isVolCtrl && update->Z >= tZ) {
                syncAll(box, particle, update);
                update->rtype |= PressReturn;
                breakFlag = true;
            }
            if (update->stepCol >= update->nextOutputStep) {
                update->nextOutputStep = update->stepCol + update->outputPeriod;
                update->rtype |= ThermoReturn;
                breakFlag = true;
            }
            buildEventList_VR(box, particle, update);
            if (isVolCtrl) {
                tmax = (pow(tVF / update->volFrac, 1.0 / DIM) - 1.0) / (update->rrateSet / update->timeUnits);
            }
            update->Tdone = false;
            update->Zdone = false;
            update->Edone = false;
            particle->isSync = false;
            if (update->rtype & HaltReturn) {
                breakFlag = true;
            }
            
            if (breakFlag)
                break;
        } else {
            doColHeapTop_VR(box, particle, update);
            updateColPair_VR(box, particle, update, atomIdx, partnerTag);
            
            colCnt++;
            if (colCnt == particle->nAtom) {
                colCnt = 0;
                stepCol++;
            }
        }
        
        if (stepCol >= imaxstep || stepCol >= 1000) {
            break;
        }
    }
    
    expd->rtimeStepReal = accRunTimeReal;
    return update->rtype;
}
int delExpandEdmd(Box* box, Particle* particle, Update* update) {
    syncAll(box, particle, update);
    if (!getExpandEdmd(box, particle, update))
        return -1;
    
    calcKinTensor(box, particle, update);
    calcPoten(box, particle, update);
    calcPressure(box, particle, update);
    
    safeFprintf(stderr, "ColCnt: %lld (step: %d); InaccurateCnt: %lld (%g %%)\n",
                (long long int)(update->colCnt + update->stepCol * particle->nAtom),
                update->stepCol, (long long int)(update->errCnt + update->stepErr * particle->nAtom),
                (double)(update->errCnt + update->stepErr * particle->nAtom) / (double)(update->colCnt + update->stepCol * particle->nAtom) * 100.0);
    
    update->errCnt = update->colCnt = update->stepCol = update->stepErr = 0;
    update->runtimeReal = 0;
    update->nextTtstep = update->nextZCstep = update->nextZMstep = 0;
    update->rrateSet = 0;
    delToolkit(&update->toolkit, "__SimRun_ExpandEDMD__");
    
    return 0;
}

//===Const-Z (both P and T) Event Driven Molecular Dynamics===
double estimateRate(Box* box, Particle* particle, Update* update, NptIsoEDMD* nptiso) {
    // update units of momentum;
    double boxWg = particle->nAtom * pow(nptiso->tauZ * update->timeUnits, 2);
    double boxPg = 0;
    int nStep = 10.0;
    double dt = nptiso->predtStep / nStep;
    double alpha = dt / (nptiso->tauZ * update->timeUnits);
    double sigma = sqrt(2.0 * alpha / boxWg);
    double b = 1.0 / (1.0 + 0.5 * alpha);                  // unit is 1
    double a = (1.0 - 0.5 * alpha) / (1.0 + 0.5 * alpha);  // unit is 1
    double gForce = particle->nAtom * (update->Z * update->Tint - nptiso->Ztarg * update->kTt);
    for (int ith = 0; ith < nStep; ith++) {
        double beta = rndStdNorm() * sigma;
        boxPg = boxPg * a + a * dt * gForce / boxWg + b * beta;
    }
    
    // expand particles
    double rrate = (exp(-boxPg * nptiso->predtStep) - 1.0) / nptiso->predtStep;
    if (isinf(rrate)) {
        //Info("Warning! rrate: %g.", rrate,boxPg,nptiso->predtStep);
        rrate = (boxPg > 0 ? -__maxRrate__ : __maxRrate__) / update->timeUnits;
    }
    update->rrateSet = rrate * update->timeUnits;
    if (fabs(rrate) < __rrateThreshold__ / update->timeUnits) {
        rrate = 0.0;
        update->rrateSet = 0.0;
        nptiso->zeroRate = true;
    } else if (fabs(rrate) >= __maxRrate__ / update->timeUnits) {
        rrate = rrate / fabs(rrate) * __maxRrate__ / update->timeUnits;
        update->rrateSet = rrate * update->timeUnits;
    }
    
    update->eventList.isValid = false;
    return rrate;
}
ReturnType npt_firstStep(Box* box, Particle* particle, Update* update, NptIsoEDMD* nptiso) {
    MinEventHeap* eventList = &update->eventList;
    
    update->rrateSet = 0.0;
    update->rtype = TrivialReturn;
    int colCnt = 0, stepCol = 0;
    if (!eventList->isValid) {
        buildEventList_VR(box, particle, update);
        particle->isSync = false;
        if (update->rtype & HaltReturn) {
            update->rtype &= HaltReturn;
        }
    }
    
    while (true) {
        Event data = eventList->data[eventList->eIdx[1]];
        if (data.ePartner != nebrPartner && update->nebrList.colRebuild >= 1.5 * particle->nAtom) {
            // adjust NebrList parameter
            // modify Event Structure to rebuild NebrList;
            data.ePartner = nebrPartner;
            data.eTime = update->currentStamp + (data.eTime - update->currentStamp) * 0.5;
            
            //====estimate rskin====
            double skinSet = update->nebrList.skinSet * 0.95;
            // estimated from spaces between particles
            NebrList* nebrList = &update->nebrList;
            double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
            double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
            double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
            // estimated from box side length
            doubleVector distPlane;
            calcDistBoxPlane(distPlane, box->boxEdge);
            double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
            for (int idim = 0; idim < DIM; idim++) {
                minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
            }
            double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
            double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
            double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
            double lmaxset = cpuMin(lmsset, maxRskinSet);
            skinSet = cpuMin(skinSet, lmaxset);
            nebrList->skinSet = skinSet;
            nebrList->compelInit = true;
        }
        int atomIdx = data.atomIdx;
        int partnerTag = data.ePartner;
        
        double dtime = data.eTime - update->currentStamp;
        double dtimeReal = dtime / update->timeUnits;
        update->runtimeReal += dtimeReal;
        update->currentStamp = data.eTime;
        update->accTimePeriod += dtime;
        
        if (partnerTag == nebrPartner) {
            if (update->nebrList.colRebuild <= 1.1 * particle->nAtom) {
                double skinSet = update->nebrList.skinSet * 1.05;
                // estimated from spaces between particles
                NebrList* nebrList = &update->nebrList;
                double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
                double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
                double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
                // estimated from box side length
                doubleVector distPlane;
                calcDistBoxPlane(distPlane, box->boxEdge);
                double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
                for (int idim = 0; idim < DIM; idim++) {
                    minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
                }
                double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
                double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
                double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
                double lmaxset = cpuMin(lmsset, maxRskinSet);
                skinSet = cpuMin(skinSet, lmaxset);
                nebrList->skinSet = skinSet;
                nebrList->compelInit = true;
            }
            
            eventList->isValid = false;
            buildEventList_VR(box, particle, update);
            particle->isSync = false;
            if (update->rtype & HaltReturn) {
                update->rtype &= HaltReturn;
            }
            
        } else {
            doColHeapTop_VR(box, particle, update);
            updateColPair_VR(box, particle, update, atomIdx, partnerTag);
            colCnt++;
            if (colCnt == particle->nAtom) {
                colCnt = 0;
                stepCol++;
                break;
            }
        }
    }
    nptiso->predtStep = update->accTimePeriod;
    
    calcKinTensor(box, particle, update);
    calcPressure(box, particle, update);
    calcPoten(box, particle, update);
    
    return update->rtype;
}
ReturnType npt_step(Box* box, Particle* particle, Update* update, NptIsoEDMD* nptiso) {
    MinEventHeap* eventList = &update->eventList;
    int colCnt = 0, stepCol = 0;
    update->rtype = TrivialReturn;
    
    // update units of momentum;
    double rrate = estimateRate(box, particle, update, nptiso);
    buildEventList_VR(box, particle, update);
    particle->isSync = false;
    
    double accDuringTime = 0;
    while (true) {
        Event data = eventList->data[eventList->eIdx[1]];
        if (data.ePartner != nebrPartner && update->nebrList.colRebuild >= 1.5 * particle->nAtom) {
            // modify Event Structure to rebuild NebrList;
            data.ePartner = nebrPartner;
            data.eTime = update->currentStamp + (data.eTime - update->currentStamp) * 0.5;
            
            //====estimate rskin====
            double skinSet = update->nebrList.skinSet * 0.95;
            // estimated from spaces between particles
            NebrList* nebrList = &update->nebrList;
            double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
            double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
            double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
            // estimated from box side length
            doubleVector distPlane;
            calcDistBoxPlane(distPlane, box->boxEdge);
            double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
            for (int idim = 0; idim < DIM; idim++) {
                minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
            }
            double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
            double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
            double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
            double lmaxset = cpuMin(lmsset, maxRskinSet);
            skinSet = cpuMin(skinSet, lmaxset);
            nebrList->skinSet = skinSet;
            nebrList->compelInit = true;
        }
        
        double dtime = data.eTime - update->currentStamp;
        double dtimeReal = 0.0;
        
        if (1.0 + rrate * dtime < 1E-10) {  // reaching zero density
            dtime = (1E-10 - 1.0) / (rrate);
            dtimeReal = log(1E-10) / (rrate * update->timeUnits);
            
            update->runtimeReal += dtimeReal;
            update->currentStamp = data.eTime;
            update->accTimePeriod += dtime;
            
            eventList->isValid = false;
            update->rtype |= ErrorReturn;
            break;
        }
        if (nptiso->zeroRate)
            dtimeReal = dtime / update->timeUnits;
        else
            dtimeReal = log(1.0 + rrate * dtime) / (rrate * update->timeUnits);
        
        update->runtimeReal += dtimeReal;
        update->currentStamp = data.eTime;
        update->accTimePeriod += dtime;
        accDuringTime += dtime;
        
        int atomIdx = data.atomIdx;
        int partnerTag = data.ePartner;
        
        if (partnerTag == nebrPartner) {
            if (update->nebrList.colRebuild <= 1.1 * particle->nAtom) {
                double skinSet = update->nebrList.skinSet * 1.05;
                // estimated from spaces between particles
                NebrList* nebrList = &update->nebrList;
                double lbox = pow(box->volume / particle->nAtom / VolUnitSphere, 1.0 / DIM);
                double lsph = pow(box->volume * update->volFrac / particle->nAtom / VolUnitSphere, 1.0 / DIM);
                double lmsset = (lbox - lsph) / (nebrList->minDiameterScale * particle->meanDiameter);
                // estimated from box side length
                doubleVector distPlane;
                calcDistBoxPlane(distPlane, box->boxEdge);
                double minAxByCz = box->boxH[spaceIdx2voigt(0, 0)];
                for (int idim = 0; idim < DIM; idim++) {
                    minAxByCz = cpuMin(minAxByCz, box->boxH[spaceIdx2voigt(idim, idim)]);
                }
                double maxRcut = nebrList->maxDiameterScale * particle->meanDiameter;
                double maxRskin = (0.5 * minAxByCz - maxRcut) * 0.5 * 0.99;
                double maxRskinSet = maxRskin / (nebrList->minDiameterScale * particle->meanDiameter);
                double lmaxset = cpuMin(lmsset, maxRskinSet);
                skinSet = cpuMin(skinSet, lmaxset);
                nebrList->skinSet = skinSet;
                nebrList->compelInit = true;
            }
            
            eventList->isValid = false;
            buildEventList_VR(box, particle, update);
            if (update->rtype & HaltReturn) {
                update->rtype &= HaltReturn;
            }
        } else {
            update->Tdone = false;
            update->Zdone = false;
            update->Edone = false;
            particle->isSync = false;
            
            doColHeapTop_VR(box, particle, update);
            updateColPair_VR(box, particle, update, atomIdx, partnerTag);
            
            colCnt++;
            if (colCnt == particle->nAtom) {
                colCnt = 0;
                stepCol++;
                break;
            }
        }
    }
    nptiso->predtStep = accDuringTime;
    
    calcKinTensor(box, particle, update);
    calcPressure(box, particle, update);
    calcPoten(box, particle, update);
    
    return update->rtype;
}

NptIsoEDMD* getNptEdmd(Box* box, Particle* particle, Update* update) {
    int whichTool = findToolkit(&update->toolkit, "__SimRun_NptIsoEDMD__");
    if (whichTool < 0)
        return NULL;
    return (NptIsoEDMD*)update->toolkit.toolkit[whichTool];
}
NptIsoEDMD* addNptEdmd(Box* box, Particle* particle, Update* update, Variable* var) {
    cmdArg* cmd = findVariable(var, "npt");
    if (cmd == NULL || cmd->cmdArgc != 2)
        Abort("--npt Ztarget tauZ");
    // Abort("--npt runtime Zstart Zstop tauZ");
    
    NptIsoEDMD* nptiso = (NptIsoEDMD*)calloc(1, sizeof(NptIsoEDMD));
    addToolkit(&update->toolkit, (void*)nptiso, NULL, "__SimRun_NptIsoEDMD__");
    
    box->isShapeFixed = true;
    particle->isSizeFixed = false;
    
    // set paramter
    nptiso->Zstart = atof(cmd->cmdArgv[0]);
    nptiso->Zstop = nptiso->Zstart;
    nptiso->tauZ = atof(cmd->cmdArgv[1]);
    nptiso->Zrate = 0.0;
    if (nptiso->Zstart <= 1.0 || nptiso->Zstop <= 1.0 || nptiso->tauZ <= 0.0 || fabs(nptiso->Zrate) >= 1.0
        /*|| nptiso->runTimeRealSet <= 0 */) {
        Abort("--npt Ztarget tauZ");
        // Abort("--npt runtime Zstart Zstop tauZ");
    }
    
    // init nptiso
    nptiso->Ztarg = nptiso->Zstart;
    nptiso->predtStep = 0;
    
    update->nebrList.isValid = false;
    update->nebrList.compelInit = true;
    update->nebrList.nBuild = 0;
    update->nebrList.maxBinHead = 0;
    update->nebrList.maxNebrPerAtom = 0;
    update->eventList.isValid = false;
    
    update->kTt = 1.0;
    update->runtimeReal = 0;
    
    update->Ttperiod = 1;
    update->nextTtstep = 0;
    
    update->colCnt = update->errCnt = 0;
    update->stepCol = update->stepErr = 0;
    update->nextZMstep = 0;
    update->nextZCstep = 0;
    
    update->Edone = update->Zdone = update->Tdone = false;
    update->nextOutputStep = 0;
    
    // estimate instant Zint.
    npt_firstStep(box, particle, update, nptiso);
    
    return nptiso;
}
ReturnType stepForward_npt(Box* box, Particle* particle, Update* update, int maxStep) {
    NptIsoEDMD* nptiso = getNptEdmd(box, particle, update);
    if (!nptiso)
        Abort("Call addNptEdmd().");
    
    int imaxstep = (maxStep < 0 ? INT_MAX : maxStep);
    int step0 = update->stepCol;
    
    do {
        npt_step(box, particle, update, nptiso);
        
        if (update->rtype & ErrorReturn) {
            break;
        }
        
        bool breakFlag = false;
        if (update->stepCol - step0 >= imaxstep) {
            breakFlag = true;
        }
        if (update->Z >= JammingZmin) {
            update->rtype |= JamReturn;
            breakFlag = true;
        }
        if (update->stepCol >= update->nextOutputStep) {
            update->nextOutputStep = update->stepCol + update->outputPeriod;
            update->rtype |= ThermoReturn;
            breakFlag = true;
        }
        if (breakFlag)
            break;
    } while (update->stepCol - step0 < 10000);
    
    return update->rtype;
}
int delNptEdmd(Box* box, Particle* particle, Update* update) {
    syncAll(box, particle, update);
    if (!getNptEdmd(box, particle, update))
        return -1;
    
    calcKinTensor(box, particle, update);
    calcPoten(box, particle, update);
    calcPressure(box, particle, update);
    
    safeFprintf(stderr, "ColCnt: %lld (step: %d); InaccurateCnt: %lld (%g %%)\n",
                (long long int)(update->colCnt + update->stepCol * particle->nAtom),
                update->stepCol, (long long int)(update->errCnt + update->stepErr * particle->nAtom),
                (double)(update->errCnt + update->stepErr * particle->nAtom) / (double)(update->colCnt + update->stepCol * particle->nAtom) * 100.0);
    
    update->errCnt = update->colCnt = update->stepCol = update->stepErr = 0;
    update->runtimeReal = 0;
    update->nextTtstep = update->nextZCstep = update->nextZMstep = 0;
    update->rrateSet = 0;
    delToolkit(&update->toolkit, "__SimRun_NptIsoEDMD__");
    
    return 0;
}
