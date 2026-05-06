#ifndef MULTI_IMPULSE_TRANSFER_OPTIMIZATION_H
#define MULTI_IMPULSE_TRANSFER_OPTIMIZATION_H

/**
 * MultiImpulseTransferOptimization.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-8-9
 */

#include "MultiImpulseTransfer.h"

class MultiImpulseTransferOptimization {
public:
    MultiImpulseTransferOptimization() = default;

    MultiImpulseTransferOptimization(int maxImpulseNum_, double startTime_, double endTime_, int targetId0,
                                     int targetId1, double dvMaxSingle_,
                                     int impulseNum_, double rMin_, Dynamics dynamicsModel_);

    MultiImpulseTransferOptimization(int maxImpulseNum_, const MultiImpulseTransfer &mit);

    ~MultiImpulseTransferOptimization() = default;

    int maxImpulseNum{5};
    double t0;
    double tf;
    std::array<double, 6> rv0;
    std::array<double, 6> rvf;
    double dRAAN;
    double dH;
    int ratioFlag, P1, S1;

    MultiImpulseTransfer transfer;

    void optimization(MultiImpulseTransfer::MITResult &result);

    void optimization_short(MultiImpulseTransfer::MITResult &result);

};


#endif //MULTI_IMPULSE_TRANSFER_OPTIMIZATION_H
