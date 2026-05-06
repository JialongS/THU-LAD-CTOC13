#include "MultiImpulseTransferOptimization.h"
#include "optimization.h"

#include <algorithm>
#include <constant.h>

MultiImpulseTransferOptimization::MultiImpulseTransferOptimization(int maxImpulseNum_, double startTime_,
                                                                   double endTime_, int targetId0, int targetId1,
                                                                   double dvMaxSingle_, int impulseNum_, double rMin_,
                                                                   Dynamics dynamicsModel_) : maxImpulseNum{
        maxImpulseNum_
}, transfer(startTime_, endTime_, targetId0, targetId1, dvMaxSingle_,
            impulseNum_, rMin_, dynamicsModel_) {
}

MultiImpulseTransferOptimization::MultiImpulseTransferOptimization(int maxImpulseNum_, const MultiImpulseTransfer &mit)
        :
        maxImpulseNum{maxImpulseNum_}, transfer(mit) {
}

void MultiImpulseTransferOptimization::optimization(MultiImpulseTransfer::MITResult &result) {
    transfer.rv0_array = rv0;
    transfer.rvf_array = rvf;
    transfer.startTime = t0;
    transfer.endTime = tf;

    // std::vector<double> fbestList(maxImpulseNum - 1);
    // std::vector<MultiImpulseTransfer::MITResult> resultList;
    // for (int i = 5; i <= maxImpulseNum; i++)
    // {
    transfer.impulseNum = 3;
    transfer.dvMaxSingle = 300.0;

    const int numVar = transfer.calcNumVars();
    std::vector<double> xbest(numVar);
    double fbest;

    transfer.dynamicModel = Dynamics::ACCURACY;
    transfer.rMin = RADIUS_EARTH + 200e3;
    PSO_parallel2(MITOAdapter, &transfer, xbest, fbest, numVar, 10 * numVar, 100, 20);

    //transfer.dynamicModel = Dynamics::ACCURACY;
    // MultiImpulseTransfer::MITResult result;
    //transfer.rMin = 0;
    transfer.getResult(xbest, result);
    //     fbestList[i - 2] = resultI.dvTotal;
    //     resultList.push_back(resultI);
    // }
    //
    // const auto minIter = std::min_element(fbestList.begin(), fbestList.end());
    // const size_t bestIndex = std::distance(fbestList.begin(), minIter);
    // result = resultList[bestIndex];
    // result.mergeNode(1e-3);
}

void MultiImpulseTransferOptimization::optimization_short(MultiImpulseTransfer::MITResult &result) {
    transfer.rv0_array = rv0;
    transfer.rvf_array = rvf;
    transfer.startTime = t0;
    transfer.endTime = tf;
    transfer.dRAAN = dRAAN;
    transfer.ratioFlag = ratioFlag;
    transfer.P1 = P1;
    transfer.S1 = S1;
    transfer.dH = dH;

    // std::vector<double> fbestList(maxImpulseNum - 1);
    // std::vector<MultiImpulseTransfer::MITResult> resultList;
    // for (int i = 5; i <= maxImpulseNum; i++)
    // {
    transfer.impulseNum = 3;
    transfer.dvMaxSingle = 300.0;

    const int numVar = 4 * transfer.impulseNum - 5;
    std::vector<double> xbest(numVar);
    double fbest;

    transfer.dynamicModel = Dynamics::ACCURACY;
    transfer.rMin = RADIUS_EARTH + 200e3;

    PSO_parallel2(MITOAdapter_short, &transfer, xbest, fbest, numVar, 10 * numVar, 100, 20);

    //transfer.dynamicModel = Dynamics::ACCURACY;
    // MultiImpulseTransfer::MITResult result;
    //transfer.rMin = 0;
    transfer.getResult_short(xbest, result);
    //     fbestList[i - 2] = resultI.dvTotal;
    //     resultList.push_back(resultI);
    // }
    //
    // const auto minIter = std::min_element(fbestList.begin(), fbestList.end());
    // const size_t bestIndex = std::distance(fbestList.begin(), minIter);
    // result = resultList[bestIndex];
    // result.mergeNode(1e-3);
}
