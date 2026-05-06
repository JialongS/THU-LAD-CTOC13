#include "MultiImpulseTransfer.h"
#include "dynamics.h"
#include "J2Lambert.h"
#include "OrbitFun.h"
#include "data_gtoc9.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <dace/dacebase.h>
#include "../InPlaneOptimization_tym.h"


Target::Target(int id_) : id{id_} {
}

Target::Target(const Target &tar) {
    this->id = tar.id;
}

void Target::getRv(double t, double *rv) const {
    getDebrisRV(id, t, rv);
}

MultiImpulseTransfer::MultiImpulseTransfer(const double startTime_, const double endTime_, const int targetId0,
                                           const int targetId1,
                                           const double dvMaxSingle_, const int impulseNum_, const double rMin_,
                                           const Dynamics dynamicsModel_) {
    startTime = startTime_;
    endTime = endTime_;
    target0 = Target(targetId0);
    target1 = Target(targetId1);
    dvMaxSingle = dvMaxSingle_;
    impulseNum = impulseNum_;
    rMin = rMin_;
    dynamicModel = dynamicsModel_;
}

MultiImpulseTransfer::MultiImpulseTransfer(const MultiImpulseTransfer &mit) {
    startTime = mit.startTime;
    endTime = mit.endTime;
    target0 = mit.target0;
    target1 = mit.target1;
    dvMaxSingle = mit.dvMaxSingle;
    impulseNum = mit.impulseNum;
    rMin = mit.rMin;
    dynamicModel = mit.dynamicModel;
    rv0_array = mit.rv0_array;
    rvf_array = mit.rvf_array;
    P1 = mit.P1;
    S1 = mit.S1;
    ratioFlag = mit.ratioFlag;
    dRAAN = mit.dRAAN;
    dH = mit.dH;
    //tImpulse = mit.tImpulse;
    //dvImpulse = mit.dvImpulse;
    //rvImpulse = mit.rvImpulse;
}

// 计算优化变量维数 N = 4 * impulseNum - 6
int MultiImpulseTransfer::calcNumVars() const {
    return 4 * impulseNum - 6;
}

void MultiImpulseTransfer::propagateTrajectory(double t0, double tf, const double *rv0, double *rvf) const {
    if (dynamicModel == Dynamics::ACCURACY) {
        int flag = orbitPropagateJ2(t0, tf, rv0, rvf);
    } else if (dynamicModel == Dynamics::APPROXIMATE) {
        orbitPropagateJ2Approximation(t0, tf, rv0, rvf);
    }
}

void MultiImpulseTransfer::solveLambert(double tof, const double *rv1, const double *rv2, double *v1, double *v2,
                                        int &flag, double mu) const {
    if (dynamicModel == Dynamics::APPROXIMATE)
        J2_Lambert_MRPLP(v1, v2, rv1, rv2, tof, flag, mu, 200, 100);
    if (dynamicModel == Dynamics::ACCURACY)
        J2_Lambert_MRPLP(v1, v2, rv1, rv2, tof, flag, mu, 200, 1e-6);
}

bool MultiImpulseTransfer::isBiggerThanRmin(const double *rv, double mu) {
    double coe[6];
    int flag;
    rv2coe(flag, coe, rv, mu);
    //std::cout << flag << std::endl;
    //std::cout << coe[0] * (1 - coe[1]) << std::endl;
    if (!flag)
        return false;
    return coe[0] * (1 - coe[1]) > rMin;
}

double MultiImpulseTransfer::transform(const std::vector<double> &x) {
    tImpulse.clear();
    tImpulse.resize(impulseNum + 1);
    dvImpulse.clear();
    dvImpulse.resize(impulseNum + 1);
    rvImpulse.clear();
    rvImpulse.resize(impulseNum + 1);

    double endTime_new = endTime + 3600.0 * (x[4 * impulseNum - 7] - 0.5);


    double t0 = startTime + (endTime - startTime) * x[0];
    double rv[6];
    for (int i = 0; i < 6; i++) {
        rv[i] = rv0_array[i];
    }
    double rvf[6];
    //target0.getRv(t0, rv);
    if (!isBiggerThanRmin(rv, MU_EARTH)) {
        //std::cout << "rv0!\n";
        return 1e4;
    }

    propagateTrajectory(startTime, t0, rv, rvf);
    std::copy_n(rvf, 6, rv);
    if (!isBiggerThanRmin(rv, MU_EARTH)) {
        //std::cout << "rv0!\n";
        return 1e4;
    }

    for (int i = 0; i < 3; i++)
        dvImpulse[0][i] = 0.0;


    tImpulse[0] = startTime;
    tImpulse[1] = t0;

    for (int i = 0; i < 6; i++)
        rvImpulse[0][i] = rv0_array[i];

    // 深空机动
    for (int i = 1; i < impulseNum - 1; i++) {
        tImpulse[i + 1] = tImpulse[i] + (endTime - tImpulse[i]) * x[4 * i - 3];
        double dvNorm = dvMaxSingle * x[4 * i - 2];
        double phi = DPI * (x[4 * i - 1] - 0.5);
        double theta = D2PI * x[4 * i];

        dvImpulse[i][0] = std::cos(phi) * std::cos(theta) * dvNorm;
        dvImpulse[i][1] = std::cos(phi) * std::sin(theta) * dvNorm;
        dvImpulse[i][2] = std::sin(phi) * dvNorm;
        for (int j = 0; j < 3; j++) {
            rv[j + 3] += dvImpulse[i][j];
        }
        if (!isBiggerThanRmin(rv, MU_EARTH)) {
            //std::cout << "rv!\n";
            return 1e4;
        }
        std::copy_n(rv, 6, rvImpulse[i].begin());

        propagateTrajectory(tImpulse[i], tImpulse[i + 1], rv, rvf);
        std::copy_n(rvf, 6, rv);
    }

    // 解Lambert问题
    tImpulse[impulseNum] = endTime;
    for (int i = 0; i < 6; i++) {
        rvf[i] = rvf_array[i];
    }
    //target1.getRv(tImpulse[impulseNum], rvf);


    double tof = tImpulse[impulseNum] - tImpulse[impulseNum - 1];
    double v1[3], v2[3];
    int lambert_flag;
    solveLambert(tof, rv, rvf, v1, v2, lambert_flag, MU_EARTH);
    if (lambert_flag > 0)
        return 1e5;
    for (int i = 0; i < 3; i++) {
        dvImpulse[impulseNum - 1][i] = v1[i] - rv[i + 3];
        dvImpulse[impulseNum][i] = rvf[i + 3] - v2[i];
    }
    for (int i = 0; i < 3; i++) {
        rvImpulse[impulseNum - 1][i] = rv[i];
        rvImpulse[impulseNum - 1][i + 3] = v1[i];
        rvImpulse[impulseNum][i] = rvf[i];
        rvImpulse[impulseNum][i + 3] = rvf[i + 3];
    }
    if (!isBiggerThanRmin(rvImpulse[impulseNum - 1].data(), MU_EARTH)) {
        //std::cout << "rv!\n";
        return 1e4;
    }

    double dv_total = 0.0;
    for (int i = 1; i < impulseNum + 1; i++) {
        dv_total += std::sqrt(
                dvImpulse[i][0] * dvImpulse[i][0] + dvImpulse[i][1] * dvImpulse[i][1] +
                dvImpulse[i][2] * dvImpulse[i][2]);
    }

    return dv_total;
}

//double MultiImpulseTransfer::transform_t(const std::vector<double> &x) {
//    tImpulse.clear();
//    tImpulse.resize(impulseNum + 1);
//    dvImpulse.clear();
//    dvImpulse.resize(impulseNum + 1);
//    rvImpulse.clear();
//    rvImpulse.resize(impulseNum + 1);
//
//    double endTime_new = endTime + 3600.0 * (x[4 * impulseNum - 7] - 0.5);
//
//
//    double t0 = startTime + (endTime - startTime) * x[0];
//    double rv[6];
//    for (int i = 0; i < 6; i++) {
//        rv[i] = rv0_array[i];
//    }
//    double rvf[6];
//    //target0.getRv(t0, rv);
//    if (!isBiggerThanRmin(rv, MU_EARTH)) {
//        //std::cout << "rv0!\n";
//        return 1e4;
//    }
//
//    propagateTrajectory(startTime, t0, rv, rvf);
//    std::copy_n(rvf, 6, rv);
//    if (!isBiggerThanRmin(rv, MU_EARTH)) {
//        //std::cout << "rv0!\n";
//        return 1e4;
//    }
//
//    for (int i = 0; i < 3; i++)
//        dvImpulse[0][i] = 0.0;
//
//
//    tImpulse[0] = startTime;
//    tImpulse[1] = t0;
//
//    for (int i = 0; i < 6; i++)
//        rvImpulse[0][i] = rv0_array[i];
//
//    // 深空机动
//    for (int i = 1; i < impulseNum - 1; i++) {
//        tImpulse[i + 1] = tImpulse[i] + (endTime - tImpulse[i]) * x[4 * i - 3];
//        double dvNorm = dvMaxSingle * x[4 * i - 2];
//        double phi = DPI * (x[4 * i - 1] - 0.5);
//        double theta = D2PI * x[4 * i];
//
//        dvImpulse[i][0] = std::cos(phi) * std::cos(theta) * dvNorm;
//        dvImpulse[i][1] = std::cos(phi) * std::sin(theta) * dvNorm;
//        dvImpulse[i][2] = std::sin(phi) * dvNorm;
//        for (int j = 0; j < 3; j++) {
//            rv[j + 3] += dvImpulse[i][j];
//        }
//        if (!isBiggerThanRmin(rv, MU_EARTH)) {
//            //std::cout << "rv!\n";
//            return 1e4;
//        }
//        std::copy_n(rv, 6, rvImpulse[i].begin());
//
//        propagateTrajectory(tImpulse[i], tImpulse[i + 1], rv, rvf);
//        std::copy_n(rvf, 6, rv);
//    }
//
//    // 解Lambert问题
//    tImpulse[impulseNum] = endTime;
//    for (int i = 0; i < 6; i++) {
//        rvf[i] = rvf_array[i];
//    }
//    //target1.getRv(tImpulse[impulseNum], rvf);
//
//
//    double tof = tImpulse[impulseNum] - tImpulse[impulseNum - 1];
//    double v1[3], v2[3];
//    int lambert_flag;
//    solveLambert(tof, rv, rvf, v1, v2, lambert_flag, MU_EARTH);
//    if (lambert_flag > 0)
//        return 1e5;
//    for (int i = 0; i < 3; i++) {
//        dvImpulse[impulseNum - 1][i] = v1[i] - rv[i + 3];
//        dvImpulse[impulseNum][i] = rvf[i + 3] - v2[i];
//    }
//    for (int i = 0; i < 3; i++) {
//        rvImpulse[impulseNum - 1][i] = rv[i];
//        rvImpulse[impulseNum - 1][i + 3] = v1[i];
//        rvImpulse[impulseNum][i] = rvf[i];
//        rvImpulse[impulseNum][i + 3] = rvf[i + 3];
//    }
//    if (!isBiggerThanRmin(rvImpulse[impulseNum - 1].data(), MU_EARTH)) {
//        //std::cout << "rv!\n";
//        return 1e4;
//    }
//
//    double dv_total = 0.0;
//    for (int i = 1; i < impulseNum + 1; i++) {
//        dv_total += std::sqrt(
//                dvImpulse[i][0] * dvImpulse[i][0] + dvImpulse[i][1] * dvImpulse[i][1] +
//                dvImpulse[i][2] * dvImpulse[i][2]);
//    }
//
//    return dv_total;
//}

void MultiImpulseTransfer::MITResult::clear() {
    dvTotal = 1e10;
    success = false;
//    dRAAN_out.clear();
//    tf_out.clear();
    timeList.clear();
    rvList.clear();
    dvList.clear();
}

/**
 * @brief 如果在一个脉冲序列中两个节点时间间隔过短，则合并之
 * @param tol 两节点视为同一个的时间间隔阈值
 */
void MultiImpulseTransfer::MITResult::mergeNode(double tol) {
    if (!success)
        return;

    std::vector<double> timeList_new;
    std::vector<std::array<double, 6>> rvList_new; // 包含了dv
    std::vector<std::array<double, 3>> dvList_new;

    for (size_t i = 0; i < timeList.size(); i++) {
        if (i == 0 || std::abs(timeList[i] - timeList[i - 1]) >= tol) {
            timeList_new.push_back(timeList[i]);
            rvList_new.push_back(rvList[i]);
            dvList_new.push_back(dvList[i]);
        } else {
            // 合并节点
            timeList_new.back() = timeList[i];
            rvList_new.back() = rvList[i];
            dvList_new.back()[0] += dvList[i][0];
            dvList_new.back()[1] += dvList[i][1];
            dvList_new.back()[2] += dvList[i][2];
        }
    }

    timeList = timeList_new;
    rvList = rvList_new;
    dvList = dvList_new;

    dvTotal = 0.0;
    for (const auto &dv: dvList) {
        dvTotal += std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
    }
}

void MultiImpulseTransfer::getResult(const std::vector<double> &x, MITResult &result) {
    result.clear();
    result.dvTotal = transform(x);
    result.success = true;
    if (result.dvTotal > 0.9e10)
        result.success = false;
    result.startId = target0.id;
    result.endId = target1.id;
    result.timeList.resize(impulseNum);
    result.rvList.resize(impulseNum);
    result.dvList.resize(impulseNum);
    std::copy_n(tImpulse.begin() + 1, impulseNum, result.timeList.begin());
    std::copy_n(rvImpulse.begin() + 1, impulseNum, result.rvList.begin());
    std::copy_n(dvImpulse.begin() + 1, impulseNum, result.dvList.begin());
    // result.timeList = tImpulse;
    // result.rvList = rvImpulse;
    // result.dvList = dvImpulse;
}

void MultiImpulseTransfer::getResult_short(const std::vector<double> &x, MITResult &result) {
    result.clear();
    result.dvTotal = transform_short(x);
    result.success = true;
    if (result.dvTotal > 0.9e10)
        result.success = false;
    result.startId = target0.id;
    result.endId = target1.id;
    result.timeList.resize(impulseNum);
    result.rvList.resize(impulseNum);
    result.dvList.resize(impulseNum);
    std::copy_n(tImpulse.begin() + 1, impulseNum, result.timeList.begin());
    std::copy_n(rvImpulse.begin() + 1, impulseNum, result.rvList.begin());
    std::copy_n(dvImpulse.begin() + 1, impulseNum, result.dvList.begin());
    result.tf_out = tf_out;
    result.dRAAN_out = dRAAN_out;
    // result.timeList = tImpulse;
    // result.rvList = rvImpulse;
    // result.dvList = dvImpulse;
}

double MultiImpulseTransfer::transform_short(const std::vector<double> &x) {

    tImpulse.clear();
    tImpulse.resize(impulseNum + 1);
    dvImpulse.clear();
    dvImpulse.resize(impulseNum + 1);
    rvImpulse.clear();
    rvImpulse.resize(impulseNum + 1);

    double endTime_new = endTime + 60 * 60.0 * (x[4 * impulseNum - 7] - 0.5);
    double rvf_[6];
    double tf_;
    double dOmega = dRAAN + (x[4 * impulseNum - 6] - 0.5) * 0.1 / 180.0 * DPI;

    double t1_ = endTime_new;
    InPlane_pso(P1, S1, t1_, &tf_, rvf_, 1, ratioFlag, dOmega, dH);
    endTime = tf_;


    double t0 = startTime + (endTime - startTime) * x[0];
    double rv[6];
    for (int i = 0; i < 6; i++) {
        rv[i] = rv0_array[i];
    }
    double rvf[6];
    //target0.getRv(t0, rv);
    if (!isBiggerThanRmin(rv, MU_EARTH)) {
        //std::cout << "rv0!\n";
        return 1e4;
    }

    propagateTrajectory(startTime, t0, rv, rvf);
    std::copy_n(rvf, 6, rv);
    if (!isBiggerThanRmin(rv, MU_EARTH)) {
        //std::cout << "rv0!\n";
        return 1e4;
    }

    for (int i = 0; i < 3; i++)
        dvImpulse[0][i] = 0.0;

    tImpulse[0] = startTime;
    tImpulse[1] = t0;

    for (int i = 0; i < 6; i++)
        rvImpulse[0][i] = rv0_array[i];

    // 深空机动
    for (int i = 1; i < impulseNum - 1; i++) {
        tImpulse[i + 1] = tImpulse[i] + (endTime - tImpulse[i]) * x[4 * i - 3];
        double dvNorm = dvMaxSingle * x[4 * i - 2];
        double phi = DPI * (x[4 * i - 1] - 0.5);
        double theta = D2PI * x[4 * i];

        dvImpulse[i][0] = std::cos(phi) * std::cos(theta) * dvNorm;
        dvImpulse[i][1] = std::cos(phi) * std::sin(theta) * dvNorm;
        dvImpulse[i][2] = std::sin(phi) * dvNorm;
        for (int j = 0; j < 3; j++) {
            rv[j + 3] += dvImpulse[i][j];
        }
        if (!isBiggerThanRmin(rv, MU_EARTH)) {
            //std::cout << "rv!\n";
            return 1e4;
        }
        std::copy_n(rv, 6, rvImpulse[i].begin());

        propagateTrajectory(tImpulse[i], tImpulse[i + 1], rv, rvf);
        std::copy_n(rvf, 6, rv);
    }

    // 解Lambert问题
    tImpulse[impulseNum] = endTime;
    for (int i = 0; i < 6; i++) {
        rvf[i] = rvf_[i];
    }
    //target1.getRv(tImpulse[impulseNum], rvf);


    double tof = tImpulse[impulseNum] - tImpulse[impulseNum - 1];
    double v1[3], v2[3];
    int lambert_flag;
    solveLambert(tof, rv, rvf, v1, v2, lambert_flag, MU_EARTH);
    if (lambert_flag > 0)
        return 1e5;
    for (int i = 0; i < 3; i++) {
        dvImpulse[impulseNum - 1][i] = v1[i] - rv[i + 3];
        dvImpulse[impulseNum][i] = rvf[i + 3] - v2[i];
    }
    for (int i = 0; i < 3; i++) {
        rvImpulse[impulseNum - 1][i] = rv[i];
        rvImpulse[impulseNum - 1][i + 3] = v1[i];
        rvImpulse[impulseNum][i] = rvf[i];
        rvImpulse[impulseNum][i + 3] = rvf[i + 3];
    }
    if (!isBiggerThanRmin(rvImpulse[impulseNum - 1].data(), MU_EARTH)) {
        //std::cout << "rv!\n";
        return 1e4;
    }

    double dv_total = 0.0;
    for (int i = 1; i < impulseNum + 1; i++) {
        dv_total += std::sqrt(
                dvImpulse[i][0] * dvImpulse[i][0] + dvImpulse[i][1] * dvImpulse[i][1] +
                dvImpulse[i][2] * dvImpulse[i][2]);
    }

    tf_out = endTime_new;
    dRAAN_out = dOmega;

    return dv_total;

    return 0;
}

double MITOAdapter(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    MultiImpulseTransfer mit(*(static_cast<MultiImpulseTransfer *>(data)));
    daceInitializeThread();
    const auto dv = mit.transform(x);
    daceCleanupThread();
    return dv;

    //return mit.transform(x);
}

double MITOAdapter_short(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    MultiImpulseTransfer mit(*(static_cast<MultiImpulseTransfer *>(data)));
    daceInitializeThread();
    const auto dv = mit.transform_short(x);
    daceCleanupThread();
    return dv;
    //return mit.transform(x);
}