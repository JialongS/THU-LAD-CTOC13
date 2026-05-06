#ifndef MULTIIMPULSETRANSFER_H
#define MULTIIMPULSETRANSFER_H

/**
 * MultiImpulseTransfer.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-7-31
 */

#include <vector>
#include <array>


class Target {
public:
    Target() = default;

    explicit Target(int id_);

    Target(const Target &tar);

    ~Target() = default;

    void getRv(double t, double *rv) const;

    int id;
};

enum class Dynamics {
    ACCURACY, APPROXIMATE
};

/**
 * @brief 多脉冲轨道交会转移优化
 */
class MultiImpulseTransfer {
public:
    MultiImpulseTransfer() = default;

    MultiImpulseTransfer(double startTime_, double endTime_, int targetId0, int targetId1, double dvMaxSingle_,
                         int impulseNum_, double rMin_, Dynamics dynamicsModel_);

    MultiImpulseTransfer(const MultiImpulseTransfer &mit);

    ~MultiImpulseTransfer() = default;

    [[nodiscard]] int calcNumVars() const;

    void propagateTrajectory(double t0, double tf, const double *rv0, double *rvf) const;

    void solveLambert(double tof, const double *rv1, const double *rv2, double *v1, double *v2, int &flag,
                      double mu) const;

    bool isBiggerThanRmin(const double *rv, double mu);

    double transform(const std::vector<double> &x); // 目标函数
    double transform_short(const std::vector<double> &x); // 用于短时间转移

    class MITResult {
        // 优化结果
    public:
        MITResult() = default;

        ~MITResult() = default;

        void clear();

        void mergeNode(double tol = 1e-8);

        double dvTotal{1e10};
        bool success{false};

        int startId{};
        int endId{};

        double tf_out{}, dRAAN_out{};
        std::vector<double> timeList;
        std::vector<std::array<double, 6>> rvList; // 包含了dv
        std::vector<std::array<double, 3>> dvList;
    };


    void getResult(const std::vector<double> &x, MITResult &result);

    void getResult_short(const std::vector<double> &x, MITResult &result);

    double startTime{}; // 约束的起始时间
    double endTime{}; // 约束的结束时间
    Target target0{}; // 出发点目标
    Target target1{}; // 交会目标
    double dvMaxSingle{}; // 单次机动允许的最大脉冲
    int impulseNum{}; // 机动脉冲数，包含出发和交会脉冲
    double rMin{};
    std::array<double, 6> rv0_array;
    std::array<double, 6> rvf_array;
    double dRAAN{};
    double dH{};
    int ratioFlag{};
    int P1{};
    int S1{};
//    double tf_out{};
//    double dRAAN_out{};
    Dynamics dynamicModel{Dynamics::ACCURACY};

private:
    double tf_out{};
    double dRAAN_out{};
    std::vector<double> tImpulse;
    std::vector<std::array<double, 3>> dvImpulse;
    std::vector<std::array<double, 6>> rvImpulse; // 包含了dv
};


double MITOAdapter(const std::vector<double> &x, std::vector<double> &grad, void *data);

double MITOAdapter_short(const std::vector<double> &x, std::vector<double> &grad, void *data);

#endif //MULTIIMPULSETRANSFER_H
