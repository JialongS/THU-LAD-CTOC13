#pragma once

/**
 * inPlaneOptimization.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-8-22
 */

#include "data_ctoc13.h"
#include <vector>
#include <array>
#include <fstream>

/**
* @brief 目标卫星结构体
*/
class Satellite
{
public:
    inline Satellite()
    = default;

    inline Satellite(int n1_, int n2_, int n3_) : n1(n1_), n2(n2_), n3(n3_)
    {
    };

    inline Satellite(const Satellite& sat) : n1(sat.n1), n2(sat.n2), n3(sat.n3)
    {
    };

    inline Satellite& operator=(const Satellite& rhs)
    {
        if (this != &rhs)
        {
            n1 = rhs.n1;
            n2 = rhs.n2;
            n3 = rhs.n3;
        }
        return *this;
    }

    int n1{};
    int n2{};
    int n3{};

    inline void getrv(double t, double* rv) const
    {
        get_constellation_RV(n1, n2, n3, t, rv);
    }
};

enum class Event { LANUCH = 0, PULSE = 1, INSPECTION = 2, END = 3 };

/**
* @brief 单目标飞越
*/
class InPlaneTransfer
{
public:
    enum class Mode { GEAR, CRICLE, TIME };

    enum class Dynamics { ACCURACY, APPROXIMATE };

    class IPTResult
    {
    public:
        double dv_total{};
        double v_flyby{};
        std::vector<double> timeList;
        std::vector<std::array<double, 3>> dvList; // 不包含飞越速度
        std::vector<std::array<double, 6>> rvList; // 不包含dv
        std::vector<int> eventList;
        std::vector<std::array<int, 3>> targetList;

        void clear();
        void print();
        void save(std::fstream& fout, int& lineID, int lauchID, int satID) const;
    };

    Satellite target;
    double startTime{};
    double endTime{};
    Mode mode{Mode::GEAR};
    Dynamics dynamics{Dynamics::ACCURACY};
    std::array<double, 6> rv0;

    InPlaneTransfer();
    InPlaneTransfer(const InPlaneTransfer& ipt);

    [[nodiscard]] int calcNumVars() const;
    void propagateTrajectory(double t0, double tf, const double* rv0, double* rvf) const;
    void solveLambert(double tof, const double* rv1, const double* rv2, double* v1, double* v2, int& flag,
                      double mu) const;
    double gearTransfer(const std::vector<double>& x);
    //double timeTransfer(const std::vector<double>& x);
    double transfer(const std::vector<double>& x);
    void getResult(const std::vector<double>& x, IPTResult& result);

    static void scaleRV(const double* rv0, double dr, double dv, double* rv);
    static double estimateTime(const double* rv1, const double* rv2, double dv);
    static double calcDvInDirection(const double* v1, const double* dv);
    static std::tuple<int, double, double, double> estimateNearestDistance(double t0, const double* rv_chase,
                                                                           const Satellite& target,
                                                                           int maxIter = 10, double rTol = 50e3,
                                                                           double vTol = 150.0);
    static std::tuple<int, double, double, double> estimateNearestDistance_Discrete(double t0, const double* rv_chase,
        const Satellite& target, double t_step, double t_end,
        double rTol = 50e3,
        double vTol = 150.0);
    static double estimateMeetingTime(const double* rv_chase, const double* rv_target, double t_span);
    static double normalizeAngle(double radians);

private:
    std::vector<double> timeList;
    std::vector<std::array<double, 3>> dvList;
    std::vector<std::array<double, 6>> rvList; // 包含了dv
    std::vector<int> eventList; // 包含了飞越速度
    std::vector<std::array<int, 3>> targetList;
    void reset(int n);
};

double IPOAdapter(const std::vector<double>& x, std::vector<double>& grad, void* data);

/**
* @brief 飞越同轨道面优化
*/
class IPTOptimization
{
public:
    std::vector<std::array<int, 3>> satelliteSequence;
    double t0{};
    int lineID{1};
    int launchID{1};
    int satID{1};
    double dv{};
    std::array<double, 6> rv0{};

    void run(std::fstream& fout);
};
