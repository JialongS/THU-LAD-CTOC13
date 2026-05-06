#include "constant.h"
#include "dynamics.h"
#include "inPlaneOptimization.h"
#include "J2Lambert.h"
#include "io.h"
#include "optimization.h"
#include "dace/dace_s.h"
#include <cmath>
#include <iostream>
#include <OrbitFun.h>
#include <VecMat.h>

InPlaneTransfer::InPlaneTransfer()
{
}

InPlaneTransfer::InPlaneTransfer(const InPlaneTransfer& ipt)
{
    target.n1 = ipt.target.n1;
    target.n2 = ipt.target.n2;
    target.n3 = ipt.target.n3;
    startTime = ipt.startTime;
    endTime = ipt.endTime;
    mode = ipt.mode;
    dynamics = ipt.dynamics;
    rv0 = ipt.rv0;
}

/**
* @brief 计算优化变量的个数
*/
int InPlaneTransfer::calcNumVars() const
{
    int num;
    if (mode == Mode::GEAR)
        num = 4;
    else if (mode == Mode::CRICLE)
        num = 5;
    else if (mode == Mode::TIME)
        num = 2;
    else
        num = -1;
    return num;
}

void InPlaneTransfer::propagateTrajectory(double t0, double tf, const double* rv0, double* rvf) const
{
    if (dynamics == Dynamics::ACCURACY)
    {
        int flag = orbitPropagateJ2(t0, tf, rv0, rvf);
    }
    else if (dynamics == Dynamics::APPROXIMATE)
    {
        orbitPropagateJ2Approximation(t0, tf, rv0, rvf);
    }
}

void InPlaneTransfer::solveLambert(double tof, const double* rv1, const double* rv2, double* v1, double* v2, int& flag,
                                   double mu) const
{
    if (mode == Mode::GEAR || mode == Mode::CRICLE)
        J2_Lambert_MRPLP(v1, v2, rv1, rv2, tof, flag, mu, 40, 1.0);
}

/**
* @brief 齿轮转移
* @param x - 优化变量: x[0] -> 机动时间, x[1] -> 飞越时间, x[2], x[3], x[4] -> 飞越位置
*/
double InPlaneTransfer::gearTransfer(const std::vector<double>& x)
{
    reset(2);
    //dvList.resize(3);

    double t0 = startTime; // + x[0] * (endTime - startTime);
    double t1 = t0 + x[0] * (endTime - t0);
    double delta = 1001.0 + x[1] * (49999.0 - 1001.0); // 距离
    double theta = x[2] * D2PI; // 经度
    double phi = -DPI / 2.0 + x[3] * DPI; // 纬度

    timeList[0] = t0;
    timeList[1] = t1;

    std::array<double, 6> rv_flyby{}, rv_tar;
    target.getrv(t1, rv_tar.data());
    rv_flyby[0] = rv_tar[0] + delta * cos(phi) * cos(theta);
    rv_flyby[1] = rv_tar[1] + delta * cos(phi) * sin(theta);
    rv_flyby[2] = rv_tar[2] + delta * sin(phi);

    // 积分
    std::array<double, 6> rv;
    propagateTrajectory(startTime, t0, rv0.data(), rv.data());

    // 解Lambert
    int lambert_flag;
    std::array<double, 3> v1, v2, dv;
    solveLambert(t1 - t0, rv.data(), rv_flyby.data(), v1.data(), v2.data(), lambert_flag, MU_EARTH);
    if (lambert_flag > 0)
        return 1e10;

    propagateTrajectory(t0, t1, rv.data(), rv_flyby.data());

    // for (int i = 0; i < 3; i++)
    // {
    //     rv_flyby[i + 3] = v2[i];
    // }

    for (int i = 0; i < 3; i++)
    {
        dv[i] = v1[i] - rv[i + 3];
        rv[i + 3] = v1[i];
    }
    dvList[0] = dv;
    rvList[0] = rv;
    rvList[1] = rv_flyby;
    eventList[0] = static_cast<int>(Event::PULSE);
    eventList[1] = static_cast<int>(Event::INSPECTION);
    targetList[0] = {0, 0, 0};
    targetList[1] = {target.n1, target.n2, target.n3};

    double J = 0.0;
    J += sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);

    for (int i = 0; i < 3; i++)
    {
        dv[i] = rv_tar[i + 3] - v2[i];
    }
    dvList[1] = dv;
    if (double err = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]); err > 149.99999)
        return 1e10 * (err - 149.99999);
    // if (double err = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]); err > 149.99999)
    // {
    //     for (int i = 0; i < 3; i++)
    //     {
    //         dv[i] = dv[i] * (err - 149.0) / err;
    //     }
    //     J += sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
    // }

    // 计算飞越速度在目标速度方向上的投影
    double dv_t = calcDvInDirection(rv_tar.data() + 3, dv.data());
    //J += fabs(dv_t - 140.0);

    return J;
}

/**
 *
 * @param x 优化变量: x[0] -> 机动时间, x[1] -> 飞越时间
 * @return 优化目标值
 */
// double InPlaneTransfer::timeTransfer(const std::vector<double>& x)
// {
//     reset(2);
//
//     double t0 = startTime + x[0] * (endTime - startTime);
//     double t1 = t0 + x[1] * (endTime - t0);
//
//     timeList[0] = t0;
//     timeList[1] = t1;
//
//     std::array<double, 6> rv_tar, rv_flyby;
//     target.getrv(t1, rv_tar.data());
//
//     // 积分
//     std::array<double, 6> rv;
//     propagateTrajectory(startTime, t0, rv0.data(), rv.data());
//
//     // 解Lambert
//     int lambert_flag;
//     std::array<double, 3> v1, v2, dv, r2_true;
//     J2_Lambert_MRPLP_best(v1.data(), v2.data(), rv.data(), rv_tar.data(), t1 - t0, r2_true, lambert_flag, MU_EARTH);
//     if (lambert_flag > 0)
//         return 1e10;
//
//     for (int i = 0; i < 3; i++)
//     {
//         rv_flyby[i] = r2_true[i];
//         rv_flyby[i + 3] = v2[i];
//     }
//
//     for (int i = 0; i < 3; i++)
//     {
//         dv[i] = v1[i] - rv[i + 3];
//         rv[i + 3] = v1[i];
//     }
//     dvList[0] = dv;
//     rvList[0] = rv;
//     rvList[1] = rv_flyby;
//     eventList[0] = static_cast<int>(Event::PULSE);
//     eventList[1] = static_cast<int>(Event::INSPECTION);
//     targetList[0] = {0, 0, 0};
//     targetList[1] = {target.n1, target.n2, target.n3};
//
//     double J = 0.0;
//     J += sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
//
//     for (int i = 0; i < 3; i++)
//     {
//         dv[i] = rv_tar[i + 3] - v2[i];
//     }
//     dvList[1] = dv;
//     if (sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]) > 149.99999)
//         return 1e10;
//
//     // 计算飞越速度在目标速度方向上的投影
//     double dv_t = calcDvInDirection(rv_tar.data() + 3, dv.data());
//     //J += fabs(dv_t - 140.0);
//
//     return J;
// }

double InPlaneTransfer::transfer(const std::vector<double>& x)
{
    double J = 0.0;
    if (mode == Mode::GEAR)
        J = gearTransfer(x);
    // if (mode == Mode::TIME)
    //     J = timeTransfer(x);

    return J;
}

void InPlaneTransfer::getResult(const std::vector<double>& x, IPTResult& result)
{
    transfer(x);
    result.timeList = timeList;
    result.rvList = rvList;
    result.dvList = dvList;
    for (int i = 0; i < 3; i++)
        result.dvList.back()[i] = 0.0;
    result.eventList = eventList;
    result.targetList = targetList;

    result.dv_total = 0;
    for (int i = 0; i < timeList.size(); i++)
        result.dv_total += sqrt(
            result.dvList[i][0] * result.dvList[i][0] + result.dvList[i][1] * result.dvList[i][1] + result.dvList[i][2]
            * result.dvList[i][2]);
    result.v_flyby = sqrt(
        dvList.back()[0] * dvList.back()[0] + dvList.back()[1] * dvList.back()[1] + dvList.back()[2] * dvList.back()[
            2]);
}

/**
* @brief 缩放rv
*/
void InPlaneTransfer::scaleRV(const double* rv0, double dr, double dv, double* rv)
{
    double r0 = sqrt(rv0[0] * rv0[0] + rv0[1] * rv0[1] + rv0[2] * rv0[2]);
    double r = r0 + dr;
    double v0 = sqrt(rv0[3] * rv0[3] + rv0[4] * rv0[4] + rv0[5] * rv0[5]);
    double v = v0 + dv;
    for (int i = 0; i < 3; i++)
    {
        rv[i] = rv0[i] * r / r0;
        rv[i + 3] = rv0[i + 3] * v / v0;
    }
}

/**
* @brief 在同一轨道面附近，具有dv速度差的物体相遇的时间的二倍
*/
double InPlaneTransfer::estimateTime(const double* rv1, const double* rv2, double dv)
{
    double r1 = sqrt(rv1[0] * rv1[0] + rv1[1] * rv1[1] + rv1[2] * rv1[2]);
    double r2 = sqrt(rv2[0] * rv2[0] + rv2[1] * rv2[1] + rv2[2] * rv2[2]);
    double dot = rv1[0] * rv2[0] + rv1[1] * rv2[1] + rv1[2] * rv2[2];
    double theta = acos(dot / r1 / r2);
    double time = 2.0 * theta * (r1 > r2 ? r1 : r2) / dv;
    return time;
}

/**
* @brief 计算dv在v1方向上的投影大小
*/
double InPlaneTransfer::calcDvInDirection(const double* v1, const double* dv)
{
    double dot = dv[0] * v1[0] + dv[1] * v1[1] + dv[2] * v1[2];
    double norm_v1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);

    return dot / norm_v1;
}

/**
 * @brief 估计当前卫星和目标间能达到的最小距离
 * @param t0 当前时间 s
 * @param rv_chase 追赶卫星的位置速度
 * @param target 目标
 * @param maxIter 最大迭代次数
 * @param rTol 收敛阈值
 * @param vTol
 * @return [flag(0 for success), time_meet(s), r_min(m)]
 */
std::tuple<int, double, double, double> InPlaneTransfer::estimateNearestDistance(double t0, const double* rv_chase,
    const Satellite& target, int maxIter,
    double rTol, double vTol)
{
    int flag = 1;
    double t = t0;
    double t_meet = t;
    double r_min_g = 1e10;
    double r_min = 1e10;
    double v_flyby = 1e10;
    double v_flyby_g = 1e10;
    double t_g = t;

    double rv_c[6];
    V_Copy(rv_c, rv_chase, 6);

    for (int i = 0; i < maxIter; i++)
    {
        double rv_t[6];
        target.getrv(t, rv_t);
        t_meet += estimateMeetingTime(rv_c, rv_t, t - t0);
        double rv_temp[6];
        orbitPropagateJ2(t, t_meet, rv_c, rv_temp);
        t = t_meet;
        V_Copy(rv_c, rv_temp, 6);
        target.getrv(t_meet, rv_t);
        r_min = sqrt(
            (rv_c[0] - rv_t[0]) * (rv_c[0] - rv_t[0]) + (rv_c[1] - rv_t[1]) * (rv_c[1] - rv_t[1]) + (rv_c[2] - rv_t[2])
            * (rv_c[2] - rv_t[2]));
        v_flyby = sqrt(
            (rv_c[3] - rv_t[3]) * (rv_c[3] - rv_t[3]) + (rv_c[4] - rv_t[4]) * (rv_c[4] - rv_t[4]) + (rv_c[5] - rv_t[5])
            * (rv_c[5] - rv_t[5]));
        if (r_min < r_min_g)
        {
            r_min_g = r_min;
            v_flyby_g = v_flyby;
            t_g = t;
        }
        if (r_min <= rTol && v_flyby <= vTol)
        {
            r_min_g = r_min;
            v_flyby_g = v_flyby;
            t_g = t;
            break;
        }
    }

    if (r_min_g <= rTol && v_flyby_g <= vTol)
        flag = 0; // success

    return std::make_tuple(flag, t_g, r_min_g, v_flyby_g);
}

std::tuple<int, double, double, double> InPlaneTransfer::estimateNearestDistance_Discrete(double t0,
    const double* rv_chase, const Satellite& target, double t_step, double t_end, double rTol, double vTol)
{
    int flag = 1;
    double t = t0;
    double t_g = t;
    double r_min = 1e10;
    double v_flyby = 1e10;
    double r_min_g = 1e10;
    double v_flyby_g = 1e10;
    int N = ceil((t_end - t0) / t_step);
    double rv_c[6];
    V_Copy(rv_c, rv_chase, 6);

    for (int i = 0; i < N; i++)
    {
        t = t0 + (i + 1) * t_step;
        double rv_t[6], rv_temp[6];
        target.getrv(t, rv_t);
        orbitPropagateJ2(t, t - t_step, rv_c, rv_temp);
        V_Copy(rv_c, rv_temp, 6);
        r_min = sqrt(
            (rv_c[0] - rv_t[0]) * (rv_c[0] - rv_t[0]) + (rv_c[1] - rv_t[1]) * (rv_c[1] - rv_t[1]) + (rv_c[2] - rv_t[2])
            * (rv_c[2] - rv_t[2]));
        v_flyby = sqrt(
            (rv_c[3] - rv_t[3]) * (rv_c[3] - rv_t[3]) + (rv_c[4] - rv_t[4]) * (rv_c[4] - rv_t[4]) + (rv_c[5] - rv_t[5])
            * (rv_c[5] - rv_t[5]));
        if (r_min < r_min_g)
        {
            r_min_g = r_min;
            v_flyby_g = v_flyby;
            t_g = t;
        }
        if (r_min <= rTol && v_flyby <= vTol)
        {
            r_min_g = r_min;
            v_flyby_g = v_flyby;
            t_g = t;
            break;
        }
    }

    if (r_min_g <= rTol && v_flyby_g <= vTol)
        flag = 0; // success

    return std::make_tuple(flag, t_g, r_min_g, v_flyby_g);
}

double InPlaneTransfer::estimateMeetingTime(const double* rv_chase, const double* rv_target, double t_span)
{
    double coe_chase[6], coe_target[6];
    int flag;
    rv2coe(flag, coe_chase, rv_chase, MU_EARTH);
    rv2coe(flag, coe_target, rv_target, MU_EARTH);
    const double n_chase = sqrt(MU_EARTH / coe_chase[0] / coe_chase[0] / coe_chase[0]);
    const double n_target = sqrt(MU_EARTH / coe_target[0] / coe_target[0] / coe_target[0]);

    double M_chase = f2E(flag, coe_chase[5] + coe_chase[4], coe_chase[1]);
    double M_target = f2E(flag, coe_chase[5] + coe_chase[4], coe_target[1]);
    M_chase = E2M(flag, M_chase, coe_chase[1]);
    M_target = E2M(flag, M_target, coe_target[1]);
    M_chase = normalizeAngle(M_chase);
    M_target = normalizeAngle(M_target);

    double t = (M_target - M_chase) / (n_chase - n_target);
    if (t < 0 && fabs(t) > t_span)
        t += D2PI / fabs(n_chase - n_target);

    return t;
}

double InPlaneTransfer::normalizeAngle(double radians)
{
    double twoPi = 2.0 * M_PI;
    double normalized = fmod(radians, twoPi);
    if (normalized < 0)
    {
        normalized += twoPi;
    }

    return normalized;
}

void InPlaneTransfer::reset(int n)
{
    timeList.clear();
    dvList.clear();
    rvList.clear();
    eventList.clear();
    targetList.clear();

    timeList.resize(n);
    dvList.resize(n);
    rvList.resize(n);
    eventList.resize(n);
    targetList.resize(n);
}

void InPlaneTransfer::IPTResult::clear()
{
    rvList.clear();
    targetList.clear();
    timeList.clear();
    eventList.clear();
    dvList.clear();
}

void InPlaneTransfer::IPTResult::print()
{
    dv_total = 0.0;
    for (int i = 0; i < timeList.size(); i++)
    {
        dv_total += sqrt(dvList[i][0] * dvList[i][0] + dvList[i][1] * dvList[i][1] + dvList[i][2] * dvList[i][2]);
    }
    std::cout << "dv_total: " << dv_total << " m/s" << std::endl;
    //std::cout << "flyby speed: " << v_flyby << " m/s" << std::endl;
    std::cout << "tof: " << (timeList.back() - timeList.front()) / 3600.0 << " h" << std::endl;
}

void InPlaneTransfer::IPTResult::save(std::fstream& fout, int& lineID, int lauchID, int satID) const
{
    for (int i = 0; i < timeList.size(); i++)
    {
        // double rv_before[6];
        // for (int iter = 0; iter < 3; iter++)
        // {
        //     rv_before[iter] = rvList[i][iter];
        //     rv_before[iter + 3] = rvList[i][iter + 3] - dvList[i][iter + 3];
        // }
        printCTOC13(fout, lineID, lauchID, satID, eventList[i], timeList[i], targetList[i][0], targetList[i][1],
                    targetList[i][2], rvList[i].data(), dvList[i].data());
        lineID++;
    }
}

/**
* @brief InPlaneTransfer类的包装器函数
*/
double IPOAdapter(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    InPlaneTransfer ipt(*(static_cast<InPlaneTransfer*>(data)));
    daceInitializeThread();
    const auto dv = ipt.transfer(x);
    daceCleanupThread();

    return dv;
}

void IPTOptimization::run(std::fstream& fout)
{
    InPlaneTransfer ipt;
    ipt.startTime = t0;
    ipt.mode = InPlaneTransfer::Mode::GEAR;
    ipt.dynamics = InPlaneTransfer::Dynamics::ACCURACY;
    Satellite target(satelliteSequence[0][0], satelliteSequence[0][1], satelliteSequence[0][2]);
    double rv_temp[6];
    ipt.rv0 = rv0;
    dv = 0.0;

    for (int i = 1; i < satelliteSequence.size(); i++)
    {
        ipt.target = Satellite(satelliteSequence[i][0], satelliteSequence[i][1], satelliteSequence[i][2]);
        ipt.target.getrv(ipt.startTime, rv_temp);
        ipt.endTime = ipt.startTime + InPlaneTransfer::estimateTime(ipt.rv0.data(), rv_temp, 150.0);
        InPlaneTransfer::IPTResult result;

        // auto [flag,t_meet,r_min, v_flyby] = InPlaneTransfer::estimateNearestDistance(
        //     ipt.startTime, ipt.rv0.data(), ipt.target);
        // auto [flag,t_meet,r_min, v_flyby] = InPlaneTransfer::estimateNearestDistance_Discrete(
        //     ipt.startTime, ipt.rv0.data(), ipt.target, 10.0, ipt.endTime);
        //
        // std::cout << "r_min: " << r_min << " m\n" << "v_flyby: " << v_flyby << " m/s" << std::endl;
        //
        // if (flag)
        {
            std::vector<double> x(ipt.calcNumVars());
            double fbest;
            DE(IPOAdapter, &ipt, x, fbest, ipt.calcNumVars(), 10 * ipt.calcNumVars(), 1000, 500);
            ipt.getResult(x, result);
        }
        // else
        // {
        //     result.dv_total = 0.0;
        //     result.v_flyby = v_flyby;
        //     result.timeList = {t_meet};
        //     std::array<double, 6> rv_meet{};
        //     orbitPropagateJ2(ipt.startTime, t_meet, ipt.rv0.data(), rv_meet.data());
        //     result.rvList.clear();
        //     result.rvList.push_back(rv_meet);
        //     result.dvList.clear();
        //     result.dvList.push_back({0.0, 0.0, 0.0});
        //     result.eventList = {static_cast<int>(Event::INSPECTION)};\
        //     result.targetList.clear();
        //     result.targetList.push_back({ipt.target.n1, ipt.target.n2, ipt.target.n3});
        // }

        result.print();
        result.save(fout, lineID, launchID, satID);
        dv += result.dv_total;
        ipt.startTime = result.timeList.back();
        ipt.rv0 = result.rvList.back();
    }
    std::cout << "dv: " << dv << " m/s" << std::endl;
    std::cout << "time: " << (ipt.startTime - t0) / 86400.0 << " day" << std::endl;
}
