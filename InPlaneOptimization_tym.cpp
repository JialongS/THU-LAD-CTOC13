#include "constant.h"
#include "dynamics.h"
#include "inPlaneOptimization_tym.h"
#include "J2Lambert.h"
#include "io.h"
#include "optimization.h"
#include <cmath>
#include <iostream>

void Constraint_InPlane_S2S_1(unsigned m, double* result, unsigned n, const double* X, double* grad, void* f_data)
{
    //优化约束函数，初始时s0在一号星附近
    //输入f_data：0~5初值，6~7飞行时间初值,8为t0时刻,9为星座编号,10为星座轨道面编号,11~13为目标三星编号
    //输入X：0~5初值位置速度改变值，6~7积分时间改变值

    double* para = (double*)f_data; //类型强制转换
    double s0[6];
    for (int i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            s0[i] = para[i] + (2.0 * X[i] - 1.0) * 50e3; //50km伸缩
        }
        else
        {
            s0[i] = para[i] + (2.0 * X[i] - 1.0) * 300.0; //500m/s伸缩
        }
    }
    double tfa = para[6] + (2.0 * X[6] - 1.0) * 20.0; //30s伸缩
    double tfb = para[7] + (2.0 * X[7] - 1.0) * 20.0;
    double t0 = para[8];
    int flag1, flag2;
    double rvf1[6], rvf2[6];
    flag1 = orbitPropagateJ2(t0, t0 + tfa, s0, rvf1, 1.0, 1e-13, 1e-13); //第一段，飞向二号星
    flag2 = orbitPropagateJ2(t0, t0 + tfa + tfb, s0, rvf2, 1.0, 1e-13, 1e-13); //第二段，飞向三号星
    int n1 = para[9];
    int n2 = para[10];
    int k1 = para[11];
    int k2 = para[12];
    int k3 = para[13];

    double rv1[6], rv2[6], rv3[6];
    get_constellation_RV(n1, n2, k1, t0, rv1);
    get_constellation_RV(n1, n2, k2, t0 + tfa, rv2);
    get_constellation_RV(n1, n2, k3, t0 + tfa + tfb, rv3);

    double dr_1[3], dr_2[3], dr_3[3], dv_1[3], dv_2[3], dv_3[3];
    for (int i = 0; i < 3; ++i)
    {
        dr_1[i] = s0[i] - rv1[i];
        dr_2[i] = rvf1[i] - rv2[i];
        dr_3[i] = rvf2[i] - rv3[i];
        dv_1[i] = s0[i + 3] - rv1[i + 3];
        dv_2[i] = rvf1[i + 3] - rv2[i + 3];
        dv_3[i] = rvf2[i + 3] - rv3[i + 3];
    }

    double dr1 = sqrt(dr_1[0] * dr_1[0] + dr_1[1] * dr_1[1] + dr_1[2] * dr_1[2]);
    double dr2 = sqrt(dr_2[0] * dr_2[0] + dr_2[1] * dr_2[1] + dr_2[2] * dr_2[2]);
    double dr3 = sqrt(dr_3[0] * dr_3[0] + dr_3[1] * dr_3[1] + dr_3[2] * dr_3[2]);
    double dv1 = sqrt(dv_1[0] * dv_1[0] + dv_1[1] * dv_1[1] + dv_1[2] * dv_1[2]);
    double dv2 = sqrt(dv_2[0] * dv_2[0] + dv_2[1] * dv_2[1] + dv_2[2] * dv_2[2]);
    double dv3 = sqrt(dv_3[0] * dv_3[0] + dv_3[1] * dv_3[1] + dv_3[2] * dv_3[2]);
    result[0] = 1.0e3 - dr1;
    result[1] = dr1 - 50e3;
    result[2] = 1.0e3 - dr2;
    result[3] = dr2 - 50e3;
    result[4] = 1.0e3 - dr3;
    result[5] = dr3 - 50e3;
    result[6] = dv1 - 150.0;
    result[7] = dv2 - 150.0;
    result[8] = dv3 - 150.0;
    /*   std::cout << " NLOPT优化不等式约束：" << std::endl;
       for (size_t i = 0; i < 9; i++)
       {
           std::cout << result[i] << std::endl;
       }*/

    return;
}

void Constraint_InPlane_S2S_2(unsigned m, double* result, unsigned n, const double* X, double* grad, void* f_data)
{
    //优化约束函数
    //输入f_data：0~5初值，6~7飞行时间初值,8为t0时刻,9为星座编号,10为星座轨道面编号,11~13为目标三星编号
    //输入X：0~2初值位置速度改变值，3~4积分时间改变值

    double* para = (double*)f_data; //类型强制转换
    double s0[6];
    for (int i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            s0[i] = para[i];
        }
        else
        {
            s0[i] = para[i] + (2.0 * X[i - 3] - 1.0) * 300.0; //500m/s伸缩
        }
    }
    double tfa = para[6] + (2.0 * X[3] - 1.0) * 30.0; //30s伸缩
    double tfb = para[7] + (2.0 * X[4] - 1.0) * 30.0;
    double t0 = para[8];
    int flag1, flag2;
    double rvf1[6], rvf2[6];
    flag1 = orbitPropagateJ2(t0, t0 + tfa, s0, rvf1, 1.0, 1e-13, 1e-13); //第一段，飞向二号星
    flag2 = orbitPropagateJ2(t0, t0 + tfa + tfb, s0, rvf2, 1.0, 1e-13, 1e-13); //第二段，飞向三号星
    int n1 = para[9];
    int n2 = para[10];
    int k1 = para[11];
    int k2 = para[12];
    int k3 = para[13];

    double rv1[6], rv2[6], rv3[6];
    get_constellation_RV(n1, n2, k1, t0, rv1);
    get_constellation_RV(n1, n2, k2, t0 + tfa, rv2);
    get_constellation_RV(n1, n2, k3, t0 + tfa + tfb, rv3);

    double dr_1[3], dr_2[3], dr_3[3], dv_1[3], dv_2[3], dv_3[3];
    for (int i = 0; i < 3; ++i)
    {
        dr_2[i] = rvf1[i] - rv2[i];
        dr_3[i] = rvf2[i] - rv3[i];
        dv_2[i] = rvf1[i + 3] - rv2[i + 3];
        dv_3[i] = rvf2[i + 3] - rv3[i + 3];
    }

    double dr2 = sqrt(dr_2[0] * dr_2[0] + dr_2[1] * dr_2[1] + dr_2[2] * dr_2[2]);
    double dr3 = sqrt(dr_3[0] * dr_3[0] + dr_3[1] * dr_3[1] + dr_3[2] * dr_3[2]);
    double dv2 = sqrt(dv_2[0] * dv_2[0] + dv_2[1] * dv_2[1] + dv_2[2] * dv_2[2]);
    double dv3 = sqrt(dv_3[0] * dv_3[0] + dv_3[1] * dv_3[1] + dv_3[2] * dv_3[2]);
    result[0] = 1.0e3 - dr2;
    result[1] = dr2 - 50e3;
    result[2] = 1.0e3 - dr3;
    result[3] = dr3 - 50e3;
    result[4] = dv2 - 150.0;
    result[5] = dv3 - 150.0;
    return;
}

double Function_InPlane_S2S_1(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
    //优化目标函数，初始时s0在一号星附近
    //输入f_data：0~5初值，6~7飞行时间初值,8为t0时刻,9为星座编号,10为星座轨道面编号,11~13为目标三星编号
    //输入X：0~5初值位置速度改变值，6~7积分时间改变值

    double* para = (double*)f_data; //类型强制转换
    double s0[6];
    double r0[3];
    for (int i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            s0[i] = para[i] + (2.0 * X[i] - 1.0) * 50e3; //50km伸缩
            r0[i] = s0[i];
        }
        else
        {
            s0[i] = para[i] + (2.0 * X[i] - 1.0) * 300.0; //500m/s伸缩
        }
    }
    double tfa = para[6] + (2.0 * X[6] - 1.0) * 20.0; //60s伸缩
    double tfb = para[7] + (2.0 * X[7] - 1.0) * 20.0;
    double t0 = para[8];
    int flag1;
    double rvf1[6], rvf2[6];
    flag1 = orbitPropagateJ2(t0, t0 + tfa, s0, rvf1, 1.0, 1e-13, 1e-13); //第一段，飞向二号星
    int n1 = para[9];
    int n2 = para[10];
    int k1 = para[11];
    int k2 = para[12];

    double rv1[6], rv2[6];
    get_constellation_RV(n1, n2, k1, t0, rv1);
    get_constellation_RV(n1, n2, k2, t0 + tfa, rv2);

    double r_1[3], r_2[3], rf_1[3];
    for (int i = 0; i < 3; ++i)
    {
        r_1[i] = rv1[i];
        r_2[i] = rv2[i];
        rf_1[i] = rvf1[i];
    }
    double cosphi1 = cosphi_obs(t0, r_1, r0);
    double cosphi2 = cosphi_obs(t0, r_2, rf_1);
    double J = -cosphi1 * 10.0 - cosphi2 * 10.0;

    //double result[9];
    //Constraint_InPlane_S2S_1(1, result, 1, X.data(), grad.data(), f_data);
    //for (int i = 0; i < 9; ++i) {
    //    J += fmax(result[i], 0) * fmax(result[i], 0) * 1e3;
    //}

    return J;
}

double Function_InPlane_S2S_2(const std::vector<double>& X, std::vector<double>& grad, void* f_data)
{
    //优化目标函数
    //输入f_data：0~5初值，6~7飞行时间初值,8为t0时刻,9为星座编号,10为星座轨道面编号,11~13为目标三星编号
    //输入X：0~2初值速度改变值，3~4积分时间改变值

    double* para = (double*)f_data; //类型强制转换
    double s0[6];
    for (int i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            s0[i] = para[i];
        }
        else
        {
            s0[i] = para[i] + (2.0 * X[i - 3] - 1.0) * 300.0; //500m/s伸缩
        }
    }
    double tfa = para[6] + (2.0 * X[3] - 1.0) * 30.0; //30s伸缩
    double tfb = para[7] + (2.0 * X[4] - 1.0) * 30.0;
    double t0 = para[8];
    int flag1;
    double rvf1[6], rvf2[6];
    flag1 = orbitPropagateJ2(t0, t0 + tfa, s0, rvf1, 1.0, 1e-13, 1e-13); //第一段，飞向二号星
    int n1 = para[9];
    int n2 = para[10];
    int k2 = para[12];

    double rv2[6];
    get_constellation_RV(n1, n2, k2, t0 + tfa, rv2);

    double r_2[3], rf_1[3], dv_mane[3];
    for (int i = 0; i < 3; ++i)
    {
        r_2[i] = rv2[i];
        rf_1[i] = rvf1[i];
        dv_mane[i] = s0[i + 3] - para[i + 3];
    }
    double dv = sqrt(dv_mane[0] * dv_mane[0] + dv_mane[1] * dv_mane[1] + dv_mane[2] * dv_mane[2]);
    double cosphi = cosphi_obs(t0, r_2, rf_1);
    double J = -cosphi * 10.0 + dv;

    return J;
}

bool InPlaneTransfer(int n1, int n2, double t0, double* rv0, double* tf, double* rvf,
                     InPlaneTransfer::IPTResult& myResult, int dir)
{
    // 输入参数
    // n1 星座编号
    // n2 星座轨道面编号
    // t0 抵达轨道面时刻
    // lambda 抵达轨道面时的真经度

    bool success = true;

    std::array<double, 6> plane;
    plane = getPlane(n1, n2, t0);

    double lambda = DPI / 2.0;
    int N = static_cast<int>(plane[3]); //星座平面内卫星数量
    std::vector<double> dlambda_t(N);
    double coe[6];
    for (int i = 1; i <= N; ++i)
    {
        get_constellation_coe(n1, n2, i, t0, coe);
        dlambda_t[i - 1] = cos(fmod(coe[4] + coe[5], D2PI) - lambda); //近地点幅角和真近点角之和
    }
    auto max_element_iter = std::max_element(dlambda_t.begin(), dlambda_t.end()); // 查找最小元素的迭代器
    int max_index = std::distance(dlambda_t.begin(), max_element_iter); // 计算最小值的索引
    //double min_value = *min_element_iter; // 解引用以获取最小值

    //距离lambda最近的卫星编号
    int n3 = max_index + 1;
    double sma = plane[0]; //星座内卫星轨道半长轴
    if (dir == -1)
    {
        double sma_test = sma * pow(1.0 - 1.0 / N, 2.0 / 3.0);
        double h_min = sma_test * 2.0 - sma - 50.0e3 - 6378e3;
        if (h_min / 1e3 <= 200)
        {
            dir = 1;
        }
    }
    int k1, k2, k3;
    double sma_ro;
    if (dir == -1)
    {
        sma_ro = sma * pow(1.0 - 1.0 / N, 2.0 / 3.0); //共振轨道半长轴,(1-1/N)共振
        //初始化
        k1 = n3;
        k2 = k1 + 1;
        if (k2 > N)
        {
            k2 += -N;
        }
        k3 = k2 + 1;
        if (k3 > N)
        {
            k3 += -N;
        }
    }
    else
    {
        sma_ro = sma * pow(1.0 + 1.0 / N, 2.0 / 3.0); //共振轨道半长轴,(1+1/N)共振
        //初始化
        k1 = n3;
        k2 = k1 - 1;
        if (k2 <= 0)
        {
            k2 += N;
        }
        k3 = k2 - 1;
        if (k3 <= 0)
        {
            k3 += N;
        }
    }
    double T_ro = D2PI * sqrt(sma_ro * sma_ro * sma_ro / MU_EARTH); //共振轨道周期

    double rv1[6];
    //double tf_0, tf_1;//预估转移时间
    get_constellation_RV(n1, n2, k1, t0, rv1);
    double s0[6];
    double kr = (1.0 + 1e3 / sqrt(rv1[0] * rv1[0] + rv1[1] * rv1[1] + rv1[2] * rv1[2])); //巡检卫星初始位置伸缩
    double kv = sqrt(MU_EARTH * (2.0 / (kr * sqrt(rv1[0] * rv1[0] + rv1[1] * rv1[1] + rv1[2] * rv1[2])) - 1.0 / sma_ro))
        / sqrt(rv1[3] * rv1[3] + rv1[4] * rv1[4] + rv1[5] * rv1[5]); //巡检卫星初始速度伸缩
    for (int i = 0; i < 6; ++i)
    {
        if (i < 3)
        {
            s0[i] = kr * rv1[i];
        }
        else
        {
            s0[i] = kv * rv1[i];
        }
    }

    //从最近的卫星开始遍历
    double f_opt = 1e8;
    double f_data[14] = {0};
    for (size_t i = 0; i < 6; i++)
        f_data[i] = s0[i];
    f_data[6] = T_ro;
    f_data[7] = T_ro;
    f_data[8] = t0;
    f_data[9] = n1;
    f_data[10] = n2;
    f_data[11] = k1;
    f_data[12] = k2;
    f_data[13] = k3;
    std::vector<double> x_opt = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; //优化参数设置

    nlopt_main_S2S_1(Function_InPlane_S2S_1, Constraint_InPlane_S2S_1, f_data, x_opt, f_opt, 8, 0, 1e4);

    std::array<double, 6> rv_b;
    std::array<double, 6> rv_a;
    double s1[6];
    for (size_t i = 0; i < 6; i++)
    {
        if (i < 3)
        {
            s1[i] = s0[i] + (2.0 * x_opt[i] - 1.0) * 50e3;
        }
        else
        {
            s1[i] = s0[i] + (2.0 * x_opt[i] - 1.0) * 300.0;
        }
        rv0[i] = s1[i];
        rv_a[i] = s0[i];
        rv_b[i] = s1[i];
    }
    double tfa = T_ro + (2.0 * x_opt[6] - 1.0) * 20.0;
    double tfb = T_ro + (2.0 * x_opt[7] - 1.0) * 20.0;
    int flag1;
    std::array<double, 6> rvf1;
    flag1 = orbitPropagateJ2(t0, t0 + tfa, s1, rvf1.data(), 1.0, 1e-13, 1e-13); //第一段，飞向二号星
    std::array<int, 3> tar0 = {0, 0, 0};
    std::array<int, 3> tar1 = {n1, n2, k1};
    std::array<int, 3> tar2 = {n1, n2, k2};
    std::array<double, 3> dv0 = {0, 0, 0};
    std::array<double, 3> dv1;

    if (t0 + tfa > 90.0 * 86400.0)
    {
        *tf = t0;
        for (size_t i = 0; i < 6; i++)
        {
            rvf[i] = s1[i];
        }
        success = false;
        return success;
    }
    //myResult.timeList.push_back(t0);
    //myResult.eventList.push_back(0);
    //myResult.rvList.push_back(rv_b);
    //myResult.dvList.push_back(dv0);
    //myResult.targetList.push_back(tar0);

    myResult.timeList.push_back(t0);
    myResult.eventList.push_back(2);
    myResult.rvList.push_back(rv_b);
    myResult.dvList.push_back(dv0);
    myResult.targetList.push_back(tar1);

    myResult.timeList.push_back(t0 + tfa);
    myResult.eventList.push_back(2);
    myResult.rvList.push_back(rvf1);
    myResult.dvList.push_back(dv0);
    myResult.targetList.push_back(tar2);
    t0 = t0 + tfa;
    *tf = t0;
    //自第二星开始遍历，直至完成
    for (int list = 1; list < N - 1; list++)
    {
        if (dir == -1)
        {
            k1 = k2;
            k2 = k1 + 1;
            if (k2 > N)
            {
                k2 += -N;
            }
            k3 = k2 + 1;
            if (k3 > N)
            {
                k3 += -N;
            }
        }
        else
        {
            k1 = k2;
            k2 = k1 - 1;
            if (k2 <= 0)
            {
                k2 += N;
            }
            k3 = k2 - 1;
            if (k3 <= 0)
            {
                k3 += N;
            }
        }

        for (int i = 0; i < 6; ++i)
        {
            s0[i] = rvf1[i]; //重置s0
        }
        f_opt = 1e8;
        for (size_t i = 0; i < 6; i++)
            f_data[i] = s0[i];
        f_data[6] = tfb;
        f_data[7] = T_ro;
        f_data[8] = t0;
        f_data[9] = n1;
        f_data[10] = n2;
        f_data[11] = k1;
        f_data[12] = k2;
        f_data[13] = k3;
        x_opt = {0.5, 0.5, 0.5, 0.5, 0.5}; //优化参数设置

        // nloptmain(Function_InPlane_S2S_2, f_data, x_opt, f_opt, 5, 10000);
        nlopt_main_S2S_2(Function_InPlane_S2S_2, Constraint_InPlane_S2S_2, f_data, x_opt, f_opt, 5, 0, 3e4);
        double result[6];
        std::vector<double> grad;
        Constraint_InPlane_S2S_2(1, result, 1, x_opt.data(), grad.data(), f_data);

        // std::cout << "con: " << n1 << " " << n2 << " " << k2 << std::endl;
        // for (int i = 0; i < 6; ++i)
        // {
        //     std::cout << std::setprecision(8) << result[i] << std::endl;
        // }

        for (size_t i = 0; i < 6; i++)
        {
            if (i < 3)
            {
                s1[i] = s0[i];
            }
            else
            {
                s1[i] = s0[i] + (2.0 * x_opt[i - 3] - 1.0) * 300.0;
                dv1[i - 3] = (2.0 * x_opt[i - 3] - 1.0) * 300.0;
            }
            rv_a[i] = s0[i];
        }
        tfa = tfb + (2.0 * x_opt[3] - 1.0) * 30.0;
        tfb = T_ro + (2.0 * x_opt[4] - 1.0) * 30.0;
        flag1 = orbitPropagateJ2(t0, t0 + tfa, s1, rvf1.data(), 1.0, 1e-13, 1e-13); //第一段，飞向二号星
        tar1 = {n1, n2, k2};

        if (t0 + tfa > 90.0 * 86400.0)
        {
            success = false;
            break;
        }

        double dv = sqrt(dv1[0] * dv1[0] + dv1[1] * dv1[1] + dv1[2] * dv1[2]);
        if (dv > 1e-6)
        {
            myResult.timeList.push_back(t0);
            myResult.eventList.push_back(1);
            myResult.rvList.push_back(rv_a);
            myResult.dvList.push_back(dv1);
            myResult.targetList.push_back(tar0);
        }
        myResult.timeList.push_back(t0 + tfa);
        myResult.eventList.push_back(2);
        myResult.rvList.push_back(rvf1);
        myResult.dvList.push_back(dv0);
        myResult.targetList.push_back(tar1);
        t0 = t0 + tfa;
        *tf = t0;
        for (size_t i = 0; i < 6; i++)
        {
            rvf[i] = rvf1[i];
        }
    }
    return success;
}
