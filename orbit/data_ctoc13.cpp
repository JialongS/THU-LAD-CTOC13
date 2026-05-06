#include "data_ctoc13.h"
#include "constant.h"
#include "OrbitFun.h"
#include "OrbitMath.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <iostream>
#include <sstream>


std::vector<std::vector<double>> GLOBAL_CONSTELLATION_DATA;
std::unordered_map<std::tuple<int, int, int>, int> GLOBAL_CONSTELLATION_MAPPING;
std::vector<std::array<double, 7>> GLOBAL_PLANE_DATA;
std::unordered_map<std::tuple<int, int>, int> GLOBAL_PLANE_MAPPING;
std::unordered_map<std::tuple<int, int>, int> GLOBAL_PLANE_Count;

bool insertPlane(int n1, int n2, double semi, double inc, double raan, double dOmega, int& index)
{
    auto sp_str = std::make_tuple(n1, n2);
    auto result = GLOBAL_PLANE_MAPPING.insert(std::make_pair(sp_str, index));
    GLOBAL_PLANE_Count[sp_str]++;
    if (result.second)
    {
        index++;
        std::array<double, 7> plane{};
        plane[0] = semi;
        plane[1] = inc;
        plane[2] = raan;
        plane[3] = dOmega;
        plane[5] = n1;
        plane[6] = n2;
        GLOBAL_PLANE_DATA.push_back(plane);
    }
    return result.second;
}

void load_data_ctoc13()
{
    std::ifstream fin;
    fin.open("../data/MegaConstellation_Data.txt");

    if (!fin.is_open())
    {
        std::cout << "Can't open data file" << std::endl;
    }

    int node_num = 30188;
    GLOBAL_CONSTELLATION_DATA.clear();
    GLOBAL_CONSTELLATION_DATA.resize(node_num);
    for (int i = 0; i < node_num; i++)
        GLOBAL_CONSTELLATION_DATA[i].resize(9);
    GLOBAL_PLANE_DATA.clear();
    GLOBAL_PLANE_MAPPING.clear();
    GLOBAL_PLANE_Count.clear();
    int plane_num = 0;

    int int1, int2, int3; // 用于存储前三列整数
    //0：轨道半长轴；1：轨道偏心率；2：轨道倾角；3：轨道升交点赤经；4：轨道近地点幅角；5：轨道平近点角
    for (int i = 0; i < node_num; i++)
    {
        fin >> int1 >> int2 >> int3;
        auto identifier = std::make_tuple(int1, int2, int3);
        fin >> GLOBAL_CONSTELLATION_DATA[i][0] >> GLOBAL_CONSTELLATION_DATA[i][1]
            >> GLOBAL_CONSTELLATION_DATA[i][2] >> GLOBAL_CONSTELLATION_DATA[i][3]
            >> GLOBAL_CONSTELLATION_DATA[i][4] >> GLOBAL_CONSTELLATION_DATA[i][5];
        GLOBAL_CONSTELLATION_DATA[i][0] *= 1e3; //转换成m
        GLOBAL_CONSTELLATION_MAPPING[identifier] = i;

        double a = GLOBAL_CONSTELLATION_DATA[i][0];
        double e = GLOBAL_CONSTELLATION_DATA[i][1];
        double inc = GLOBAL_CONSTELLATION_DATA[i][2];
        double OMEGA = GLOBAL_CONSTELLATION_DATA[i][3];
        double omega = GLOBAL_CONSTELLATION_DATA[i][4];
        double f0 = GLOBAL_CONSTELLATION_DATA[i][5];

        int flag;
        double E0 = f2E(flag, f0, e);
        double M0 = E2M(flag, E0, e);
        GLOBAL_CONSTELLATION_DATA[i][5] = M0;

        double n = sqrt(MU_EARTH / (a * a * a));
        double c1 = pow((MU_EARTH / a), 3.5);
        double c2 = 1.0 / (1.0 - e * e);
        double c3 = J2_EARTH * RADIUS_EARTH * RADIUS_EARTH;
        double c4 = MU_EARTH * MU_EARTH * MU_EARTH;

        double dOmega = -c1 * c2 * c2 * 1.5 * c3 / c4 * cos(inc);
        double dw = -c1 * c2 * c2 * 1.5 * c3 / c4 * (2.5 * sin(inc) * sin(inc) - 2.0);
        double dM = n * (1.0 + 0.75 * c3 / a / a / pow(1.0 - e * e, 1.5) * (2.0 - 3.0 * sin(inc) * sin(inc)));
        GLOBAL_CONSTELLATION_DATA[i][6] = dOmega;
        GLOBAL_CONSTELLATION_DATA[i][7] = dw;
        GLOBAL_CONSTELLATION_DATA[i][8] = dM;
        insertPlane(int1, int2, a, inc, OMEGA, dOmega, plane_num);
    }
    fin.close();
    for (const auto& pair : GLOBAL_PLANE_Count)
    {
        GLOBAL_PLANE_DATA[GLOBAL_PLANE_MAPPING[pair.first]][4] = pair.second;
    }
}

void get_constellation_RV(int n1, int n2, int n3, double t, double* rv)
{
    //计算t时刻卫星id的coe
    auto identifier = std::make_tuple(n1, n2, n3);
    int id = GLOBAL_CONSTELLATION_MAPPING[identifier];
    double coe[6];
    int flag;
    for (int i = 0; i < 3; i++)
        coe[i] = GLOBAL_CONSTELLATION_DATA[id][i];

    double dt = t;
    double no = GLOBAL_CONSTELLATION_DATA[id][3] + GLOBAL_CONSTELLATION_DATA[id][6] * dt;
    double nw = GLOBAL_CONSTELLATION_DATA[id][4] + GLOBAL_CONSTELLATION_DATA[id][7] * dt;
    coe[3] = fmod(no, D2PI);
    coe[4] = fmod(nw, D2PI);

    double nm = GLOBAL_CONSTELLATION_DATA[id][5] + GLOBAL_CONSTELLATION_DATA[id][8] * dt;
    nm = fmod(nm, D2PI);


    double E = M2E(flag, nm, GLOBAL_CONSTELLATION_DATA[id][1]);
    coe[5] = E2f(flag, E, GLOBAL_CONSTELLATION_DATA[id][1]);

    //由coe计算rv
    coe2rv(flag, rv, coe, MU_EARTH);
}

void get_constellation_coe(int n1, int n2, int n3, double t, double* coe) // watch out 8.27
{
    //计算t时刻卫星id的coe
    auto identifier = std::make_tuple(n1, n2, n3);
    int id = GLOBAL_CONSTELLATION_MAPPING[identifier];
    //double coe[6];
    int flag;
    for (int i = 0; i < 3; i++)
        coe[i] = GLOBAL_CONSTELLATION_DATA[id][i];

    double dt = t;
    double no = GLOBAL_CONSTELLATION_DATA[id][3] + GLOBAL_CONSTELLATION_DATA[id][6] * dt;
    double nw = GLOBAL_CONSTELLATION_DATA[id][4] + GLOBAL_CONSTELLATION_DATA[id][7] * dt;
    coe[3] = fmod(no, D2PI);
    coe[4] = fmod(nw, D2PI);

    double nm = GLOBAL_CONSTELLATION_DATA[id][5] + GLOBAL_CONSTELLATION_DATA[id][8] * dt;
    nm = fmod(nm, D2PI);
    double E = M2E(flag, nm, GLOBAL_CONSTELLATION_DATA[id][1]);
    coe[5] = E2f(flag, E, GLOBAL_CONSTELLATION_DATA[id][1]);
}

/**
 *
 * @param n1 星座编号
 * @param n2 轨道面编号
 * @param t 时间
 * @return 0: 半长轴, 1: 倾角, 2: 升交点经度, 3: 卫星数量
 */
std::array<double, 6> getPlane(const int n1, const int n2, const double t)
{
    std::array<double, 6> plane{};
    auto sp_str = std::make_tuple(n1, n2);
    int index = GLOBAL_PLANE_MAPPING[sp_str];
    plane[0] = GLOBAL_PLANE_DATA[index][0];
    plane[1] = GLOBAL_PLANE_DATA[index][1];
    plane[2] = GLOBAL_PLANE_DATA[index][2] + GLOBAL_PLANE_DATA[index][3] * t;
    plane[3] = GLOBAL_PLANE_DATA[index][4];
    plane[4] = n1;
    plane[5] = n2;
    return plane;
}

void get_sun_rv(double t, double* rv)
{
    double a = 149703636.934590 * 1e3;
    double e = 0.017350;
    double inc = 23.436 / 180.0 * DPI;
    double Omega = 359.996 / 180.0 * DPI;
    double w = 281.913 / 180.0 * DPI;
    double f0 = 358.275 / 180.0 * DPI;
    double coe0[6] = {a, e, inc, Omega, w, f0};
    double coe1[6];
    int flag;
    coe02coef(flag, coe1, coe0, t, MU_SUN);
    coe2rv(flag, rv, coe1, MU_SUN);
}

double cosphi_obs(double t, const double* rv_tar, const double* rv_sat)
{
    double rv_sun[6];
    get_sun_rv(t, rv_sun);
    double r_obs[3] = {rv_sat[0] - rv_tar[0], rv_sat[1] - rv_tar[1], rv_sat[2] - rv_tar[2]};
    double r_sun[3] = {rv_sun[0] - rv_tar[0], rv_sun[1] - rv_tar[1], rv_sun[2] - rv_tar[2]};
    double cosphi = (r_obs[0] * r_sun[0] + r_obs[1] * r_sun[1] + r_obs[2] * r_sun[2]) / sqrt(
        r_obs[0] * r_obs[0] + r_obs[1] * r_obs[1] + r_obs[2] * r_obs[2]) / sqrt(
        r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]);
    return cosphi;
}
