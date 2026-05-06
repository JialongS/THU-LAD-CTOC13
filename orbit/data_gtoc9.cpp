#include "data_gtoc9.h"
#include "constant.h"
#include "OrbitFun.h"

#include <cmath>
#include <fstream>

std::vector<std::vector<double>> GLOBAL_DEBRIS_DATA;

void loadDebris()
{
    std::ifstream fin;
    fin.open("../data/MegaConstellation_Data.txt");

    int node_num = 123;
    GLOBAL_DEBRIS_DATA.clear();
    GLOBAL_DEBRIS_DATA.resize(node_num);
    for (int i = 0; i < node_num; i++)
        GLOBAL_DEBRIS_DATA[i].resize(10);

    double temp;
    //0：轨道半长轴；1：轨道偏心率；2：轨道倾角；3：轨道升交点赤经；4：轨道近地点幅角；5：轨道平近点角
    for (int i = 0; i < node_num; i++)
    {
        fin >> temp >> GLOBAL_DEBRIS_DATA[i][9] //输入数据
            >> GLOBAL_DEBRIS_DATA[i][0] >> GLOBAL_DEBRIS_DATA[i][1]
            >> GLOBAL_DEBRIS_DATA[i][2] >> GLOBAL_DEBRIS_DATA[i][3]
            >> GLOBAL_DEBRIS_DATA[i][4] >> GLOBAL_DEBRIS_DATA[i][5];
        double a = GLOBAL_DEBRIS_DATA[i][0];
        double e = GLOBAL_DEBRIS_DATA[i][1];
        double inc = GLOBAL_DEBRIS_DATA[i][2];
        double OMEGA = GLOBAL_DEBRIS_DATA[i][3];
        double omega = GLOBAL_DEBRIS_DATA[i][4];
        double M0 = GLOBAL_DEBRIS_DATA[i][5];
        double p = a * (1 - e * e);
        double c2 = (RADIUS_EARTH / p) * (RADIUS_EARTH / p);
        double ci = cos(inc);
        double n = sqrt(MU_EARTH / (a * a * a));
        double dOmega = -1.5 * J2_EARTH * c2 * n * ci;
        double dw = 0.75 * J2_EARTH * c2 * n * (5 * ci * ci - 1);
        GLOBAL_DEBRIS_DATA[i][6] = n;
        GLOBAL_DEBRIS_DATA[i][7] = dOmega;
        GLOBAL_DEBRIS_DATA[i][8] = dw;
        GLOBAL_DEBRIS_DATA[i][9] *= 86400.0;
    }
}

void getDebrisRV(int id, double t, double* rv)
{
    //计算t时刻碎片id的coe
    double coe[6];
    int flag;
    for (int i = 0; i < 3; i++)
        coe[i] = GLOBAL_DEBRIS_DATA[id][i];

    double dt = t - GLOBAL_DEBRIS_DATA[id][9];
    double no = GLOBAL_DEBRIS_DATA[id][3] + GLOBAL_DEBRIS_DATA[id][7] * dt;
    double nw = GLOBAL_DEBRIS_DATA[id][4] + GLOBAL_DEBRIS_DATA[id][8] * dt;
    coe[3] = fmod(no, D2PI);
    coe[4] = fmod(nw, D2PI);

    double nm = GLOBAL_DEBRIS_DATA[id][5] + GLOBAL_DEBRIS_DATA[id][6] * dt;
    nm = fmod(nm, D2PI);
    double E = M2E(flag, nm, GLOBAL_DEBRIS_DATA[id][1]);
    coe[5] = E2f(flag, E, GLOBAL_DEBRIS_DATA[id][1]);

    //由coe计算rv
    coe2rv(flag, rv, coe, MU_EARTH);
}
