#ifndef DATA_GTOC9_H
#define DATA_GTOC9_H

/**
 * data_gtoc9.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-8-9
 */

#include <vector>


/**
 碎片数据
 0: 轨道半长轴 a
 1：轨道偏心率 e
 2：轨道倾角 inc
 3：升交点赤经 Omega
 4：近地点幅角 omega
 5：平近点角 M
 6：轨道平均角速度 n
 7：升交点赤经变化率 dOmega
 8：近地点幅角变化率 dw
 9：轨道历元时刻
 */
extern std::vector<std::vector<double>> GLOBAL_DEBRIS_DATA;

void loadDebris(); // 载入碎片数据

void getDebrisRV(int id, double t, double* rv);

#endif //DATA_GTOC9_H
