#ifndef DATA_CTOC13_H
#define DATA_CTOC13_H

/**
 * data_ctoc13.h
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
#include <string>
#include <unordered_map>
#include <array>
#include <functional>

/**
 星座卫星数据
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

namespace std
{
 template <>
 struct hash<std::tuple<int, int, int>>
 {
  size_t operator()(const std::tuple<int, int, int>& k) const noexcept
  {
   const auto h1 = std::hash<int>()(std::get<0>(k));
   const auto h2 = std::hash<int>()(std::get<1>(k));
   const auto h3 = std::hash<int>()(std::get<2>(k));
   return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
 };
}

namespace std
{
 template <>
 struct hash<std::tuple<int, int>>
 {
  size_t operator()(const std::tuple<int, int>& k) const noexcept
  {
   const auto h1 = std::hash<int>()(std::get<0>(k));
   const auto h2 = std::hash<int>()(std::get<1>(k));
   return h1 ^ (h2 << 1);
  }
 };
}

extern std::vector<std::vector<double>> GLOBAL_CONSTELLATION_DATA;
extern std::unordered_map<std::tuple<int, int, int>, int> GLOBAL_CONSTELLATION_MAPPING;

/**
 * 星座轨道平面数据
 * 0: 轨道半长轴 a
 * 1: 轨道倾角 i
 * 2: 升交点精度 RAAN
 * 3: 升交点赤经变化率 dOmega
 * 4: 包含的卫星数
 * 5: n1
 * 6: n2
 */
extern std::vector<std::array<double, 7>> GLOBAL_PLANE_DATA; // 存储了轨道平面的信息
extern std::unordered_map<std::tuple<int, int>, int> GLOBAL_PLANE_MAPPING;
extern std::unordered_map<std::tuple<int, int>, int> GLOBAL_PLANE_Count;

bool insertPlane(int n1, int n2, double semi, double inc, double raan, double dOmega, int& index);

void load_data_ctoc13(); // 载入碎片数据

void get_constellation_RV(int n1, int n2, int n3, double t, double* rv);

void get_constellation_coe(int n1, int n2, int n3, double t, double* coe);

std::array<double, 6> getPlane(int n1, int n2, double t);

void get_sun_rv(double t, double* rv);

double cosphi_obs(double t, const double* rv_tar, const double* rv_sat);
#endif //DATA_CTOC13_H
