#ifndef CONSTANT_H
#define CONSTANT_H

/**
 * constant.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-7-28
 */
#include <numbers>

//PI
constexpr long double DPI = 3.1415926535897932384626433832795;
constexpr long double D2PI = 6.283185307179586476925286766559;
constexpr long double D2R = 0.017453292519943295769236907684886;
constexpr long double R2D = 57.295779513082320876798154814105;
//Time
constexpr double JD2S = 86400.0;
//Constant
constexpr double MU_EARTH = 398600.436233e9; //m^3/s^2
constexpr double MU_SUN = 1.32712428e20;
constexpr double RADIUS_EARTH = 6378.1370e3; //m
constexpr double EPSILON = 1.0e-14;
//J2
constexpr double J2_EARTH = 1.081874e-3; //CTOC13 data
//Pro

constexpr double r_min = 6578.1370e3; //m

#endif //CONSTANT_H
