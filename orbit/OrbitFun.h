#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_

#include "constant.h"

/**************************************************************************************************************************************/
/****************************************************经典轨道根数三种角度关系**********************************************************/
/**************************************************************************************************************************************/
double E2f(int& flag, double E, double e);
// E2M 根据偏近点角和偏心率求平近点角
double E2M(int& flag, double E, double e);
// f2E 根据真近点角和偏心率求偏近点角
double f2E(int& flag, double f, double e);
// M2E 根据平近点角和偏心率求偏近点角
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = EPSILON);
// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu = MU_EARTH, int MaxIter = 100, double epsilon = EPSILON);
// f0ft2dt 根据初始真近点角和最终真近点角求演化时间
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = MU_EARTH);
/**************************************************************************************************************************************/
/*******************************************************春分点轨道根数角度关系*********************************************************/
/**************************************************************************************************************************************/
double L0dt2Lt(int& flag, double L0, double dt, const double* ee, double mu = MU_EARTH);
double L0Lt2dt(int& flag, double L0, double Lt, const double* ee, double mu = MU_EARTH);
/**************************************************************************************************************************************/
/*********************************************经典轨道根数、直角坐标、改进春分点轨道根数的转换*****************************************/
/**************************************************************************************************************************************/
// coe2rv 根据经典轨道根数求地心惯性直角坐标系下的位置和速度分量
void coe2rv(int& flag, double* rv, const double* coe, double mu = MU_EARTH);
// rv2coe 根据地心惯性直角坐标系下的位置和速度分量求经典轨道根数
void rv2coe(int& flag, double* coe, const double* RV, double mu = MU_EARTH);
// coe2ee 根据经典轨道根数求改进春分点轨道根数，对轨道倾角180度奇异
void coe2ee(int& flag, double* ee, const double* coe, double mu = MU_EARTH);
// ee2coe 根据改进春分点根数求经典轨道根数
void ee2coe(int& flag, double* coe, const double* ee, double mu = MU_EARTH);
// ee2rv 根据改进春分点轨道根数求地心惯性直角坐标系下的位置和速度分量
void ee2rv(int& flag, double* rv, const double* ee, double mu = MU_EARTH);
// rv2ee 根据地心惯性直角坐标系下的位置和速度分量求改进春分点轨道根数，对轨道倾角180度时奇异
void rv2ee(int& flag, double* ee, const double* RV, double mu = MU_EARTH);
/**************************************************************************************************************************************/
/************************************************已知初值和时间，求末端状态************************************************************/
/**************************************************************************************************************************************/
//根据初始时刻状态coe0求末端时刻dt的状态coe1，按二体推进。若计算成功,flag返回1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu = MU_EARTH);
//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu = MU_EARTH);
//根据初始时刻状态ee0求末端时刻dt的状态ee1，按二体推进。若计算成功,flag返回1
void ee02eef(int& flag, double* ee1, const double* ee0, double dt, double mu = MU_EARTH);
/**************************************************************************************************************************************/
/*********************************************************考虑j2摄动的轨道根数转换*****************************************************/
/**************************************************************************************************************************************/
void coe2soe(const double* coe, double* soe);
void soe2coe(const double* soe, double* coe);
//输入为平均经典轨道根数，输出为瞬时轨道根数
void M2O(const double* soe_m, double* soe_o, double j2 = J2_EARTH);
void O2M(const double* soe_o, double* soe_m, double j2 = J2_EARTH);
int myfun(int n, const double* x, double* fvec, int iflag, const double* para);
/**************************************************************************************************************************************/
/************************************************************考虑j2摄动的轨道递推******************************************************/
/**************************************************************************************************************************************/
void j2mcoe02mcoef(const double* me0, const double dt, double* mef);
void j2ocoe02ocoef(double* ocoef, const double* ocoe0, const double dt);
void j2rv02rvf(const double* rv0, const double dt, double* rvf);

void j2mcoe02mcoef_eJ2(const double* me0, const double dt, const double eJ2, double* mef);
void j2ocoe02ocoef_eJ2(double* ocoef, const double* ocoe0, const double dt, const double eJ2);
void j2rv02rvf_eJ2(const double* rv0, const double dt, const double eJ2, double* rvf);
#endif
