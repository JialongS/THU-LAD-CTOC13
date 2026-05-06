#ifndef _ORBITBASE_H_
#define _ORBITBASE_H_
#include <array>

//求解J2lambert问题
void J2lambert(double* v1, double* v2, const double* rv1, const double* rv2, double tf,
               int& flag, double mu, int Maxiter = 200, double tol = 1.0e-11);

void J2_Lambert_MRPLP(double* v1, double* v2, const double* rv1, const double* rv2, double tof, int& flag, double MU,
                      int Maxiter = 40, double ftol = 1.0);

// void J2_Lambert_MRPLP_best(double* v1, double* v2, const double* rv1, const double* rv2, double tof,
//                            std::array<double, 3>& r2_true,
//                            int& flag, double MU, double rmin = 1e3, double rmax = 50e3, int Maxiter = 40,
//                            double ftol = 1.0);

//求解J2lambert问题时的子函数，用于打靶
//输入：n为打靶变量个数，x为t0时刻输入的dv，fvec为打靶误差，iflag=1表示该段程序正常运行，
//para[0:5]表示rv0，para[6:11]表示rvf，para[12]表示t0，para[13]表示tf，para[14]表示J2同伦参数
int SubFun_J2lambert(int n, const double* x, double* fvec, int iflag, const double* para);

//求解J2lambert问题（其中积分改用轨道根数递推）
void J2lambert_Approximation(double* v1, double* v2, const double* rv1, const double* rv2, double tf,
                             int& flag, double mu, int Maxiter = 200, double tol = 1.0e-11);

//求解J2lambert问题（其中积分改用轨道根数递推）
void J2lambert_Approximation_1(double* v1, double* v2, const double* rv1, const double* rv2, double tf,
                               int& flag, double mu, int Maxiter = 200, double tol = 1.0e-11);

//求解J2lambert问题（其中积分改用轨道根数递推）
void J2lambert_Approximation_400(double* v1, double* v2, const double* rv1, const double* rv2, double tf,
                                 int& flag, double mu, int Maxiter = 200, double tol = 1.0e-11);

//求解J2lambert问题时的子函数，用于打靶（其中积分改用轨道根数递推）
//输入：n为打靶变量个数，x为t0时刻输入的dv，fvec为打靶误差，iflag=1表示该段程序正常运行，
//para[0:5]表示rv0，para[6:11]表示rvf，para[12]表示t0，para[13]表示tf，para[14]表示J2同伦参数
int SubFun_J2lambert_Approximation(int n, const double* x, double* fvec, int iflag, const double* para);

//利用轨道根数递推进行积分
int Propagate_Approximation(const double* rv0, double t0, double t1, double* rv1, const double* para);

#endif
