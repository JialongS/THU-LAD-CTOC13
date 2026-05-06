/****************************************************************************
% * Copyright (C++), 2020-2031 清华大学航天航空学院动力学与控制实验室
% * 作者: 张楠
% * 文件名: nlopt_main.h
% * 内容简述：引入求解连续变量优化问题的nlopt算法，头文件
% * 文件历史：
% * 版本号     日期         作者       说明
% * 01a       2021-02-24    张楠     创建该文件
% */

#ifndef _NLOPT_MAIN_H_
#define _NLOPT_MAIN_H_

#include<vector>
#include"nlopt.hpp"

//D变量个数
void nloptmain(double (*ObjFun)(const std::vector<double>& x, std::vector<double>& grad, void* f_data), void* f_data,
               std::vector<double>& xbest, double& fbest, int num_variable, int ItMax = 10000);

void nlopt_main_S2S_1(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
                      void (*inequality_constraint)(unsigned m, double* result, unsigned n, const double* X,
                                                    double* grad, void* f_data),
                      void* f_data, std::vector<double>& X, double& f, int num_variable, int Np, int maxstep);

void nlopt_main_S2S_2(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
                      void (*inequality_constraint)(unsigned m, double* result, unsigned n, const double* X,
                                                    double* grad, void* f_data),
                      void* f_data, std::vector<double>& X, double& f, int num_variable, int Np, int maxstep);


#endif
