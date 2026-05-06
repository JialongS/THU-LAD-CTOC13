/****************************************************************************
% * Copyright (C++), 2020-2031 清华大学航天航空学院动力学与控制实验室
% * 作者: 张楠
% * 文件名: nlopt_main.h
% * 内容简述：引入求解连续变量优化问题的nlopt算法，源文件
% * 文件历史：
% * 版本号     日期         作者       说明
% * 01a       2021-02-24    张楠     创建该文件
% */

#include"nlopt_main.h"

#include <iostream>

using namespace nlopt;

//D变量个数
void nloptmain(double (*ObjFun)(const std::vector<double>& x, std::vector<double>& grad, void* f_data), void* f_data,
               std::vector<double>& xbest, double& fbest, int num_variable, int ItMax)
{
	nlopt::opt opter(nlopt::LN_SBPLX, num_variable); //局部优化
	//nlopt::opt opter(nlopt::GN_CRS2_LM, num_variable);			//全局优化
	opter.set_min_objective(ObjFun, f_data); //指标
	//opter.set_population(2 * num_variable * num_variable);

	std::vector<double> lb(num_variable), rb(num_variable), dx(num_variable);
	//xbest.resize(num_variable);
	fbest = 1.0e20;
	for (int i = 0; i < num_variable; i++)
	{
		//xbest[i] = rand() / RAND_MAX;								//赋初值
		lb[i] = 0.0; //下界
		rb[i] = 1.0; //上界
		dx[i] = 1.0e-4; //初始步长
	}
	opter.set_initial_step(dx); //设置初始步长
	opter.set_lower_bounds(lb); //设置下届
	opter.set_upper_bounds(rb); //设置上届
	double tol = 1e-5;
	opter.set_xtol_rel(tol);
	//opter.set_force_stop(tol);


	opter.set_maxeval(ItMax); //优化该次数后停止

	nlopt::result res = opter.optimize(xbest, fbest); //进行优化
}

void nlopt_main_S2S_1(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
                      void (*inequality_constraint)(unsigned m, double* result, unsigned n, const double* X,
                                                    double* grad, void* f_data),
                      void* f_data, std::vector<double>& X, double& f, int num_variable, int Np, int maxstep)
{
	//nlopt::opt opter(nlopt::LD_MMA, num_variable);
	nlopt::opt opter(nlopt::LN_COBYLA, num_variable); //定义一个优化器，使用LN_BOBYQA算法不处理非线性约束，使用二次优化，速度更快
	//使用LN_COBYLA算法(不需要梯度且能处理非线性约束)，
	// GN_DIRECT_L
	//测试效果来看 LN_SBPLX、 GN_CRS2_LM  比较好
	// nlopt::opt opter(nlopt::LD_AUGLAG, num_variable);
	// nlopt::opt local_opt(nlopt::LN_NELDERMEAD, num_variable);
	// opter.set_local_optimizer(local_opt);
	int ineq_num = 9; //不等式约束几个维度
	double* para = (double*)f_data;

	std::vector<double> tol_ineq(ineq_num);
	for (int i = 0; i < tol_ineq.size(); i++) tol_ineq[i] = 1.0e-8;
	opter.add_inequality_mconstraint(inequality_constraint, f_data, tol_ineq);


	//double tol = 1e-12;
	opter.set_xtol_abs(1e-10);
	opter.set_ftol_abs(1e-6);
	opter.set_force_stop(1e-12);
	opter.set_maxeval(maxstep); //优化该次数后停止

	std::vector<double> lb(X.size()), ub(X.size()), dx(X.size()); //下界 //上界
	for (int i = 0; i < X.size(); i++)
	{
		lb[i] = 0.0;
		ub[i] = 1.0;
		dx[i] = 1e-3;
	}

	opter.set_lower_bounds(lb);
	opter.set_upper_bounds(ub);
	opter.set_initial_step(dx); //设置初始步长,也可以不设置
	opter.set_min_objective(ObjFun, f_data); //指标

	nlopt::result res = opter.optimize(X, f);
	//std::cout << " NLOPT返回值：" << res << std::endl;
	//std::cout << " NLOPT返回fbest：" << std::setprecision(14) << f << std::endl;
	//int numers = opter.get_numevals();
	//std::cout << " NLOPT计算次数：" << numers << std::endl;
}

void nlopt_main_S2S_2(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
                      void (*inequality_constraint)(unsigned m, double* result, unsigned n, const double* X,
                                                    double* grad, void* f_data),
                      void* f_data, std::vector<double>& X, double& f, int num_variable, int Np, int maxstep)
{
	//nlopt::opt opter(nlopt::LN_SBPLX, num_variable);
	//nlopt::opt opter(nlopt::LN_COBYLA, num_variable); //定义一个优化器，使用LN_BOBYQA算法不处理非线性约束，使用二次优化，速度更快

	/**
	*nlopt::opt opter(nlopt::LD_AUGLAG, num_variable);
	nlopt::opt local_opt(nlopt::LN_NELDERMEAD, num_variable);
	是可用的，对于带不等式约束的问题
	 */
	nlopt::opt opter(nlopt::LD_AUGLAG, num_variable);
	nlopt::opt local_opt(nlopt::LN_NELDERMEAD, num_variable);
	opter.set_local_optimizer(local_opt);
	//使用LN_COBYLA算法(不需要梯度且能处理非线性约束)，
	// GN_DIRECT_L
	//测试效果来看 LN_SBPLX、 GN_CRS2_LM  比较好
	int ineq_num = 6; //不等式约束几个维度
	double* para = (double*)f_data;

	std::vector<double> tol_ineq(ineq_num);
	for (int i = 0; i < tol_ineq.size(); i++) tol_ineq[i] = 1.0e-10;
	opter.add_inequality_mconstraint(inequality_constraint, f_data, tol_ineq);


	//double tol = 1e-12;
	opter.set_xtol_abs(1e-8);
	opter.set_ftol_abs(1e-4);
	opter.set_force_stop(1e-14);
	opter.set_maxeval(maxstep); //优化该次数后停止

	std::vector<double> lb(X.size()), ub(X.size()), dx(X.size()); //下界 //上界
	for (int i = 0; i < X.size(); i++)
	{
		lb[i] = 0.0;
		ub[i] = 1.0;
		dx[i] = 1e-2;
	}

	opter.set_lower_bounds(lb);
	opter.set_upper_bounds(ub);
	opter.set_initial_step(dx); //设置初始步长,也可以不设置
	opter.set_min_objective(ObjFun, f_data); //指标

	nlopt::result res = opter.optimize(X, f);
	//std::cout << " NLOPT return: " << res << std::endl;
	//std::cout << " NLOPT返回fbest：" << std::setprecision(14) << f << std::endl;
	//int numers = opter.get_numevals();
	//std::cout << " NLOPT计算次数：" << numers << std::endl;
}
