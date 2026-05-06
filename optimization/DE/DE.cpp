#include "DE.h"
#include "rand_safe.h"
#include<iomanip>
#include<iostream>
#include<omp.h>

#include <algorithm>
#include <execution>
#include <mutex>

void DE(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
        std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut, double CR)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;
	if (CR < 0.0)
		CR = (18.0 + num_variable) / 40.0;
	if (CR > 0.9)
		CR = 0.9;

	int ct, r0, r1, r2;
	double F, GBVAL, lambda;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[num_variable];

	std::vector<std::vector<double>> pop(Np);
	for (int i = 0; i < Np; i++)
		pop[i].resize(num_variable);
	std::vector<std::vector<double>> u(Np);
	for (int i = 0; i < Np; i++)
		u[i].resize(num_variable);

	for (int j = 0; j < Np; j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j < num_variable; j++)
		{
			pop[i][j] = realRand(0.0, 1.0);
		}
		PBVAL[i] = ObjFun(pop[i], grad, f_data);
#pragma omp critical
		if (PBVAL[i] < GBVAL)
		{
			GBVAL = PBVAL[i];
			for (int j = 0; j < num_variable; j++)
				GBPOS[j] = pop[i][j];
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
		lambda = exp(1 - ItMax / (ItMax + 1 - ct));
		F = 0.5 * pow(2, lambda);

#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < Np; i++)
		{
			r0 = intRand(0, Np - 1) % Np;
			r1 = intRand(0, Np - 1) % Np;
			r2 = intRand(0, Np - 1) % Np;
			for (int j = 0; j < num_variable; j++)
			{
				if (realRand(0.0, 1.0) <= CR)
					u[i][j] = pop[r0][j] + F * (pop[r1][j] - pop[r2][j]);
				else
					u[i][j] = pop[i][j];
				if ((u[i][j] > 1.0) || (u[i][j] < 0.0))
				{
					u[i][j] = realRand(0.0, 1.0);
				}
			}
			double val = ObjFun(u[i], grad, f_data);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0; j < num_variable; j++)
					pop[i][j] = u[i][j];
			}
		}

		for (int i = 0; i < Np; i++)
		{
			if (PBVAL[i] < GBVAL)
			{
				GBVAL = PBVAL[i];
				for (int j = 0; j < num_variable; j++)
					GBPOS[j] = pop[i][j];
			}
		}

		for (int j = 0; j < num_variable; j++)
			xbest[j] = GBPOS[j];
		fbest = GBVAL;
		if (ct % ItOut == 0)
		{
			std::cout << "No. of iteration=" << ct << std::endl;
			for (int i = 0; i < num_variable; i++)
				std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
			std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
		}
		ct = ct + 1;
	}

	delete[] PBVAL;
	delete[] GBPOS;
}

// void DE_std(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
//             std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut, double CR)
// {
// 	static std::mutex update_mutex;
//
// 	std::vector<double> grad;
// 	if (Np == 0)
// 		Np = 10 * num_variable;
// 	if (CR < 0.0)
// 		CR = (18.0 + num_variable) / 40.0;
// 	if (CR > 0.9)
// 		CR = 0.9;
//
// 	double F, GBVAL, lambda;
// 	double* PBVAL = new double[Np];
// 	double* GBPOS = new double[num_variable];
//
// 	std::vector<std::vector<double>> pop(Np);
// 	for (int i = 0; i < Np; i++)
// 		pop[i].resize(num_variable);
// 	std::vector<std::vector<double>> u(Np);
// 	for (int i = 0; i < Np; i++)
// 		u[i].resize(num_variable);
//
// 	for (int j = 0; j < Np; j++)
// 		PBVAL[j] = 1.0e10;
// 	GBVAL = 1.0e10;
//
// 	std::for_each(std::execution::par, pop.begin(), pop.end(), [&](std::vector<double>& individual)
// 	{
// 		for (int j = 0; j < num_variable; ++j)
// 		{
// 			individual[j] = realRand(0.0, 1.0);
// 		}
// 		const size_t index = &individual - &pop[0]; // 计算当前解的索引
// 		PBVAL[index] = ObjFun(individual, grad, f_data);
//
// 		// 同步更新全局最佳解
// 		{
// 			std::lock_guard<std::mutex> lock(update_mutex);
// 			if (PBVAL[index] < GBVAL)
// 			{
// 				GBVAL = PBVAL[index];
// 				for (int j = 0; j < num_variable; ++j)
// 					GBPOS[j] = individual[j];
// 			}
// 		}
// 	});
//
// 	int ct = 1;
// 	while (ct <= ItMax)
// 	{
// 		lambda = exp(1 - ItMax / (ItMax + 1 - ct));
// 		F = 0.5 * pow(2, lambda);
//
// 		std::for_each(std::execution::par, u.begin(), u.end(), [&](std::vector<double>& u_i)
// 		{
// 			size_t index = &u_i - &u[0]; // 计算当前解的索引
// 			int r0 = intRand(0, Np - 1) % Np;
// 			int r1 = intRand(0, Np - 1) % Np;
// 			int r2 = intRand(0, Np - 1) % Np;
// 			for (int j = 0; j < num_variable; j++)
// 			{
// 				if (realRand(0.0, 1.0) <= CR)
// 					u_i[j] = pop[r0][j] + F * (pop[r1][j] - pop[r2][j]);
// 				else
// 					u_i[j] = pop[index][j];
// 				if ((u_i[j] > 1.0) || (u_i[j] < 0.0))
// 				{
// 					u_i[j] = realRand(0.0, 1.0);
// 				}
// 			}
// 			double val = ObjFun(u_i, grad, f_data);
// 			{
// 				std::lock_guard<std::mutex> lock(update_mutex);
// 				if (val < PBVAL[index])
// 				{
// 					PBVAL[index] = val;
// 					for (int j = 0; j < num_variable; j++)
// 						pop[index][j] = u[index][j];
// 				}
// 			}
// 		});
//
// 		for (int i = 0; i < Np; i++)
// 		{
// 			if (PBVAL[i] < GBVAL)
// 			{
// 				GBVAL = PBVAL[i];
// 				for (int j = 0; j < num_variable; j++)
// 					GBPOS[j] = pop[i][j];
// 			}
// 		}
//
// 		for (int j = 0; j < num_variable; j++)
// 			xbest[j] = GBPOS[j];
// 		fbest = GBVAL;
// 		if (ct % ItOut == 0)
// 		{
// 			std::cout << "No. of iteration=" << ct << std::endl;
// 			for (int i = 0; i < num_variable; i++)
// 				std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
// 			std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
// 		}
// 		ct = ct + 1;
// 	}
//
// 	delete[] PBVAL;
// 	delete[] GBPOS;
// }
