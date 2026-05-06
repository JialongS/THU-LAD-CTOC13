#include"pso.h"
#include<cstdio>
#include<cstdlib>
#include<iostream>
//#include<stdlib.h>
//#include<time.h>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>

#include "rand_safe.h"
//#include<process.h>

//using namespace std;

void PSO(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
         std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut,
         double OmegaMin, double OmegaMax, double C1Min, double C1Max, double C2Min, double C2Max, double Vmax)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;

	int ct;
	double Omega, C1, C2, GBVAL, rand0;
	double* PBVAL = new double[Np];
	double* GBPOS = new double[num_variable];
	double** popvel = NULL;
	popvel = new double*[Np];
	for (int i = 0; i < Np; i++)
		popvel[i] = new double[num_variable];

	//double** pop=new double*[Np];
	//for(int i=0;i<Np;i++)
	//	pop[i]=new double[num_variable];
	std::vector<std::vector<double>> pop(Np);
	for (int i = 0; i < Np; i++)
		pop[i].resize(num_variable);
	double** PBPOS = new double*[Np];
	for (int i = 0; i < Np; i++)
		PBPOS[i] = new double[num_variable];

	for (int j = 0; j < Np; j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;
	//srand( (unsigned)time( NULL ) );
	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j < num_variable; j++)
		{
			rand0 = (double)rand() / RAND_MAX; //第一个丢掉
			pop[i][j] = (double)rand() / RAND_MAX;
			popvel[i][j] = 0.0;
		}
	}

	ct = 1;
	while (ct <= ItMax)
	{
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < Np; i++)
		{
			double val = ObjFun(pop[i], grad, f_data);
			if (val < PBVAL[i])
			{
				PBVAL[i] = val;
				for (int j = 0; j < num_variable; j++)
					PBPOS[i][j] = pop[i][j];
			}
#pragma omp critical
			if (val < GBVAL)
			{
				GBVAL = val;
				for (int j = 0; j < num_variable; j++)
					GBPOS[j] = pop[i][j];
			}
		}
		Omega = OmegaMax - (OmegaMax - OmegaMin) / ItMax * ct;
		C1 = -(C1Max - C1Min) * ct / ItMax + C1Max;
		C2 = (C2Max - C2Min) * ct / ItMax + C2Min;
		for (int i = 0; i < Np; i++)
		{
			for (int j = 0; j < num_variable; j++)
			{
				rand0 = (double)rand() / RAND_MAX;
				popvel[i][j] = Omega * popvel[i][j] + C1 * (PBPOS[i][j] - pop[i][j]) * (double)rand() / RAND_MAX
					+ C2 * (GBPOS[j] - pop[i][j]) * (double)rand() / RAND_MAX;
				if (popvel[i][j] > Vmax) popvel[i][j] = Vmax; //
				else if (popvel[i][j] < -Vmax) popvel[i][j] = -Vmax; //
				pop[i][j] += popvel[i][j];
				if ((pop[i][j] > 1.0) || (pop[i][j] < 0.0))
				{
					rand0 = (double)rand() / RAND_MAX;
					pop[i][j] = (double)rand() / RAND_MAX;
					popvel[i][j] = 0.0;
				}
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
	for (int i = 0; i < Np; i++)
		delete[] popvel[i];
	delete[] popvel;
	//for(int i=0;i<Np;i++)
	//	delete[] pop[i];
	//delete[] pop;
	for (int i = 0; i < Np; i++)
		delete[] PBPOS[i];
	delete[] PBPOS;
}


//PSO并行版，进一步优化了并行效率。并行不稳定版本，要求Np不小于线程数，ObjFun的时间复杂度远大于加减乘除效率
//输入：	最小化目标函数
//			1. ObjFun，目标函数。其输入量X为优化变量，grad为梯度（不必输入），f_data为常值参数
//			2. f_data为常值阐述
//			3. xbest为得到的最优优化变量，每个分量均在[0.0, 1.0]之间
//			4. fbest为得到的最优目标数值
//			5. num_variable为优化变量的个数
//			6. Np为PSO算法中种群的个数，默认输入零时，Np在算法中自动设置为10*num_variable
//			7. ItMax为PSO算法最大迭代次数，默认为1000
//			8. ItOut为PSO算法中每隔ItOut代输出一次最优结果，默认为50
void PSO_parallel2(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
                   void* f_data,
                   std::vector<double>& xbest, double& fbest, int num_variable, int Np, int ItMax, int ItOut,
                   double OmegaMin, double OmegaMax, double C1Min, double C1Max, double C2Min, double C2Max,
                   double Vmax)
{
	std::vector<double> grad;
	if (Np == 0)
		Np = 10 * num_variable;

	double GBVAL;
	double* PBVAL = new double[Np]();
	double* GBPOS = new double[num_variable]();
	double** popvel = NULL;
	popvel = new double*[Np];
	for (int i = 0; i < Np; i++)
		popvel[i] = new double[num_variable];

	std::vector<std::vector<double>> pop(Np);
	for (int i = 0; i < pop.size(); i++)
		pop[i].resize(num_variable);
	double** PBPOS = new double*[Np];
	for (int i = 0; i < Np; i++)
		PBPOS[i] = new double[num_variable]();

	for (int j = 0; j < Np; j++)
		PBVAL[j] = 1.0e10;
	GBVAL = 1.0e10;
	for (int i = 0; i < Np; i++)
	{
		for (int j = 0; j < num_variable; j++)
		{
			pop[i][j] = realRand(0.0, 1.0);
			popvel[i][j] = 0.0;
		}
	}

	std::vector<std::vector<double>> pop_large(Np * ItMax);
	for (int i = 0; i < pop_large.size(); i++)
		pop_large[i].resize(num_variable);
	double** popvel_large = NULL;
	popvel_large = new double*[Np * ItMax];
	for (int i = 0; i < Np * ItMax; i++)
		popvel_large[i] = new double[num_variable]();

#pragma omp parallel for schedule(dynamic)
	for (int iter = 0; iter < pop_large.size(); iter++)
	{
		int ct = iter / Np + 1;
		int iter_i = iter % Np;
		double Omega = OmegaMax - (OmegaMax - OmegaMin) / ItMax * ct;
		double C1 = -(C1Max - C1Min) * ct / ItMax + C1Max;
		double C2 = (C2Max - C2Min) * ct / ItMax + C2Min;
		for (int j = 0; j < num_variable; j++)
		{
			popvel_large[iter][j] = Omega * popvel[iter_i][j]
				+ C1 * (PBPOS[iter_i][j] - pop[iter_i][j]) * realRand(0.0, 1.0)
				+ C2 * (GBPOS[j] - pop[iter_i][j]) * realRand(0.0, 1.0);
			if (popvel_large[iter][j] > Vmax) popvel_large[iter][j] = Vmax; //
			else if (popvel_large[iter][j] < -Vmax) popvel_large[iter][j] = -Vmax; //
			pop_large[iter][j] = pop[iter_i][j] + popvel_large[iter][j];
			if ((pop_large[iter][j] > 1.0) || (pop_large[iter][j] < 0.0))
			{
				pop_large[iter][j] = realRand(0.0, 1.0);
				popvel_large[iter][j] = 0.0;
			}
		}
		for (int j = 0; j < num_variable; j++)
		{
			pop[iter_i][j] = pop_large[iter][j];
			popvel[iter_i][j] = popvel_large[iter][j];
		}

		double val = ObjFun(pop_large[iter], grad, f_data);
#pragma omp critical
		{
			if (val < GBVAL)
			{
				GBVAL = val;
				for (int j = 0; j < num_variable; j++)
					GBPOS[j] = pop_large[iter][j];
			}
			if (val < PBVAL[iter_i])
			{
				PBVAL[iter_i] = val;
				for (int j = 0; j < num_variable; j++)
					PBPOS[iter_i][j] = pop_large[iter][j];
			}
			//for (int j = 0;j < num_variable;j++)
			//{
			//	pop[iter_i][j] = pop_large[iter][j];
			//	popvel[iter_i][j] = popvel_large[iter][j];
			//}
			if (ct % ItOut == 0 && iter_i == Np - 1)
			{
				std::cout << "No. of iteration=" << ct << std::endl;
				for (int i = 0; i < num_variable; i++)
					std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
				std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;
			}
		}
	}

	for (int j = 0; j < num_variable; j++)
		xbest[j] = GBPOS[j];
	fbest = GBVAL;

	//for (int i = 0;i < num_variable;i++)
	//	std::cout << "xbest(" << i + 1 << ")=" << std::setprecision(15) << GBPOS[i] << std::endl;
	//std::cout << "fbest=" << std::setprecision(15) << GBVAL << std::endl << std::endl;

	delete[] PBVAL;
	delete[] GBPOS;
	for (int i = 0; i < Np; i++)
		delete[] popvel[i];
	delete[] popvel;
	for (int i = 0; i < Np * ItMax; i++)
		delete[] popvel_large[i];
	delete[] popvel_large;
	for (int i = 0; i < Np; i++)
		delete[] PBPOS[i];
	delete[] PBPOS;
}

