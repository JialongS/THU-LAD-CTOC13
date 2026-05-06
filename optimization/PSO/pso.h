#ifndef _PSO_H_
#define _PSO_H_
#include<vector>


//num_variable变量个数，Np粒子个数
void PSO(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
         std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000,
         int ItOut = 50, double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5,
         double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.8);

void PSO_parallel2(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),
                   void* f_data,
                   std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000,
                   int ItOut = 50,
                   double OmegaMin = 0.4, double OmegaMax = 0.9, double C1Min = 0.5, double C1Max = 2.5,
                   double C2Min = 0.5, double C2Max = 2.5, double Vmax = 0.8);

#endif
