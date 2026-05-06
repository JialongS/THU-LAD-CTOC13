#ifndef _DE_H_
#define _DE_H_

#include<vector>

void DE(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
        std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000, int ItOut = 50,
        double CR = -1.0);

// void DE_std(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data), void* f_data,
//             std::vector<double>& xbest, double& fbest, int num_variable, int Np = 0, int ItMax = 1000, int ItOut = 50,
//             double CR = -1.0);

#endif
