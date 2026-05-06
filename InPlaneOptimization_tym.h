#include "constant.h"
#include "dynamics.h"
#include "inPlaneOptimization.h"
#include "J2Lambert.h"
#include "io.h"
#include "optimization.h"
#include <cmath>
#include <iostream>


double Function_InPlane_S2S_1(const std::vector<double>& X, std::vector<double>& grad, void* f_data);

void Constraint_InPlane_S2S_1(unsigned m, double* result, unsigned n, const double* X, double* grad, void* f_data);

double Function_InPlane_S2S_2(const std::vector<double>& X, std::vector<double>& grad, void* f_data);

void Constraint_InPlane_S2S_2(unsigned m, double* result, unsigned n, const double* X, double* grad, void* f_data);

bool InPlaneTransfer(int n1, int n2, double t0, double* rv0, double* tf, double* rvf,
                     InPlaneTransfer::IPTResult& myResult, int dir);

