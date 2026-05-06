#ifndef _LAMBERT_H_
#define _LAMBERT_H_

// Inputs:
//           R1 = Position array at departure
//           R2 = Position array at arrival
//           TOF = Transfer time (scalar)
//           mu = gravitational parameter (scalar, units have to be consistent with r1,t units)
//           way = 0 for counterclockwise transfer, 1 for clockwise transfer, default 0
//           N = number of revolutions, default 0
//           branch = 0 if the left branch is selected in a problem where N is not 0 (multirevolution), 1 for the right one, default 0
//           Maxiter = the maximum iteration, default 60
//           tol=tolerance, default 1.0e-11
//
// Outputs:
//           v1vec = Velocity at departure        (consistent units)
//           v2vec = Velocity at arrival
//			 return flag (int): 1, Successfully; 2, poor precision because transfer time too small
//								no solution: -1, negative time; -2, the rectilinear case; -3, TOF<Tmin; -4, not to be converging
int lambert(double* v1vec, double* v2vec, double& a, double& e, const double* R1, const double* R2, double TOF,
            double mu,
            int way = 0, int N = 0, int branch = 0, int Maxiter = 60, double tol = 1.0e-11);

// Inputs:
//           R1 = Position array at departure
//           R2 = Position array at arrival
//           TOF = Transfer time (scalar)
//           mu = gravitational parameter (scalar, units have to be consistent with r1,t units)
//           way = 0 for counterclockwise transfer, 1 for clockwise transfer, default 0
//           N = number of revolutions, default 0
//           branch = 0 if the left branch is selected in a problem where N is not 0 (multirevolution), 1 for the right one, default 0
//           Maxiter = the maximum iteration, default 60
//           tol=tolerance, default 1.0e-11
//
// Outputs:
//			 PD = Partial (v1vec, v2vec) / Partial (R1, R2, t1, t2), 6 * 8 dimension array, where t1 and t2 the epochs of R1 and R2, respectively.
//           v1vec = Velocity at departure        (consistent units)
//           v2vec = Velocity at arrival
//			 return flag (int): 1, Successfully; 2, poor precision because transfer time too small
//								no solution: -1, negative time; -2, the rectilinear case; -3, TOF<Tmin; -4, not to be converging
int lambert(double* PD, double* v1vec, double* v2vec, double& a, double& e, const double* R1, const double* R2,
            double TOF, double mu,
            int way = 0, int N = 0, int branch = 0, int Maxiter = 60, double tol = 1.0e-11);

//**************************************External Function End************************************************//

//**************************************Internal Function Begin**********************************************//
// Note: Users do not have to understand Internal Functions below

//this routine returns the W functionand its derivatives with respect to k up to order 1 - 4
void lambert_getD4W(double& W0, double& W1, double& W2, double& W3, double& W4, double k, int N);

//only calculate W0, only used in getinitialK Function
void lambert_getW(double& W0, double k, int N);

//if k is big, this is a precision saving version evaluation of the TOF/S equation and its partials
//the normal approach works fine, but leads to poor precision due to the (tau+pW) part of the TOF eq,
//so we use a series solution for the whole TOF eq.
void lambert_getFuncAndKdersBigK(double& Func0, double& Func1, double& Func2, double& Func3, double k, int N,
                                 double tau, double tofbyS);

//the TOF function and its partials wrt k
void lambert_getFuncAndKders(double& Func0, double& Func1, double& Func2, double& Func3, double& W0, double& W1,
                             double& W2, double& W3, double& W4, double k, double p, int N, double tau, double tofbyS);

//ith derivative of the function wrt the value
//ith order (isolated to that order only) correction
void lambert_getCorrection(double& dval0, double& dval1, double& dval2, double Func0, double Func1, double Func2,
                           double Func3);

//calculate the exact value of k of the minimum times for i-rev solutions.
double lambert_getkbi(double k, int N, double tau);

//initial guess generation
int lambert_getinitialK(double& k, double& kmin, double& kmax, int direction, int N, double tau, double tofbyS,
                        int branch);

//**************************************Internal Function End************************************************//


/***********************************Izzo组的Lambert程序**********************************/
//Subfunction that evaluates the time of flight as a function of x
double x2tof(double x, double s, double c, int lw, int N);

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
             double tf, const double* unith, int& flag, double mu, int way = 0, int N = 0, int branch = 0,
             int Maxiter = 200, double tol = 1.0e-11);

//多圈lambert遍历求解
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
              const double* rv0, const double* rv1, double t, double GM);

#endif
