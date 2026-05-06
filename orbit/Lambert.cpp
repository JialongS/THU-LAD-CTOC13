//*********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************
// Originally programmed with Fortran source code by Ryan P. Russell.	Date: Release 1.06, 10-07-2019 see http://russell.ae.utexas.edu/index_files/lambert.html
//
// Translated into C++ source code by Zhang Nan.	Date: 01-12-2022, please send questions/comments/bugs to n-zhang19@mails.tsinghua.edu.cn
//
// The translated C++ source code accompanies the paper below.  Please cite appropriately when using the code.
// [1] Nitin Arora, and Ryan P. Russell, A Fast and Robust Multiple Revolution Lambert Algorithm Using a Cosine Transformation. In:AAS/AIAA Astrodynamics Specialist Conference, AAS 13-728, Hilton Head, SC (2013)
// [2] Ryan P. Russell, On the Solution to Every Lambert Problem, Celestial Mechanics and Dynamical Astronomy, 2019, https://dx.doi.org/10.1007/s10569-019-9927-z
// [3] Nitin Arora, Ryan P. Russell, Nathan Strange, and David Ottesen, Partial Derivatives of the Solution to the Lambert Boundary Value Problem, Journal of Guidance, Control, and Dynamics, 2015, https://doi.org/10.2514/1.G001030
//
// This code contains a classic Lambert problem multi-rev solver using the vercosine Lambert formulation and the initial guess.
// The initial guess adopts Ref[1] and other parts of the solver adopt Ref[2].
// This code is mainly different from Ref[1] and [2] in two respects.
// First, this code accurately calculates the minimum time Tmin in the multi-rev lambert problem in order to determine whether the input transfer time TOF is appropriate.
// Second, this code adopt both third-order householder update and binary search in order to ensure the convergence of the iteration.
//
// Note: lambert problems with low normalized transfer time suffer from poor precision. The rectilinear case can not be handled.
//*********************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************

/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院航天动力学与控制实验室
* 联系人: 武迪 1522620129@qq.com
* 文件名: Lambert.cpp
* 内容简述：包括张楠翻译的Ryan P. Russell组的lambert求解程序（效率较高，无法求解共线情况）及其梯度求解程序
*           蒋老师翻译Izzo组的lambert求解程序（较为通用）
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2022-09-15    武迪     创建文件并添加部分注释
****************************************************************************/

#include "Lambert.h"
#include <cmath>
#include <VecMat.h>

#include "OrbitMath.h"

/***********************************Ryan P. Russell组的Lambert和梯度程序**********************************/

const double ivLamThresh_parabolaBand = 0.02;
//nominal value 0.02  for double precision compiling; used to avoid singularity, for the series solution near parabola, we truncate at order 8, for double the max distance is .02 to meet 16 digits, for quad it is .0002 to meet ~33 digits, see mapleKlam_v2b.mw
const double ivLamThresh_zeroBand = 0.02;
//nominal value 0.02  for double precision compiling; used to save precision, not to avoid singularity, for minor precision problem of acos(1 + k ^ 2) when k is small, we truncate at order 8, for double the max distance is .02 to meet 16 digits, for quad it is .0002 to meet ~33 digits, see mapleKlam_v2b.mw
const double ivLamThresh_bigKiterate = 1000.0;
//nominal value 1.e3  for double precision compiling; used to save precision, not to avoid singularity, series solution for TOF / S in terms of 1 / k, for quad it should be ~1.d6, see kHugeFuncSeries.mw
const double ivLamThresh_littlePiterate = 0.1;
//nominal value 0.1	  for double precision compiling; used to save precision, not to avoid singularity, for iterating on the p(p = 1 - k * tau) variable instead of k when p is small, for quad i have not tested or tuned for best value
const double ivLamThresh_alternateTau = 0.01;
//nominal value 0.01  for double precision compiling; used to save precision, not to avoid singularity, for computing tau to higher precision when theta is close to pi(and only useful for low Tbar cases, i.e.N = 0 hyperbola); i.e.when geom% OnePctheta < ivLamThresh_alternateTau, for quad i have not tested or tuned for best value
const double ivLamThresh_kClosetoMsqrt2 = 0.01;
//nominal value 0.01  for double precision compiling; used to save precision, not to avoid singularity, for computing W to higher precision when k is close to - sqrt(2), it doesn't seem to affect accuracy of output (unsure why), but it does improve smoothness of W, for quad i have not tested or tuned for best value, see kHugeFuncSeries.mw
const bool ivLamParam_usePrecisionTricks = true;
//nominal value T; if set, then all the available tricks are used, otherwise not.  Mainly useful for N=0.  Very little, if any, performance penalty for using.  Nominally set to T, mainly for debugging/developing

// ivLam_Ebi0vec is the first 20 values of DeltaEbi when tau is zero. DeltaEbi is the change in eccentric anomaly of the minimum times for i-rev solutions.
const double ivLam_Ebi0vec[20] = {
	2.848574, 2.969742, 3.019580, 3.046927, 3.064234,
	3.076182, 3.084929, 3.091610, 3.096880, 3.101145,
	3.104666, 3.107623, 3.110142, 3.112312, 3.114203,
	3.115864, 3.117335, 3.118646, 3.119824, 3.120886
};

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
int lambert(double* v1vec, double* v2vec, double& a, double& e, const double* R1, const double* R2
            , double TOF, double mu, int way, int N, int branch, int Maxiter, double tol)
{
	//1, Successfully; 2, poor precision because transfer time too small
	//no solution: -1, negative time; -2, the rectilinear case; -3, TOF<Tmin; -4, not to be converging
	int flag = 1;
	if (TOF <= 0.0)
	{
		//		cout<<"Negative time as input"<<endl;
		flag = -1;
		return flag;
	}

	double r1vec[3], r2vec[3];
	for (int i = 0; i < 3; i++)
	{
		r1vec[i] = R1[i];
		r2vec[i] = R2[i];
	}

	//Non dimensional units
	double R = sqrt(r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1] + r1vec[2] * r1vec[2]);
	double V = sqrt(mu / R);
	double T = R / V;

	//working with non-dimensional radii and time-of-flight
	for (int i = 0; i < 3; i++)
	{
		r1vec[i] /= R;
		r2vec[i] /= R;
	}
	TOF /= T;

	//direction = 1 when theta < pi; direction = -1 when theta > pi;
	int direction = 1;
	if (r1vec[0] * r2vec[1] - r1vec[1] * r2vec[0] < 0.0)
		direction = -direction;
	if (way == 1)
		direction = -direction;

	double S, tofbyS, tau;
	double r1, r2, OneByr1, OneByr2, r1pr2, r1r2, ctheta, OnePctheta, OneMctheta;

	r1 = sqrt(r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1] + r1vec[2] * r1vec[2]);
	r2 = sqrt(r2vec[0] * r2vec[0] + r2vec[1] * r2vec[1] + r2vec[2] * r2vec[2]);
	OneByr1 = 1.0 / r1;
	OneByr2 = 1.0 / r2;
	r1pr2 = r1 + r2;
	r1r2 = r1 * r2;
	ctheta = (r1vec[0] * r2vec[0] + r1vec[1] * r2vec[1] + r1vec[2] * r2vec[2]) * OneByr1 * OneByr2;
	OnePctheta = ctheta + 1.0;
	OneMctheta = 1.0 - ctheta;

	//compute tau, there are two ways, the second way is nominal, cheaperand doesn't divide by zero, it works fine for all cases, but looses some precision when theta is close to pi
	if (ivLamParam_usePrecisionTricks && (OnePctheta < ivLamThresh_alternateTau))
	//alternate way to save precision when theta is close to pi
	{
		double r1crossr2[3];
		r1crossr2[0] = r1vec[1] * r2vec[2] - r1vec[2] * r2vec[1];
		r1crossr2[1] = r1vec[2] * r2vec[0] - r1vec[0] * r2vec[2];
		r1crossr2[2] = r1vec[0] * r2vec[1] - r1vec[1] * r2vec[0];
		double sthetar1r2 = sqrt(
			r1crossr2[0] * r1crossr2[0] + r1crossr2[1] * r1crossr2[1] + r1crossr2[2] * r1crossr2[2]);
		//sinTheta * r1 * r2 = norm(r1 cross r2)
		tau = direction * sqrt(1.0 / (OneMctheta * r1r2)) / r1pr2 * sthetar1r2;
	}
	else //nominal way, see note above
		tau = direction * sqrt(r1r2 * OnePctheta) / r1pr2;

	if (fabs(tau) < 1.0e-12 || OnePctheta < 1.0e-12)
	{
		//		cout<<"the rectilinear case"<<endl;
		flag = -2;
		return flag;
	}

	S = r1pr2 * sqrt(r1pr2);

	//note this has a mu in it but here we assume its unity
	tofbyS = TOF / S;

	double ksol, kmin, kmax;
	int getinitialKflag = lambert_getinitialK(ksol, kmin, kmax, direction, N, tau, tofbyS, branch);
	if (getinitialKflag == 0)
	{
		//		cout<<"TOF is less than the minimum time"<<endl;
		flag = -3;
		return flag;
	}

	int iters = 0;
	double dvar, dvar0, dvar1, dvar2, Func0, Func1, Func2, Func3, p, pmin, pmax, W0, W1, W2, W3, W4;

	p = 1.0 - ksol * tau;
	dvar = 1.0;
	Func0 = 1.0;

	if (ivLamParam_usePrecisionTricks && p < ivLamThresh_littlePiterate)
	{
		// the root-solve variable can be switched to p instead of k to order to fix the case of p is close to zero

		if (tau > 0.0)
		{
			pmax = 1.0 - kmin * tau;
			pmin = 1.0 - kmax * tau;
		}
		else
		{
			pmax = 1.0 - kmax * tau;
			pmin = 1.0 - kmin * tau;
		}

		while ((fabs(Func0 / tofbyS) > tol && fabs(dvar / p) > tol) && (iters < Maxiter))
		{
			iters++;

			lambert_getFuncAndKders(Func0, Func1, Func2, Func3, W0, W1, W2, W3, W4, ksol, p, N, tau, tofbyS);

			Func1 = Func1 / (-tau);
			Func2 = Func2 / (tau * tau);
			Func3 = Func3 / (-tau * tau * tau);

			if (Func0 * Func1 > 0.0)
				pmax = p;
			else
				pmin = p;

			lambert_getCorrection(dvar0, dvar1, dvar2, Func0, Func1, Func2, Func3);
			dvar = dvar0 + dvar1 + dvar2;

			if (!(fabs(Func0 / tofbyS) > tol && fabs(dvar / p) > tol))
				break;

			p = p + dvar;
			ksol = (1.0 - p) / tau;
			//not possible to divide by zero, as iterateOnP is false if tau = 0, since p = 1 int that case

			// if p is not in (pmin, pmax), binary search is adopted
			if (p <= pmin || p >= pmax)
			{
				p = (pmax + pmin) / 2.0;
				ksol = (1.0 - p) / tau;
				dvar = (pmax - pmin) / 2.0;
			}
		}
	}
	else
	{
		while ((fabs(Func0) > tol && fabs(dvar) > tol) && (iters < Maxiter))
		{
			iters++;

			lambert_getFuncAndKders(Func0, Func1, Func2, Func3, W0, W1, W2, W3, W4, ksol, p, N, tau, tofbyS);
			if (Func0 * Func1 > 0.0)
				kmax = ksol;
			else
				kmin = ksol;


			lambert_getCorrection(dvar0, dvar1, dvar2, Func0, Func1, Func2, Func3);
			dvar = dvar0 + dvar1 + dvar2;

			if (!(fabs(Func0) > tol && fabs(dvar) > tol))
				break;

			ksol = ksol + dvar;
			p = 1.0 - ksol * tau;

			// if ksol is not in (kmin, kmax), binary search is adopted
			if (ksol <= kmin || ksol >= kmax)
			{
				ksol = (kmax + kmin) / 2.0;
				p = 1.0 - ksol * tau;
				dvar = (kmax - kmin) / 2.0;
			}
		}
	}


	if (iters >= Maxiter)
	{
		//		cout<<"Solution does not seem to be converging"<<endl;
		flag = -4;
		return flag;
	}

	//getVelocityFromP
	double sqrtp, pr12, f, g, gdot, onebyg;
	sqrtp = sqrt(p);
	pr12 = p * r1pr2;
	f = 1.0 - pr12 * OneByr1;
	g = S * tau * sqrtp;
	gdot = 1.0 - pr12 * OneByr2;
	onebyg = 1.0 / g;

	for (int i = 0; i < 3; i++)
	{
		v1vec[i] = (r2vec[i] - f * r1vec[i]) * onebyg;
		v2vec[i] = (gdot * r2vec[i] - r1vec[i]) * onebyg;
	}

	if (fabs(p) < 1.0e-6)
		flag = 2;

	double Semilatus_rectum, sqrt2;
	sqrt2 = sqrt(2);
	Semilatus_rectum = r1 * r2 * OneMctheta / p / r1pr2;

	if (fabs(ksol - sqrt2) < 1.0e-14)
	{
		a = 0.5 * Semilatus_rectum * R;
		e = 1.0;
	}
	else if (ksol < sqrt2)
	{
		a = p * r1pr2 / (2.0 - ksol * ksol);
		e = sqrt(1.0 - Semilatus_rectum / a);
		a = a * R;
	}
	else
	{
		a = p * r1pr2 / (ksol * ksol - 2.0);;
		e = sqrt(1.0 + Semilatus_rectum / a);
		a = a * R;
	}

	for (int i = 0; i < 3; i++)
	{
		v1vec[i] *= V;
		v2vec[i] *= V;
	}

	return flag;
}

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
            double TOF, double mu, int way, int N, int branch, int Maxiter, double tol)
{
	//1, Successfully; 2, poor precision because transfer time too small
	//no solution: -1, negative time; -2, the rectilinear case; -3, TOF<Tmin; -4, not to be converging
	int flag = 1;
	if (TOF <= 0.0)
	{
		//		cout<<"Negative time as input"<<endl;
		flag = -1;
		return flag;
	}

	double r1vec[3], r2vec[3];
	for (int i = 0; i < 3; i++)
	{
		r1vec[i] = R1[i];
		r2vec[i] = R2[i];
	}

	//Non dimensional units
	double R = sqrt(r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1] + r1vec[2] * r1vec[2]);
	double V = sqrt(mu / R);
	double T = R / V;

	//working with non-dimensional radii and time-of-flight
	for (int i = 0; i < 3; i++)
	{
		r1vec[i] /= R;
		r2vec[i] /= R;
	}
	TOF /= T;

	//direction = 1 when theta < pi; direction = -1 when theta > pi;
	int direction = 1;
	if (r1vec[0] * r2vec[1] - r1vec[1] * r2vec[0] < 0.0)
		direction = -direction;
	if (way == 1)
		direction = -direction;

	double S, tofbyS, tau;
	double r1, r2, OneByr1, OneByr2, r1pr2, r1r2, ctheta, stheta, OnePctheta, OneMctheta;

	r1 = sqrt(r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1] + r1vec[2] * r1vec[2]);
	r2 = sqrt(r2vec[0] * r2vec[0] + r2vec[1] * r2vec[1] + r2vec[2] * r2vec[2]);
	OneByr1 = 1.0 / r1;
	OneByr2 = 1.0 / r2;
	r1pr2 = r1 + r2;
	r1r2 = r1 * r2;
	ctheta = (r1vec[0] * r2vec[0] + r1vec[1] * r2vec[1] + r1vec[2] * r2vec[2]) * OneByr1 * OneByr2;
	stheta = sqrt(1.0 - ctheta * ctheta) * direction;
	OnePctheta = ctheta + 1.0;
	OneMctheta = 1.0 - ctheta;

	//compute tau, there are two ways, the second way is nominal, cheaperand doesn't divide by zero, it works fine for all cases, but looses some precision when theta is close to pi
	if (ivLamParam_usePrecisionTricks && (OnePctheta < ivLamThresh_alternateTau))
	//alternate way to save precision when theta is close to pi
	{
		double r1crossr2[3];
		r1crossr2[0] = r1vec[1] * r2vec[2] - r1vec[2] * r2vec[1];
		r1crossr2[1] = r1vec[2] * r2vec[0] - r1vec[0] * r2vec[2];
		r1crossr2[2] = r1vec[0] * r2vec[1] - r1vec[1] * r2vec[0];
		double sthetar1r2 = sqrt(
			r1crossr2[0] * r1crossr2[0] + r1crossr2[1] * r1crossr2[1] + r1crossr2[2] * r1crossr2[2]);
		//sinTheta * r1 * r2 = norm(r1 cross r2)
		tau = direction * sqrt(1.0 / (OneMctheta * r1r2)) / r1pr2 * sthetar1r2;
	}
	else //nominal way, see note above
		tau = direction * sqrt(r1r2 * OnePctheta) / r1pr2;

	if (fabs(tau) < 1.0e-12 || OnePctheta < 1.0e-12)
	{
		//		cout<<"the rectilinear case"<<endl;
		flag = -2;
		return flag;
	}

	S = r1pr2 * sqrt(r1pr2);

	//note this has a mu in it but here we assume its unity
	tofbyS = TOF / S;

	double ksol, kmin, kmax;
	int getinitialKflag = lambert_getinitialK(ksol, kmin, kmax, direction, N, tau, tofbyS, branch);
	if (getinitialKflag == 0)
	{
		//		cout<<"TOF is less than the minimum time"<<endl;
		flag = -3;
		return flag;
	}

	int iters = 0;
	double dvar, dvar0, dvar1, dvar2, Func0, Func1, Func2, Func3, p, pmin, pmax, W0, W1, W2, W3, W4;

	p = 1.0 - ksol * tau;
	dvar = 1.0;
	Func0 = 1.0;

	if (ivLamParam_usePrecisionTricks && p < ivLamThresh_littlePiterate)
	{
		// the root-solve variable can be switched to p instead of k to order to fix the case of p is close to zero

		if (tau > 0.0)
		{
			pmax = 1.0 - kmin * tau;
			pmin = 1.0 - kmax * tau;
		}
		else
		{
			pmax = 1.0 - kmax * tau;
			pmin = 1.0 - kmin * tau;
		}

		while ((fabs(Func0 / tofbyS) > tol && fabs(dvar / p) > tol) && (iters < Maxiter))
		{
			iters++;

			lambert_getFuncAndKders(Func0, Func1, Func2, Func3, W0, W1, W2, W3, W4, ksol, p, N, tau, tofbyS);

			Func1 = Func1 / (-tau);
			Func2 = Func2 / (tau * tau);
			Func3 = Func3 / (-tau * tau * tau);

			if (Func0 * Func1 > 0.0)
				pmax = p;
			else
				pmin = p;

			lambert_getCorrection(dvar0, dvar1, dvar2, Func0, Func1, Func2, Func3);
			dvar = dvar0 + dvar1 + dvar2;

			if (!(fabs(Func0 / tofbyS) > tol && fabs(dvar / p) > tol))
				break;

			p = p + dvar;
			ksol = (1.0 - p) / tau;
			//not possible to divide by zero, as iterateOnP is false if tau = 0, since p = 1 int that case

			// if p is not in (pmin, pmax), binary search is adopted
			if (p <= pmin || p >= pmax)
			{
				p = (pmax + pmin) / 2.0;
				ksol = (1.0 - p) / tau;
				dvar = (pmax - pmin) / 2.0;
			}
		}
	}
	else
	{
		while ((fabs(Func0) > tol && fabs(dvar) > tol) && (iters < Maxiter))
		{
			iters++;

			lambert_getFuncAndKders(Func0, Func1, Func2, Func3, W0, W1, W2, W3, W4, ksol, p, N, tau, tofbyS);
			if (Func0 * Func1 > 0.0)
				kmax = ksol;
			else
				kmin = ksol;

			lambert_getCorrection(dvar0, dvar1, dvar2, Func0, Func1, Func2, Func3);
			dvar = dvar0 + dvar1 + dvar2;

			if (!(fabs(Func0) > tol && fabs(dvar) > tol))
				break;

			ksol = ksol + dvar;
			p = 1.0 - ksol * tau;

			// if ksol is not in (kmin, kmax), binary search is adopted
			if (ksol <= kmin || ksol >= kmax)
			{
				ksol = (kmax + kmin) / 2.0;
				p = 1.0 - ksol * tau;
				dvar = (kmax - kmin) / 2.0;
			}
		}
	}


	if (iters >= Maxiter)
	{
		//		cout<<"Solution does not seem to be converging"<<endl;
		flag = -4;
		return flag;
	}

	//getVelocityFromP
	double sqrtp, pr12, f, g, gdot, onebyg;
	sqrtp = sqrt(p);
	pr12 = p * r1pr2;
	f = 1.0 - pr12 * OneByr1;
	g = S * tau * sqrtp;
	gdot = 1.0 - pr12 * OneByr2;
	onebyg = 1.0 / g;

	for (int i = 0; i < 3; i++)
	{
		v1vec[i] = (r2vec[i] - f * r1vec[i]) * onebyg;
		v2vec[i] = (gdot * r2vec[i] - r1vec[i]) * onebyg;
	}

	if (fabs(p) < 1.0e-6)
		flag = 2;

	double Semilatus_rectum, sqrt2;
	sqrt2 = sqrt(2);
	Semilatus_rectum = r1 * r2 * OneMctheta / p / r1pr2;

	if (fabs(ksol - sqrt2) < 1.0e-14)
	{
		a = 0.5 * Semilatus_rectum * R;
		e = 1.0;
	}
	else if (ksol < sqrt2)
	{
		a = p * r1pr2 / (2.0 - ksol * ksol);
		e = sqrt(1.0 - Semilatus_rectum / a);
		a = a * R;
	}
	else
	{
		a = p * r1pr2 / (ksol * ksol - 2.0);;
		e = sqrt(1.0 + Semilatus_rectum / a);
		a = a * R;
	}

	for (int i = 0; i < 3; i++)
	{
		v1vec[i] *= V;
		v2vec[i] *= V;
	}

	// Partial Derivatives calculation
	// x = [r1; r2; theta; TOF];
	// P = [r1vec; r2vec; TOF];
	// q = [f; g; dg];

	double PfPr1, PfPr2, PgdotPr2, PfPtau, PgdotPtau, PgPtau, PgPS, OneByr1pr2, tauBytwo,
	       tauByr1pr2, PSPr1, PSPr2, PtauPr1, PtauPr2, PtauPtheta, c, gamma, OneByPTOFPk, TOFByPTOFPk,
	       TOFByPTOFPk_S, TOFByPTOFPk_tau_gamma, PfPk, PgdotPk, PgPk, onebyg_T, PgPSmPSPr,
	       PkPx[4], r2vecUnit[3], PthetaPP[6], PkPP[7], PtauPP[6], DfDP[7], DgdotDP[7], DgDP[7];

	// PkPx calculation
	OneByr1pr2 = 1.0 / r1pr2;
	tauBytwo = 0.5 * tau;
	tauByr1pr2 = tau * OneByr1pr2;
	PSPr1 = 1.5 * S * OneByr1pr2;
	PSPr2 = PSPr1;
	PtauPr1 = tauBytwo * OneByr1 - tauByr1pr2;
	PtauPr2 = tauBytwo * OneByr2 - tauByr1pr2;
	PtauPtheta = -tauBytwo * stheta / OnePctheta;
	c = p / tau;
	if (ivLamParam_usePrecisionTricks && ksol > ivLamThresh_bigKiterate)
		lambert_getD4W(W0, W1, W2, W3, W4, ksol, N);
	gamma = (1.0 + c * W0) / (1.0 - 0.5 * ksol / c * (3 * c * W0 + 1.0));
	OneByPTOFPk = 1.0 / (-0.5 * TOF / c + S * tau * sqrt(c * tau) * (W1 * c - W0));
	TOFByPTOFPk = TOF * OneByPTOFPk;
	TOFByPTOFPk_S = TOFByPTOFPk / S;
	TOFByPTOFPk_tau_gamma = TOFByPTOFPk / tau / gamma;
	PkPx[0] = -TOFByPTOFPk_S * PSPr1 - TOFByPTOFPk_tau_gamma * PtauPr1;
	PkPx[1] = -TOFByPTOFPk_S * PSPr2 - TOFByPTOFPk_tau_gamma * PtauPr2;
	PkPx[2] = -TOFByPTOFPk_tau_gamma * PtauPtheta;
	PkPx[3] = OneByPTOFPk;

	// PkPP, PSPP, PtauPP
	r2vecUnit[0] = r2vec[0] / r2;
	r2vecUnit[1] = r2vec[1] / r2;
	r2vecUnit[2] = r2vec[2] / r2;
	for (int i = 0; i < 3; i++)
	{
		PthetaPP[i] = (r1vec[i] * ctheta - r2vecUnit[i]) * OneByr1 / stheta;
		PthetaPP[3 + i] = (r2vecUnit[i] * ctheta - r1vec[i]) * OneByr2 / stheta;
	}
	PkPP[0] = PkPx[0] * r1vec[0] + PkPx[2] * PthetaPP[0];
	PkPP[1] = PkPx[0] * r1vec[1] + PkPx[2] * PthetaPP[1];
	PkPP[2] = PkPx[0] * r1vec[2] + PkPx[2] * PthetaPP[2];
	PkPP[3] = PkPx[1] * r2vecUnit[0] + PkPx[2] * PthetaPP[3];
	PkPP[4] = PkPx[1] * r2vecUnit[1] + PkPx[2] * PthetaPP[4];
	PkPP[5] = PkPx[1] * r2vecUnit[2] + PkPx[2] * PthetaPP[5];
	PkPP[6] = PkPx[3];
	PtauPP[0] = PtauPr1 * r1vec[0] + PtauPtheta * PthetaPP[0];
	PtauPP[1] = PtauPr1 * r1vec[1] + PtauPtheta * PthetaPP[1];
	PtauPP[2] = PtauPr1 * r1vec[2] + PtauPtheta * PthetaPP[2];
	PtauPP[3] = PtauPr2 * r2vecUnit[0] + PtauPtheta * PthetaPP[3];
	PtauPP[4] = PtauPr2 * r2vecUnit[1] + PtauPtheta * PthetaPP[4];
	PtauPP[5] = PtauPr2 * r2vecUnit[2] + PtauPtheta * PthetaPP[5];

	// DfDP, DgdotDP, DgDP
	PfPr1 = (1.0 - f) * (OneByr1 - OneByr1pr2);
	PfPr2 = -(1.0 - f) * OneByr1pr2;
	PgdotPr2 = (1.0 - gdot) * (OneByr2 - OneByr1pr2);
	PfPtau = ksol * r1pr2 / r1;
	PgdotPtau = ksol * r1pr2 / r2;
	PgPtau = g / tau * (1.0 - 0.5 * ksol / c);
	PgPS = g / S;
	PfPk = tau * r1pr2 / r1;
	PgdotPk = tau * r1pr2 / r2;
	PgPk = -0.5 * g / c;
	DfDP[0] = PfPr1 * r1vec[0] + PfPtau * PtauPP[0] + PfPk * PkPP[0];
	DfDP[1] = PfPr1 * r1vec[1] + PfPtau * PtauPP[1] + PfPk * PkPP[1];
	DfDP[2] = PfPr1 * r1vec[2] + PfPtau * PtauPP[2] + PfPk * PkPP[2];
	DfDP[3] = PfPr2 * r2vecUnit[0] + PfPtau * PtauPP[3] + PfPk * PkPP[3];
	DfDP[4] = PfPr2 * r2vecUnit[1] + PfPtau * PtauPP[4] + PfPk * PkPP[4];
	DfDP[5] = PfPr2 * r2vecUnit[2] + PfPtau * PtauPP[5] + PfPk * PkPP[5];
	DfDP[6] = PfPk * PkPP[6];
	DgdotDP[3] = PgdotPr2 * r2vecUnit[0] + PgdotPtau * PtauPP[3] + PgdotPk * PkPP[3];
	DgdotDP[4] = PgdotPr2 * r2vecUnit[1] + PgdotPtau * PtauPP[4] + PgdotPk * PkPP[4];
	DgdotDP[5] = PgdotPr2 * r2vecUnit[2] + PgdotPtau * PtauPP[5] + PgdotPk * PkPP[5];
	DgdotDP[6] = PgdotPk * PkPP[6];
	PgPSmPSPr = PgPS * PSPr1;
	DgDP[0] = PgPSmPSPr * r1vec[0] + PgPtau * PtauPP[0] + PgPk * PkPP[0];
	DgDP[1] = PgPSmPSPr * r1vec[1] + PgPtau * PtauPP[1] + PgPk * PkPP[1];
	DgDP[2] = PgPSmPSPr * r1vec[2] + PgPtau * PtauPP[2] + PgPk * PkPP[2];
	DgDP[3] = PgPSmPSPr * r2vecUnit[0] + PgPtau * PtauPP[3] + PgPk * PkPP[3];
	DgDP[4] = PgPSmPSPr * r2vecUnit[1] + PgPtau * PtauPP[4] + PgPk * PkPP[4];
	DgDP[5] = PgPSmPSPr * r2vecUnit[2] + PgPtau * PtauPP[5] + PgPk * PkPP[5];
	DgDP[6] = PgPk * PkPP[6];

	//PD
	onebyg_T = onebyg / T;

	PD[0] = (-r1vec[0] * DfDP[0] - f - v1vec[0] / V * DgDP[0]) * onebyg_T;
	PD[1] = (-r1vec[0] * DfDP[1] - v1vec[0] / V * DgDP[1]) * onebyg_T;
	PD[2] = (-r1vec[0] * DfDP[2] - v1vec[0] / V * DgDP[2]) * onebyg_T;
	PD[3] = (1.0 - r1vec[0] * DfDP[3] - v1vec[0] / V * DgDP[3]) * onebyg_T;
	PD[4] = (-r1vec[0] * DfDP[4] - v1vec[0] / V * DgDP[4]) * onebyg_T;
	PD[5] = (-r1vec[0] * DfDP[5] - v1vec[0] / V * DgDP[5]) * onebyg_T;
	PD[6] = -(-r1vec[0] * DfDP[6] * V - v1vec[0] * DgDP[6]) * onebyg_T;
	PD[7] = -PD[6];
	PD[8] = PD[1];
	PD[9] = (-r1vec[1] * DfDP[1] - f - v1vec[1] / V * DgDP[1]) * onebyg_T;
	PD[10] = (-r1vec[1] * DfDP[2] - v1vec[1] / V * DgDP[2]) * onebyg_T;
	PD[11] = (-r1vec[1] * DfDP[3] - v1vec[1] / V * DgDP[3]) * onebyg_T;
	PD[12] = (1.0 - r1vec[1] * DfDP[4] - v1vec[1] / V * DgDP[4]) * onebyg_T;
	PD[13] = (-r1vec[1] * DfDP[5] - v1vec[1] / V * DgDP[5]) * onebyg_T;
	PD[14] = -(-r1vec[1] * DfDP[6] * V - v1vec[1] * DgDP[6]) * onebyg_T;
	PD[15] = -PD[14];
	PD[16] = PD[2];
	PD[17] = PD[10];
	PD[18] = (-r1vec[2] * DfDP[2] - f - v1vec[2] / V * DgDP[2]) * onebyg_T;
	PD[19] = (-r1vec[2] * DfDP[3] - v1vec[2] / V * DgDP[3]) * onebyg_T;
	PD[20] = (-r1vec[2] * DfDP[4] - v1vec[2] / V * DgDP[4]) * onebyg_T;
	PD[21] = (1.0 - r1vec[2] * DfDP[5] - v1vec[2] / V * DgDP[5]) * onebyg_T;
	PD[22] = -(-r1vec[2] * DfDP[6] * V - v1vec[2] * DgDP[6]) * onebyg_T;
	PD[23] = -PD[22];
	PD[24] = -PD[3];
	PD[25] = -PD[11];
	PD[26] = -PD[19];
	PD[27] = -(-r2vec[0] * DgdotDP[3] - gdot + v2vec[0] / V * DgDP[3]) * onebyg_T;
	PD[28] = -(-r2vec[0] * DgdotDP[4] + v2vec[0] / V * DgDP[4]) * onebyg_T;
	PD[29] = -(-r2vec[0] * DgdotDP[5] + v2vec[0] / V * DgDP[5]) * onebyg_T;
	PD[30] = (-r2vec[0] * DgdotDP[6] * V + v2vec[0] * DgDP[6]) * onebyg_T;
	PD[31] = -PD[30];
	PD[32] = -PD[4];
	PD[33] = -PD[12];
	PD[34] = -PD[20];
	PD[35] = PD[28];
	PD[36] = -(-r2vec[1] * DgdotDP[4] - gdot + v2vec[1] / V * DgDP[4]) * onebyg_T;
	PD[37] = -(-r2vec[1] * DgdotDP[5] + v2vec[1] / V * DgDP[5]) * onebyg_T;
	PD[38] = (-r2vec[1] * DgdotDP[6] * V + v2vec[1] * DgDP[6]) * onebyg_T;
	PD[39] = -PD[38];
	PD[40] = -PD[5];
	PD[41] = -PD[13];
	PD[42] = -PD[21];
	PD[43] = PD[29];
	PD[44] = PD[37];
	PD[45] = -(-r2vec[2] * DgdotDP[5] - gdot + v2vec[2] / V * DgDP[5]) * onebyg_T;
	PD[46] = (-r2vec[2] * DgdotDP[6] * V + v2vec[2] * DgDP[6]) * onebyg_T;
	PD[47] = -PD[46];

	return flag;
}

//this routine returns the W functionand its derivatives with respect to k up to order 1 - 4
void lambert_getD4W(double& W0, double& W1, double& W2, double& W3, double& W4, double k, int N)
{
	double ksqm1, m, nu, t2, t4, t6, t8, t10, t12, t14,
	       onebym, ksq, kps2, tNp, tb1, tb2, tb3, tb9, tb10;
	double kregion;

	double ivLamThresh_zeroBand2 = ivLamThresh_zeroBand * ivLamThresh_zeroBand;
	double SQRT2 = sqrt(2.0);
	double pi = 3.14159265358979323846264338327950;
	double twopi = 2.0 * pi;
	double twopiN = twopi * N;

	double t3, t18, tnpp;

	tnpp = twopiN + pi;

	W0 = 0.0;

	ksq = k * k;
	nu = k - SQRT2;

	kregion = 0;
	if (ksq <= ivLamThresh_zeroBand2)
		kregion = 0; //close to zero series(ellipse still)
	else if (N == 0)
		if (fabs(nu) < ivLamThresh_parabolaBand)
			kregion = 2; //close to sqrt2 series(parab)
		else if (ksq > 2.0)
			kregion = 3; //hyperbola not picked up by parab
		else
			kregion = 1; //ellipse not picked up by parab
	else
		kregion = 1; //ellipse not picked up by parab

	if (kregion == 2)
	{
		t2 = nu * nu; //2
		t4 = t2 * nu; //3
		t6 = t2 * t2; //4
		t8 = t6 * nu; //5
		t10 = t6 * t2; //6
		t12 = t6 * t4; //7
		t14 = t6 * t6; //8

		W0 = 0.47140452079103168293389624140323e0 - 0.20000000000000000000000000000000e0 * nu
			+ 0.80812203564176859931525069954840e-1 * t2 - 0.31746031746031746031746031746032e-1 * t4
			+ 0.12244273267299524232049253023461e-1 * t6 - 0.46620046620046620046620046620047e-2 * t8
			+ 0.17581520588942906589609183828558e-2 * t10 - 0.65816536404771698889345948169478e-3 * t12
			+ 0.24494378529487021564470999141955e-3 * t14;
		W1 = -0.20000000000000000000000000000000e0 + 0.16162440712835371986305013990967e0 * nu
			- 0.95238095238095238095238095238095e-1 * t2 + 0.48977093069198096928197012093843e-1 * t4
			- 0.23310023310023310023310023310023e-1 * t6 + 0.10548912353365743953765510297135e-1 * t8
			- 0.46071575483340189222542163718634e-2 * t10 + 0.19595502823589617251576799313564e-2 * t12;
		W2 = -0.19047619047619047619047619047619e0 + 0.29386255841518858156918207256306e0 * nu
			- 0.27972027972027972027972027972028e0 * t2 + 0.21097824706731487907531020594271e0 * t4
			- 0.13821472645002056766762649115590e0 * t6 + 0.82301111859076392456622557116969e-1 * t8;
		W3 = 0.29386255841518858156918207256306e0 - 0.55944055944055944055944055944056e0 * nu
			+ 0.63293474120194463722593061782812e0 * t2 - 0.55285890580008227067050596462361e0 * t4
			+ 0.41150555929538196228311278558484e0 * t6;
		W4 = 0.0; //not used
	}
	else
	{
		ksqm1 = ksq - 1.0;
		m = 1.0 - ksqm1;
		onebym = 1.0 / m;

		if ((kregion == 1))
		{
			if (k > 0.0)
				W0 = (twopiN + acos(ksqm1)) * sqrt(onebym * onebym * onebym) - k * onebym; //erase
			else
			{
				kps2 = k + SQRT2;
				//see comments where the thresh is defined, this is a smoothness / precision saving effort in computing W in extreme case when close to boundaryand kps2 is close to zero
				if (ivLamParam_usePrecisionTricks && (kps2 < ivLamThresh_kClosetoMsqrt2))
				{
					tNp = twopi + twopiN;
					tb1 = kps2 * kps2;
					tb2 = sqrt(kps2);
					tb3 = tb2 * tb1;
					tb9 = tb1 * tb1;
					tb10 = tb2 * tb9;
					W0 = 0.12110150049603174603174603174603e-7 / tb3 * (-0.38926398009946925989672338336519e8 * tb3
						- 0.16515072000000000000000000e8 * tb2 * kps2 * tb1 - 0.1976320000000000000000000e7 * tb10
						* (kps2 + 0.24532575164897006338798206469711e1) - 0.18246749067162621557658908595243e7 * tb10
						+ 0.25959796716951899525909665607350e6 * (tb9 + 0.64646464646464646464646464646465e1 * tb1
							+ 0.35463203463203463203463203463203e2) * tNp * tb1 + 0.66750357442839860425810740303391e6 *
						tNp
						* (tb9 + 0.60952380952380952380952380952381e1 * tb1 + 0.26006349206349206349206349206349e2) *
						kps2
						- 0.645120e6 * tb2 * kps2 * tb9);
				}
				else
					W0 = (twopi + twopiN - acos(ksqm1)) * sqrt(onebym * onebym * onebym) - k * onebym;
			}
		}
		else if (kregion == 0) //ellipse close to zero k
		{
			t3 = ksq;
			t6 = t3 * k;
			t8 = t3 * t3;
			t18 = t8 * t8;
			W0 = 0.3535533905932738e0 * tnpp - 0.1e1 * k + 0.2651650429449553e0 * tnpp * t3
				- 0.6666666666666667e0 * t6 + 0.1657281518405971e0 * tnpp * t8
				- 0.4000000000000000e0 * t8 * k + 0.9667475524034829e-1 * tnpp * t8 * t3
				- 0.2285714285714286e0 * t8 * t6 + 0.5437954982269591e-1 * tnpp * t18;
		}
		else if (kregion == 3) //full hyperbola
			W0 = (-log(ksqm1 + sqrt(ksqm1 * ksqm1 - 1.0))) * sqrt(-onebym * onebym * onebym) - k * onebym;

		t2 = 3.e0 * W0;
		W1 = (t2 * k - 2.e0) * onebym;
		W2 = (5.e0 * W1 * k + t2) * (onebym);
		W3 = (7.e0 * W2 * k + 8.e0 * W1) * (onebym);
		W4 = (9.e0 * W3 * k + 15.e0 * W2) * (onebym);
	}
}

//only calculate W0, only used in getinitialK Function
void lambert_getW(double& W0, double k, int N)
{
	double ksqm1, m, nu, t2, t4, t6, t8, t10, t12, t14,
	       onebym, ksq, kps2, tNp, tb1, tb2, tb3, tb9, tb10;
	int kregion;

	double ivLamThresh_zeroBand2 = ivLamThresh_zeroBand * ivLamThresh_zeroBand;
	double SQRT2 = sqrt(2.0);
	double pi = 3.14159265358979323846264338327950;
	double twopi = 2.0 * pi;
	double twopiN = twopi * N;

	double t3, t18, tnpp;

	tnpp = twopiN + pi;

	W0 = 0.0;

	ksq = k * k;
	nu = k - SQRT2;

	kregion = 0;
	if (ksq <= ivLamThresh_zeroBand2)
		kregion = 0; //close to zero series(ellipse still)
	else if (N == 0)
		if (fabs(nu) < ivLamThresh_parabolaBand)
			kregion = 2; //close to sqrt2 series(parab)
		else if (ksq > 2.0)
			kregion = 3; //hyperbola not picked up by parab
		else
			kregion = 1; //ellipse not picked up by parab
	else
		kregion = 1; //ellipse not picked up by parab

	if (kregion == 2)
	{
		t2 = nu * nu; //2
		t4 = t2 * nu; //3
		t6 = t2 * t2; //4
		t8 = t6 * nu; //5
		t10 = t6 * t2; //6
		t12 = t6 * t4; //7
		t14 = t6 * t6; //8

		W0 = 0.47140452079103168293389624140323e0 - 0.20000000000000000000000000000000e0 * nu
			+ 0.80812203564176859931525069954840e-1 * t2 - 0.31746031746031746031746031746032e-1 * t4
			+ 0.12244273267299524232049253023461e-1 * t6 - 0.46620046620046620046620046620047e-2 * t8
			+ 0.17581520588942906589609183828558e-2 * t10 - 0.65816536404771698889345948169478e-3 * t12
			+ 0.24494378529487021564470999141955e-3 * t14;
	}
	else
	{
		ksqm1 = ksq - 1.0;
		m = 1.0 - ksqm1;
		onebym = 1.0 / m;

		if ((kregion == 1))
		{
			if (k > 0.0)
				W0 = (twopiN + acos(ksqm1)) * sqrt(onebym * onebym * onebym) - k * onebym; //erase
			else
			{
				kps2 = k + SQRT2;
				//see comments where the thresh is defined, this is a smoothness / precision saving effort in computing W in extreme case when close to boundaryand kps2 is close to zero
				if (ivLamParam_usePrecisionTricks && (kps2 < ivLamThresh_kClosetoMsqrt2))
				{
					tNp = twopi + twopiN;
					tb1 = kps2 * kps2;
					tb2 = sqrt(kps2);
					tb3 = tb2 * tb1;
					tb9 = tb1 * tb1;
					tb10 = tb2 * tb9;
					W0 = 0.12110150049603174603174603174603e-7 / tb3 * (-0.38926398009946925989672338336519e8 * tb3
						- 0.16515072000000000000000000e8 * tb2 * kps2 * tb1 - 0.1976320000000000000000000e7 * tb10
						* (kps2 + 0.24532575164897006338798206469711e1) - 0.18246749067162621557658908595243e7 * tb10
						+ 0.25959796716951899525909665607350e6 * (tb9 + 0.64646464646464646464646464646465e1 * tb1
							+ 0.35463203463203463203463203463203e2) * tNp * tb1 + 0.66750357442839860425810740303391e6 *
						tNp
						* (tb9 + 0.60952380952380952380952380952381e1 * tb1 + 0.26006349206349206349206349206349e2) *
						kps2
						- 0.645120e6 * tb2 * kps2 * tb9);
				}
				else
					W0 = (twopi + twopiN - acos(ksqm1)) * sqrt(onebym * onebym * onebym) - k * onebym;
			}
		}
		else if (kregion == 0) //ellipse close to zero k
		{
			t3 = ksq;
			t6 = t3 * k;
			t8 = t3 * t3;
			t18 = t8 * t8;
			W0 = 0.3535533905932738e0 * tnpp - 0.1e1 * k + 0.2651650429449553e0 * tnpp * t3
				- 0.6666666666666667e0 * t6 + 0.1657281518405971e0 * tnpp * t8
				- 0.4000000000000000e0 * t8 * k + 0.9667475524034829e-1 * tnpp * t8 * t3
				- 0.2285714285714286e0 * t8 * t6 + 0.5437954982269591e-1 * tnpp * t18;
		}
		else if (kregion == 3) //full hyperbola
			W0 = (-log(ksqm1 + sqrt(ksqm1 * ksqm1 - 1.0))) * sqrt(-onebym * onebym * onebym) - k * onebym;
	}
}

//if k is big, this is a precision saving version evaluation of the TOF/S equation and its partials
//the normal approach works fine, but leads to poor precision due to the (tau+pW) part of the TOF eq,
//so we use a series solution for the whole TOF eq.
void lambert_getFuncAndKdersBigK(double& Func0, double& Func1, double& Func2, double& Func3, double k, int N,
                                 double tau, double tofbyS)
{
	double t1, t3, t4, t5, t8, t11, t16, t32, t72, t13, t15, t14, t73, t74, t75, t70, t42, t43, t102, t103, t106, t113,
	       t115, t107;
	double dQ[4];
	double log2 = log(2.0);

	double p = 1.0 - k * tau;
	double minp = -p;

	t1 = k * k;
	t3 = minp + 1.0;
	t5 = minp * (t1 + 3.0);
	t4 = 1.0 / k;
	t8 = log(t4);
	t11 = log2;
	t13 = t1 * t1;
	t14 = t1 * k;
	t15 = t14 * tau;
	t16 = 2.0 * t15;

	t42 = t4 * t4;
	t72 = t42 * t42;
	dQ[0] = t72 * t4 * (t11 * t5 - 2.0 * t8 * t5 + 2.0 * t1 + t13 - t16 - 5.0 * t3 + 5.0);

	t102 = -minp;
	t103 = sqrt(t102); //sqrt(p)
	Func0 = t103 * dQ[0] - tofbyS;

	//-------------- first oder below -------------------
	t32 = 6.0 * t15;
	t73 = (-t16 + 3.0 * t1 - 12.0 * t3 + 15.0);
	t70 = (-2.0 * t8 + t11);
	t43 = t72 * t42;
	dQ[1] = t43 * (t70 * t73 - t13 + t32 - 8.0 * t1 + 26.0 * t3 - 31.0);

	t32 = dQ[0] * tau;
	t106 = 1.0 / t103; //1 / sqrt(p)
	Func1 = -t106 * t32 * 0.5 + t103 * dQ[1];

	//-------------- second oder below -------------------
	t74 = (t32 - 12.0 * t1 + 60.0 * t3 - 90.0);
	dQ[2] = t43 * t4 * (t70 * t74 + 2.0 * t13 - 22.0 * t15 + 38.0 * t1 - 154.0 * t3 + 216.0);

	t107 = t106 * t106; //1 / p
	t113 = t106 * t107; //1 / p ^ (3 / 2)
	t115 = tau * tau;
	Func2 = -t113 * dQ[0] * t115 * 0.250 - t106 * dQ[1] * tau + t103 * dQ[2];

	//-------------- third oder below -------------------
	t75 = (-24.0 * t15 + 60.0 * t1 - 360.0 * t3 + 630.0);
	dQ[3] = t72 * t72 * (t70 * t75 - 6.0 * t13 + 100.0 * t15 - 214.0 * t1 + 1044.0 * t3 - 1692.0);
	Func3 = -0.3750 * t113 * t107 * t32 * t115 - 0.750 * t113 * dQ[1] * t115 - 1.50 * t106 * dQ[2] * tau + t103 * dQ[3];
}

//the TOF function and its partials wrt k
void lambert_getFuncAndKders(double& Func0, double& Func1, double& Func2, double& Func3, double& W0, double& W1,
                             double& W2, double& W3, double& W4, double k, double p, int N, double tau, double tofbyS)
{
	//if k is big, this is a precision saving version evaluation of the TOF/S equation and its partials
	//the normal approach works fine, but leads to poor precision due to the (tau+pW) part of the TOF eq,
	//so we use a series solution for the whole TOF eq.
	if (ivLamParam_usePrecisionTricks && k > ivLamThresh_bigKiterate)
	//this is a series solution in 1/k, note it only happens when k is huge,
	//so its always a hyperbola and W is not computed, its hyperbola form is how the series is computed
		lambert_getFuncAndKdersBigK(Func0, Func1, Func2, Func3, k, N, tau, tofbyS);
	else
	{
		//this is the regular way to compute the TOF/S function and its partials wrt k
		double t7, p3, onebyRootp, onebyp, onebyp32, p3dw0, t7tau;
		double tau2, tau3;

		tau2 = tau * tau;
		tau3 = tau2 * tau;

		//double p = 1.0 - k * tau;
		double sqrtp = sqrt(p);

		//compute the W(k) function and its partials
		//double W0, W1, W2, W3, W4;
		lambert_getD4W(W0, W1, W2, W3, W4, k, N);

		//below compute the TOF/S function
		Func0 = sqrtp * (p * W0 + tau) - tofbyS;

		//now the function partials wrt k
		t7 = p * p;
		p3 = 3.0 * p;
		onebyRootp = 1.0 / sqrtp;
		Func1 = (-p3 * tau * W0 + 2.0 * t7 * W1 - tau2) * onebyRootp * 0.5;

		onebyp = onebyRootp * onebyRootp;
		onebyp32 = onebyRootp * onebyp;
		p3dw0 = p3 * W0;
		t7tau = t7 * tau;
		Func2 = (p3dw0 * tau2 + 4.0 * t7 * p * W2 - 12.0 * t7tau * W1 - tau3) * onebyp32 * 0.25;

		Func3 = (p3dw0 * tau3 + 18.0 * t7 * tau2 * W1 - 36.0 * t7tau * p * W2
			+ 8.0 * t7 * t7 * W3 - 3.0 * tau2 * tau2) * onebyp32 * onebyp * 0.125;
	}
}

//ith derivative of the function wrt the value
//ith order (isolated to that order only) correction
void lambert_getCorrection(double& dval0, double& dval1, double& dval2, double Func0, double Func1, double Func2,
                           double Func3)
{
	double dkA, monebydf, dvm;

	//Because we are operating at least eps away from the bottom of the curve per the fit, the dfunc can never be zero, it only is zero at the bottom,and we are operating above it
	monebydf = -1.0 / Func1;
	dval0 = Func0 * monebydf;

	//NOTE: no time penalty for ivLamParam_orderCorrection checks, since it is a compile time constant
	dvm = dval0 * monebydf;
	dkA = dval0 * dvm;
	dval1 = 0.5 * dkA * Func2;

	dval2 = dkA * dval0 * Func3 / 6.0 + dval1 * Func2 * dvm;
}

//calculate the exact value of k of the minimum times for i-rev solutions.
double lambert_getkbi(double k, int N, double tau)
{
	double tau2 = tau * tau;
	double G, dG, ddG, dddG, G11, G02;
	double W0, W1, W2, W3, W4;
	double detk = 1.0;

	int iter = 0;

	while (fabs(detk) > 1.0e-12 && iter < 10)
	{
		iter++;
		double ktau2mtau = k * tau2 - tau;
		double twok2tau2m4ktaup2 = 2.0 * k * k * tau2 - 4.0 * k * tau + 2.0;
		lambert_getD4W(W0, W1, W2, W3, W4, k, N);
		G = 3.0 * ktau2mtau * W0 + twok2tau2m4ktaup2 * W1 - tau2;
		dG = 3.0 * tau2 * W0 + 7.0 * ktau2mtau * W1 + twok2tau2m4ktaup2 * W2;
		ddG = 10.0 * tau2 * W1 + 11.0 * ktau2mtau * W2 + twok2tau2m4ktaup2 * W3;
		dddG = 21.0 * tau2 * W2 + 15.0 * ktau2mtau * W3 + twok2tau2m4ktaup2 * W4;
		G11 = dG * dG;
		G02 = G * ddG;
		detk = -G * (G11 - G02 / 2.0) / (dG * (G11 - G02) + dddG * G * G / 6.0);
		k = k + detk;
	}

	return k;
}

//initial guess generation
int lambert_getinitialK(double& k, double& kmin, double& kmax, int direction, int N, double tau, double tofbyS,
                        int branch)
{
	double sqrt2 = sqrt(2);
	double x, kn, km, ki;
	double F, F0, F1, Fi;
	double pii, W0, Z, alpha;
	if (N == 0) //Hyperbplic and Elliptical zero revolution
	{
		double Tp = sqrt(1.0 - sqrt2 * tau) * (tau + sqrt2) / 3.0; //parabolic TOF
		if (tofbyS <= Tp) //Hyperbplic
		{
			if (direction == 1)
			{
				kn = sqrt2;
				km = 1.0 / tau;
				ki = (kn + km) / 2.0;
				F = tofbyS;
				F0 = Tp;
				F1 = 0.0;
				pii = 1 - ki * tau;
				lambert_getW(W0, ki, N);
				Fi = sqrt(pii) * (tau + pii * W0);
				Z = sqrt2 / 2.0;
				alpha = 0.5;
				double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
				x = temp * temp;
				k = kn + (km - kn) * x;
				kmin = kn;
				kmax = km;
			}
			else
			{
				double p20 = 1.0 - 20.0 * tau;
				double TOF20 = sqrt(p20) * (tau + 0.04940968903 * p20);
				if (tofbyS >= TOF20) // k is in [sqrt2, 20]
				{
					kn = sqrt2;
					km = 20.0;
					ki = 7.609475708248731;
					F = tofbyS;
					F0 = Tp;
					F1 = TOF20;
					pii = 1 - ki * tau;
					Fi = sqrt(pii) * (tau + pii * 0.12479152768);
					Z = 1.0 / 3.0;
					alpha = 1.0;
					x = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
					k = kn + (km - kn) * x;
					kmin = kn;
					kmax = km;
				}
				else // k is in (20, inf)
				{
					double p100 = 1.0 - 100.0 * tau;
					double T1 = sqrt(p100) * (tau + 0.00999209404 * p100);
					double T0 = TOF20;
					double temp = (T1 * (T0 - tofbyS) * 10.0 - T0 * 4.472135954999580 * (T1 - tofbyS)) / tofbyS / (T0 -
						T1);
					k = temp * temp;
					kmin = 20.0;
					kmax = 1.0e10;
				}
			}
		}
		else //Elliptical zero revolution
		{
			double pn1 = 1.0 + 1.0 * tau;
			double TOFn1 = sqrt(pn1) * (tau + 5.712388981 * pn1);
			if (tofbyS <= TOFn1) // k is in [-1, sqrt2)
			{
				double p0 = 1.0;
				double TOF0 = 1.110720734539592 + tau;
				if (tofbyS <= TOF0) // k is in [0, sqrt2)
				{
					kn = 0.0;
					km = sqrt2;
					ki = sqrt2 / 2.0;
					F = tofbyS;
					F0 = TOF0;
					F1 = Tp;
					pii = 1 - ki * tau;
					Fi = sqrt(pii) * (tau + pii * 0.6686397730);
					Z = 0.5;
					alpha = 1.0;
					x = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
					k = kn + (km - kn) * x;
					kmin = kn;
					kmax = km;
				}
				else // k is in [-1, 0)
				{
					kn = 0.0;
					km = -1.0;
					ki = -0.5;
					F = tofbyS;
					F0 = TOF0;
					F1 = TOFn1;
					pii = 1 - ki * tau;
					Fi = sqrt(pii) * (tau + pii * 1.954946607);
					Z = 0.5;
					alpha = 1.0;
					x = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
					k = kn + (km - kn) * x;
					kmin = km;
					kmax = kn;
				}
			}
			else // k is in [-sqrt2, -1)
			{
				double c1, c2, c3, c4, gamma1, gamma2, gamma3;
				double pn138 = 1.0 + 1.38 * tau;
				double TOFn138 = sqrt(pn138) * (tau + 212.087279879 * pn138);
				if (tofbyS <= TOFn138) // k is in [-1.38, -1)
				{
					kn = -1.0;
					km = -sqrt2;
					ki = -1.38;
					c1 = 1.730076800000000e+02;
					c2 = 256.0;
					c3 = 1.0;
					c4 = 1.0;
					alpha = 16.0;
					F1 = 1.0 / TOFn1;
					Fi = 1.0 / TOFn138;
					F = 1.0 / tofbyS;
					gamma1 = Fi * (F - F1);
					gamma2 = F * (F1 - Fi);
					gamma3 = F1 * (F - Fi);
					k = -c4 * pow(((gamma1 * c1 - c3 * gamma3) * c2 + c3 * c1 * gamma2)
					              / (gamma3 * c1 - c3 * gamma1 - gamma2 * c2), 1.0 / alpha);
					kmin = km;
					kmax = kn;
				}
				else // k is in [-sqrt2, -1.38)
				{
					kn = -1.38;
					km = -sqrt2;
					ki = -1.41;
					c1 = 1.820725082227725;
					c2 = 3.759624518075655;
					c3 = 0.009786288064068;
					c4 = 1.406527249683143;
					alpha = 243.0;
					double pn141 = 1.0 + 1.41 * tau;
					double TOFn141 = sqrt(pn141) * (tau + 4839.684497246 * pn141);
					F1 = 1.0 / TOFn138;
					Fi = 1.0 / TOFn141;
					F = 1.0 / tofbyS;
					gamma1 = Fi * (F - F1);
					gamma2 = F * (F1 - Fi);
					gamma3 = F1 * (F - Fi);
					k = -c4 * pow(((gamma1 * c1 - c3 * gamma3) * c2 + c3 * c1 * gamma2)
					              / (gamma3 * c1 - c3 * gamma1 - gamma2 * c2), 1.0 / alpha);
					kmin = km;
					kmax = kn;
				}
			}
		}
	}
	else //multiple revolution
	{
		double abstau, sgntau, Ebi0, v1, v2, Ebi, kbi, pbi, Tmin;

		//guess the value of k of the minimum times for i - rev solutions.
		abstau = fabs(tau);
		sgntau = 1.0;
		if (tau < 0.0)
			sgntau = -1.0;
		Ebi0 = 3.1415926535897932384626433832795;
		if (N <= 20)
			Ebi0 = ivLam_Ebi0vec[N - 1];
		v2 = Ebi0;
		v1 = 8.0 * abstau / v2 / (sqrt2 - 2.0 * abstau);
		Ebi = v2 * (1.0 - sgntau) + v2 * sgntau * pow(1.0 / (1.0 + v1), 0.25);
		kbi = sqrt2 * fabs(cos(Ebi / 2.0));
		if (3.1415926535897932384626433832795 - Ebi < 0.0)
			kbi = -kbi;

		//calculate the exact value of k of the minimum times for i-rev solutions.
		kbi = lambert_getkbi(kbi, N, tau);
		lambert_getW(W0, kbi, N);
		pbi = 1 - kbi * tau;
		Tmin = sqrt(pbi) * (tau + pbi * W0); //the minimum times for i-rev solutions
		k = kbi;
		kmin = -sqrt2;
		kmax = sqrt2;

		if (tofbyS < Tmin)
			return 0;

		if (kbi >= 0.0)
		{
			if (branch == 0) // k is in [-sqrt2, kbi)
			{
				double p0 = 1.0;
				double TOF0 = 1.110720734539592 + 2.221441469079183 * N + tau;
				if (tofbyS <= TOF0) // k is in [0, kbi)
				{
					kn = 0.0;
					km = kbi;
					ki = kbi / 2.0;
					F = tofbyS;
					F0 = TOF0;
					F1 = Tmin;
					pii = 1 - ki * tau;
					lambert_getW(W0, ki, N);
					Fi = sqrt(pii) * (tau + pii * W0);
					Z = 0.435275281648062;
					alpha = 1.2;
					x = pow(Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F)),
					        0.833333333333333);
					k = kn + (km - kn) * x;
					kmin = kn;
					kmax = km;
				}
				else
				{
					double pn1 = 1.0 + 1.0 * tau;
					double TOFn1 = sqrt(pn1) * (tau + (5.712388981 + 6.283185307179586 * N) * pn1);
					if (tofbyS <= TOFn1) // k is in [-1.0, 0)
					{
						kn = -1.0;
						km = 0.0;
						ki = -0.5;
						F = tofbyS;
						F0 = TOFn1;
						F1 = TOF0;
						pii = 1 - ki * tau;
						Fi = sqrt(pii) * (tau + pii * (1.954946607 + 2.71408094 * N));
						Z = 0.5;
						alpha = 1.0;
						x = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
						k = kn + (km - kn) * x;
						kmin = kn;
						kmax = km;
					}
					else // k is in [-sqrt2, -1.0)
					{
						kn = -1.0;
						km = -sqrt2;
						ki = -(1.0 + 2.0 * sqrt2) / 3.0;
						F = 1.0 / tofbyS;
						F0 = 1.0 / TOFn1;
						F1 = 0.0;
						pii = 1 - ki * tau;
						Fi = 1.0 / (sqrt(pii) * (tau + pii * (27.25239909 + 27.75304668 * N)));
						Z = 0.444444444444444;
						alpha = 2.0;
						double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
						x = sqrt(temp);
						k = kn + (km - kn) * x;
						kmin = km;
						kmax = kn;
					}
				}
			}
			else // k is in [kbi, sqrt2)
			{
				if (kbi >= 1.0) // k is in [kbi, sqrt2)
				{
					kn = kbi;
					km = sqrt2;
					ki = (kbi + sqrt2) / 2.0;
					F = 1.0 / tofbyS;
					F0 = 1.0 / Tmin;
					F1 = 0.0;
					pii = 1 - ki * tau;
					lambert_getW(W0, ki, N);
					Fi = 1.0 / (sqrt(pii) * (tau + pii * W0));
					Z = 0.25;
					alpha = 2.0;
					double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
					x = sqrt(temp);
					k = kn + (km - kn) * x;
					kmin = kn;
					kmax = km;
				}
				else
				{
					double p1 = 1.0 - 1.0 * tau;
					double TOF1 = sqrt(p1) * (tau + (0.57079632 + 6.283185307179586 * N) * p1);
					if (tofbyS >= TOF1) // k is in [1.0, sqrt2)
					{
						kn = 1.0;
						km = sqrt2;
						ki = (1.0 + 2.0 * sqrt2) / 3.0;
						F = 1.0 / tofbyS;
						F0 = 1.0 / TOF1;
						F1 = 0.0;
						pii = 1 - ki * tau;
						Fi = 1.0 / (sqrt(pii) * (tau + pii * (0.50064759 + 27.75304668 * N)));
						Z = 0.444444444444444;
						alpha = 2.0;
						double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
						x = sqrt(temp);
						k = kn + (km - kn) * x;
						kmin = kn;
						kmax = km;
					}
					else // k is in [kbi, 1.0)
					{
						kn = kbi;
						km = 1;
						ki = (1.0 + kbi) / 2.0;
						F = tofbyS;
						F0 = Tmin;
						F1 = TOF1;
						pii = 1 - ki * tau;
						lambert_getW(W0, ki, N);
						Fi = sqrt(pii) * (tau + pii * W0);
						Z = 0.25;
						alpha = 2.0;
						double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
						x = sqrt(temp);
						k = kn + (km - kn) * x;
						kmin = kn;
						kmax = km;
					}
				}
			}
		}
		else
		{
			if (branch == 0) // k is in [-sqrt2, kbi)
			{
				if (kbi <= -1.0) // k is in [-sqrt2, kbi)
				{
					kn = kbi;
					km = -sqrt2;
					ki = (kbi - sqrt2) / 2.0;
					F = 1.0 / tofbyS;
					F0 = 1.0 / Tmin;
					F1 = 0.0;
					pii = 1 - ki * tau;
					lambert_getW(W0, ki, N);
					Fi = 1.0 / (sqrt(pii) * (tau + pii * W0));
					Z = 0.25;
					alpha = 2.0;
					double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
					x = sqrt(temp);
					k = kn + (km - kn) * x;
					kmin = km;
					kmax = kn;
				}
				else
				{
					double pn1 = 1.0 + 1.0 * tau;
					double TOFn1 = sqrt(pn1) * (tau + (5.712388981 + 6.283185307179586 * N) * pn1);
					{
						if (tofbyS >= TOFn1) // k is in [-sqrt2, -1.0]
						{
							kn = -1.0;
							km = -sqrt2;
							ki = -(1.0 + 2.0 * sqrt2) / 3.0;
							F = 1.0 / tofbyS;
							F0 = 1.0 / TOFn1;
							F1 = 0.0;
							pii = 1 - ki * tau;
							Fi = 1.0 / (sqrt(pii) * (tau + pii * (27.25239909 + 27.75304668 * N)));
							Z = 0.444444444444444;
							alpha = 2.0;
							double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
							x = sqrt(temp);
							k = kn + (km - kn) * x;
							kmin = km;
							kmax = kn;
						}
						else // k is in [-1.0, kbi]
						{
							kn = kbi;
							km = -1.0;
							ki = (-1.0 + kbi) / 2.0;
							F = tofbyS;
							F0 = Tmin;
							F1 = TOFn1;
							pii = 1 - ki * tau;
							lambert_getW(W0, ki, N);
							Fi = sqrt(pii) * (tau + pii * W0);
							Z = 0.25;
							alpha = 2.0;
							double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
							x = sqrt(temp);
							k = kn + (km - kn) * x;
							kmin = km;
							kmax = kn;
						}
					}
				}
			}
			else // k is in [kbi, sqrt2]
			{
				double p0 = 1.0;
				double TOF0 = 1.110720734539592 + 2.221441469079183 * N + tau;
				if (tofbyS <= TOF0) // k is in [kbi, 0.0]
				{
					kn = kbi;
					km = 0.0;
					ki = kbi / 2.0;
					F = tofbyS;
					F0 = Tmin;
					F1 = TOF0;
					pii = 1 - ki * tau;
					lambert_getW(W0, ki, N);
					Fi = sqrt(pii) * (tau + pii * W0);
					Z = 0.25;
					alpha = 2;
					double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
					x = sqrt(temp);
					k = kn + (km - kn) * x;
					kmin = kn;
					kmax = km;
				}
				else // k is in (0.0, sqrt2]
				{
					double p1 = 1.0 - tau;
					double TOF1 = sqrt(p1) * (tau + (0.57079632 + 6.283185307179586 * N) * p1);
					if (tofbyS <= TOF1) // k is in (0.0, 1.0]
					{
						kn = 0.0;
						km = 1.0;
						ki = 0.5;
						F = tofbyS;
						F0 = TOF0;
						F1 = TOF1;
						pii = 1 - ki * tau;
						Fi = sqrt(pii) * (tau + pii * (0.75913433 + 2.71408094 * N));
						Z = 0.435275281648062;
						alpha = 1.2;
						x = pow(Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F)),
						        0.833333333333333);
						k = kn + (km - kn) * x;
						kmin = kn;
						kmax = km;
					}
					else // k is in (1.0, sqrt2]
					{
						kn = 1.0;
						km = sqrt2;
						ki = (1.0 + 2 * sqrt2) / 3.0;
						F = 1.0 / tofbyS;
						F0 = 1.0 / TOF1;
						F1 = 0.0;
						pii = 1 - ki * tau;
						Fi = 1.0 / (sqrt(pii) * (tau + pii * (0.50064759 + 27.75304668 * N)));
						Z = 0.444444444444444;
						alpha = 2.0;
						double temp = Z * (F0 - F) * (F1 - Fi) / ((Fi - F) * (F1 - F0) * Z + (F0 - Fi) * (F1 - F));
						x = sqrt(temp);
						k = kn + (km - kn) * x;
						kmin = kn;
						kmax = km;
					}
				}
			}
		}
	}
	return 1;
}


/***********************************Izzo组的Lambert程序**********************************/
double x2tof(double x, double s, double c, int lw, int N)
{
	double t = 0.0;
	double am = s / 2.0;
	double a = am / (1.0 - x * x);
	double beta = 0.0, alfa = 0.0;
	double temp;
	if (x < 1.0) //ELLISSE
	{
		temp = (s - c) / 2.0 / a;
		if (temp < 0.0) temp = 0.0;
		beta = 2.0 * asin(sqrt(temp));
		if (lw)
			beta = -beta;
		alfa = 2.0 * acos(x);
	}
	else //IPERBOLE
	{
		alfa = 2.0 * log(x + sqrt(x * x - 1.0)); //acosh(x);
		temp = (s - c) / (-2.0 * a);
		if (temp < 0.0) temp = 0.0;
		temp = sqrt(temp);
		beta = 2.0 * log(temp + sqrt(temp * temp + 1.0)); //asinh(sqrt((s-c)/(-2.0*a)));
		if (lw)
			beta = -beta;
	}
	if (a > 0.0)
		t = a * sqrt(a) * ((alfa - sin(alfa)) - (beta - sin(beta)) + N * 6.283185307179586476925286766559);
	else
		t = -a * sqrt(-a) * ((sinh(alfa) - alfa) - (sinh(beta) - beta));
	return t;
}

void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
             double tf, const double* unith, int& flag, double mu, int way, int N, int branch, int Maxiter, double tol)
{
	//Originally programmed with Matlab source code by:  Dario Izzo (Advanced Concept Team) Date:
	//28/05/2004 Revision:              1 Tested by:             ----------
	//Translated into C++ source code by Jiang FH
	//Inputs:
	//           R1=Position array at departure
	//           R2=Position array at arrival
	//           tf=Transfer time (scalar)
	//           mu=gravitational parameter (scalar, units have to be
	//           consistent with r1,t units), default 132712439940.0 km^3/s^2
	//           lw=0 for counterclockwise transfer, 1 for clockwise transfer, default 0
	//           N=number of revolutions, default 0
	//           branch=0 if the left branch is selected in a problem where N
	//           is not 0 (multirevolution), 1 for the right one, default 0 (好像是lw=0时，branch=1较好)
	//           Maxiter=the maximum iteration, default 60
	//           tol=tolerance, default 1.0e-11
	//			 unith,R1和R2共线时，需提供的轨道法向单位矢量
	//
	//Outputs:
	//           v1=Velocity at departure        (consistent units)
	//           v2=Velocity at arrival
	//           a=semi major axis of the solution
	//           e=eccentricity of the solution
	flag = 1; //-1, negative time;0, not to be converged; //2, R1 and R2 collinear
	int i;
	if (tf <= 0.0)
	{
		//		cout<<"Negative time as input"<<endl;
		flag = -1;
		return;
	}

	//double tol=1.0e-11;  //Increasing the tolerance does not bring any advantage as the
	//precision is usually greater anyway (due to the rectification of the tof
	//graph) except near particular cases such as parabolas in which cases a
	//lower precision allow for usual convergence.


	double r1[3], r2[3];
	for (i = 0; i < 3; i++)
	{
		r1[i] = R1[i];
		r2[i] = R2[i];
	}
	//Non dimensional units
	double R = sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
	double V = sqrt(mu / R);
	double T = R / V;

	//working with non-dimensional radii and time-of-flight
	for (i = 0; i < 3; i++)
	{
		r1[i] /= R;
		r2[i] /= R;
	}
	tf /= T;

	//Evaluation of the relevant geometry parameters in non dimensional units
	double r2mod = sqrt(r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]); //ArrayNorm2(r2, 3);
	double theta = (r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) / r2mod;
	if (theta > 1.0 - 1.0e-14)
		theta = 0.0;
	else if (theta < -1.0 + 1.0e-14)
		theta = 3.1415926535897932384626433832795;
	else
		theta = acos(theta); //the real command is useful when theta is very
	//close to pi and the acos function could return complex numbers
	double crsprod[3];
	crsprod[0] = r1[1] * r2[2] - r1[2] * r2[1];
	crsprod[1] = r1[2] * r2[0] - r1[0] * r2[2];
	crsprod[2] = r1[0] * r2[1] - r1[1] * r2[0];
	if (crsprod[2] < 0.0)
		theta = 6.283185307179586476925286766559 - theta;
	if (way == 1)
		theta = 6.283185307179586476925286766559 - theta;
	int lw = 0;
	if (theta > 3.1415926535897932384626433832795 + 1.0e-14)
		lw = 1;

	double c = sqrt(1.0 + r2mod * r2mod - 2.0 * r2mod * cos(theta)); //non dimensional chord
	double s = (1.0 + r2mod + c) / 2.0; //non dimensional semi-perimeter
	double am = s / 2.0; //minimum energy ellipse semi major axis
	double lambda = sqrt(r2mod) * cos(theta / 2.0) / s; //lambda parameter defined in BATTIN's book


	if (c <= 1.0e-14)
	{
		a = tf / (2.0 * (1.0 + N) * 3.1415926535897932384626433832795);
		a = pow(a * a * mu, 1.0 / 3.0);
		a *= R;
		e = 0.0;
		v1[0] = unith[1] * R1[2] - unith[2] * R1[1];
		v1[1] = unith[2] * R1[0] - unith[0] * R1[2];
		v1[2] = unith[0] * R1[1] - unith[1] * R1[0];
		v2[0] = v1[0];
		v2[1] = v1[1];
		v2[2] = v1[2];
		return;
	}

	//We start finding the log(x+1) value of the solution conic:
	////NO MULTI REV --> (1 SOL)
	double x = 0.0, inn1 = 0.0, inn2 = 0.0, x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
	double xnew = 0.0, ynew = 0.0;
	int iter = 0;
	if (N == 0)
	{
		inn1 = -0.5233; //first guess point
		inn2 = 0.5233; //second guess point
		x1 = log(1.0 + inn1);
		x2 = log(1.0 + inn2);
		y1 = log(x2tof(inn1, s, c, lw, N)) - log(tf);
		y2 = log(x2tof(inn2, s, c, lw, N)) - log(tf);

		//Newton iterations
		double err = 1.0;
		int i = 0;
		while ((err > tol) && y1 != y2) //(fabs(y1-y2)>1.0e-14))，20211017更改测试替换掉y1!=y2，避免卡死无效
		{
			//20211017晚更改，避免卡死
			if (i > 100000)
				break;
			//更改结束

			i++;
			xnew = (x1 * y2 - y1 * x2) / (y2 - y1);
			ynew = log(x2tof(exp(xnew) - 1.0, s, c, lw, N)) - log(tf);
			x1 = x2;
			y1 = y2;
			x2 = xnew;
			y2 = ynew;
			err = fabs(x1 - xnew);
		}
		iter = i;
		x = exp(xnew) - 1.0;
	}

	////MULTI REV --> (2 SOL) SEPARATING RIGHT AND LEFT BRANCH
	else
	{
		if (branch == 0)
		{
			inn1 = -0.5234;
			inn2 = -0.2234;
		}
		else
		{
			inn1 = 0.7234;
			inn2 = 0.5234;
		}
		x1 = tan(inn1 * 3.1415926535897932384626433832795 / 2.0);
		x2 = tan(inn2 * 3.1415926535897932384626433832795 / 2.0);
		y1 = x2tof(inn1, s, c, lw, N) - tf;
		y2 = x2tof(inn2, s, c, lw, N) - tf;
		double err = 1.0;
		int i = 0;

		//Newton Iteration
		while ((err > tol) && (i < Maxiter) && y1 != y2) //(fabs(y1-y2)>1.0e-14))
		{
			i++;
			xnew = (x1 * y2 - y1 * x2) / (y2 - y1);
			ynew = x2tof(atan(xnew) * 2.0 / 3.1415926535897932384626433832795, s, c, lw, N) - tf;
			x1 = x2;
			y1 = y2;
			x2 = xnew;
			y2 = ynew;
			err = fabs(x1 - xnew);
		}
		x = atan(xnew) * 2.0 / 3.1415926535897932384626433832795;
		iter = i;
	}
	if (iter >= Maxiter)
	{
		//		cout<<"Solution does not seem to be converging"<<endl;
		flag = 0;

		//20211017晚更改，求解失败速度增量返回很大的数
		v1[0] = 10000.0;
		v1[1] = 10000.0;
		v1[2] = 10000.0;
		v2[0] = 10000.0;
		v2[1] = 10000.0;
		v2[2] = 10000.0;
		//更改结束

		return;
	}
	//The solution has been evaluated in terms of log(x+1) or tan(x*pi/2), we
	//now need the conic. As for transfer angles near to pi the lagrange
	//coefficient technique goes singular (dg approaches a zero/zero that is
	//numerically bad) we here use a different technique for those cases. When
	//the transfer angle is exactly equal to pi, then the ih unit vector is not
	//determined. The remaining equations, though, are still valid.
	a = am / (1.0 - x * x);
	double beta = 0.0, alfa = 0.0, psi = 0.0, eta = 0.0, eta2 = 0.0; //solution semimajor axis
	double temp;
	//calcolo psi
	if (x < 1.0) //ellisse
	{
		temp = (s - c) / 2.0 / a;
		if (temp < 0.0) temp = 0.0;
		beta = 2.0 * asin(sqrt(temp));
		if (lw)
			beta = -beta;
		alfa = 2.0 * acos(x);
		psi = (alfa - beta) / 2.0;
		eta2 = 2.0 * a * sin(psi) * sin(psi) / s;
		eta = sqrt(eta2);
	}
	else //iperbole
	{
		temp = (c - s) / 2.0 / a;
		if (temp < 0.0) temp = 0.0;
		temp = sqrt(temp);
		beta = 2.0 * log(temp + sqrt(temp * temp + 1.0)); //asinh(sqrt((c-s)/2.0/a));
		if (lw)
			beta = -beta;
		alfa = 2.0 * log(x + sqrt(x * x - 1.0)); //	alfa=2.0*acosh(x);
		psi = (alfa - beta) / 2.0;
		eta2 = -2.0 * a * sinh(psi) * sinh(psi) / s;
		eta = sqrt(eta2);
	}
	//	if(eta2==0) {a=a*R;
	double p = r2mod / am / eta2 * sin(theta / 2) * sin(theta / 2); //parameter of the solution
	double sigma1 = 1.0 / eta / sqrt(am) * (2.0 * lambda * am - (lambda + x * eta));
	e = 1.0 - p / a;
	e = (e >= 0.0) ? e : 0.0;
	e = sqrt(e);

	temp = sqrt(crsprod[0] * crsprod[0] + crsprod[1] * crsprod[1] + crsprod[2] * crsprod[2]);
	//	if(temp<=0.0) {a=a*R;cout<<"Cannot determine velocity vector"<<endl; flag=2;return;}
	//	for(i=0;i<3;i++)
	//		crsprod[i]/=temp;
	if (temp < 1.0e-14)
	{
		//		flag=2;
		for (i = 0; i < 3; i++) crsprod[i] = unith[i];
	}
	else
		for (i = 0; i < 3; i++) crsprod[i] /= temp;


	if (lw)
		for (i = 0; i < 3; i++) crsprod[i] = -crsprod[i];
	double vr1 = sigma1;
	double vt1 = sqrt(p);
	double v1tan[3] = {0.0}, v2tan[3] = {0.0};
	v1tan[0] = crsprod[1] * r1[2] - crsprod[2] * r1[1];
	v1tan[1] = crsprod[2] * r1[0] - crsprod[0] * r1[2];
	v1tan[2] = crsprod[0] * r1[1] - crsprod[1] * r1[0];
	v2tan[0] = crsprod[1] * r2[2] - crsprod[2] * r2[1];
	v2tan[1] = crsprod[2] * r2[0] - crsprod[0] * r2[2];
	v2tan[2] = crsprod[0] * r2[1] - crsprod[1] * r2[0];
	for (i = 0; i < 3; i++) v2tan[i] /= r2mod;

	double vt2 = vt1 / r2mod;
	double vr2 = -vr1 + (vt1 - vt2) / tan(theta / 2.0);
	for (i = 0; i < 3; i++)
	{
		v1[i] = V * (vr1 * r1[i] + vt1 * v1tan[i]);
		v2[i] = V * (vr2 / r2mod * r2[i] + vt2 * v2tan[i]);
	}
	a = a * R;
	return;
}

//多圈lambert遍历求解
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
              const double* rv0, const double* rv1, double t, double GM)
{
	double vt0[6], am, vt1[3], v0[3], v1[3], a, e, unith[3], dtemp0, dtemp1, dvmin = 1.0e10, dv;
	int i, Nmax, n, flag0;
	flag = 0;
	//	unith[0] = rv0[1] * rv0[5] - rv0[2] * rv0[4];
	//	unith[1] = rv0[2] * rv0[3] - rv0[0] * rv0[5];
	//	unith[2] = rv0[0] * rv0[4] - rv0[1] * rv0[3];
	//	dtemp0 = sqrt(unith[0] * unith[0] + unith[1] * unith[1] + unith[2] * unith[2]);
	V_Cross(unith, rv0, &rv0[3]);
	dtemp0 = V_Norm2(unith, 3);
	if (dtemp0 > 0.0) for (i = 0; i < 3; i++) unith[i] /= dtemp0;
	else
	{
		unith[0] = 0.0;
		unith[1] = 0.0;
		unith[2] = 1.0;
	}
	for (i = 0; i < 3; i++) vt0[i] = rv0[i] - rv1[i];
	am = (V_Norm2(rv0, 3) + V_Norm2(rv1, 3) + V_Norm2(vt0, 3)) / 4.0; //最小椭圆转移轨道半长轴
	//	am = (sqrt(rv0[0] * rv0[0] + rv0[1] * rv0[1] + rv0[2] * rv0[2]) + sqrt(rv1[0] * rv1[0] + rv1[1] * rv1[1] + rv1[2] * rv1[2])
	//		+ sqrt(vt0[0] * vt0[0] + vt0[1] * vt0[1] + vt0[2] * vt0[2])) / 4.0;//最小椭圆转移轨道半长轴
	Nmax = (int)floor(t / (6.283185307179586476925286766559 * sqrt(am * am * am / GM))); //最多可能的转移圈数
	for (n = 0; n <= Nmax; n++) //圈数从0到最大可能数穷举，保留特征速度最小的结果
	{
		for (int j = 0; j < 2; j++) //对于多圈问题，有左右分枝两个解
		{
			lambert(v0, v1, a, e, rv0, rv1, t, unith, flag0, GM, 0, n, j);
			//	if (flag0 != 1 || e >= 1.0) continue;//排除双曲和抛物轨道
			for (i = 0; i < 3; i++)
			{
				vt0[i] = v0[i] - rv0[3 + i];
				vt1[i] = rv1[3 + i] - v1[i];
			}
			dtemp0 = sqrt(vt0[0] * vt0[0] + vt0[1] * vt0[1] + vt0[2] * vt0[2]); // V_Norm2(vt0, 3);
			dtemp1 = sqrt(vt1[0] * vt1[0] + vt1[1] * vt1[1] + vt1[2] * vt1[2]); // V_Norm2(vt1, 3);
			dv = dtemp0 + dtemp1;
			if (dv < dvmin)
			{
				flag = 1;
				N = n;
				branch = j;
				dvmin = dv;
				V_Copy(dv0, vt0, 3);
				V_Copy(dv1, vt1, 3);
				//	dv0[0] = vt0[0]; dv0[1] = vt0[1]; dv0[2] = vt0[2];
				//	dv1[0] = vt1[0]; dv1[1] = vt1[1]; dv1[2] = vt1[2];
				Mdv0 = dtemp0;
				Mdv1 = dtemp1;
			}
			if (n == 0) break;
		}
	}
}
