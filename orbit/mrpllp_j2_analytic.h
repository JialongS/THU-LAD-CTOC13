#ifndef MRPLLP_J2_ANALYTIC_H
#define MRPLLP_J2_ANALYTIC_H

/**
 * mrpllp_j2_analytic.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-7-30
 *
 * Translated from Multi-Revolution Perturbed Lambert Problem with Python (MRPLPy)
 * by Andrea Bellome, Rome, Italy, 2024.
 * https://github.com/andreabellome/MRPLP_DACE_python
 *
 * This code is based on the paper:
 * Armellin, R., Gondelach, D., & San Juan, J. F. (2018). Multiple revolution perturbed Lambert problem solvers.
 * Journal of Guidance, Control, and Dynamics, 41(9), 2019-2032. https://doi.org/10.2514/1.G003531.
 */

#include <dace/dace_s.h>
#include <vector>
#include <tuple>
#include <eigen3/Eigen/Dense>

class MultiRevolutionPerturbedLambertSolver
{
public:
    MultiRevolutionPerturbedLambertSolver() = default;
    ~MultiRevolutionPerturbedLambertSolver() = default;

    static DACE::DA true2eccAnomaly(const DACE::DA& theta, const DACE::DA& e);
    static DACE::DA true2meanAnomaly(const DACE::DA& theta, const DACE::DA& e);
    static DACE::DA mean2eccAnomaly(const DACE::DA& M, const DACE::DA& e);
    static DACE::DA ecc2trueAnomaly(const DACE::DA& E, const DACE::DA& e);
    static DACE::DA mean2trueAnomaly(const DACE::DA& M, const DACE::DA& e);
    static DACE::vectorDA cart2kep(const DACE::vectorDA& x0, double mu);
    static DACE::vectorDA kep2Hill(const DACE::vectorDA& kep, double mu);
    static DACE::vectorDA hill2cart(const DACE::vectorDA& hill, double mu);
    static DACE::vectorDA osculating2meanHill(const DACE::vectorDA& hillOsc, double mu, double J2, double rE,
                                              const DACE::DA& cont);
    static DACE::vectorDA mean2osculatingHill(const DACE::vectorDA& hillMean, double mu, double J2, double rE,
                                              const DACE::DA& cont);
    static DACE::vectorDA hill2kep(const DACE::vectorDA& hill, double mu);
    static DACE::vectorDA kep2hill(const DACE::vectorDA& kep, double mu);
    static DACE::vectorDA kep2delaunay(const DACE::vectorDA& kep, double mu);
    static DACE::vectorDA averagedJ2rhs(const DACE::vectorDA& xxx, double mu, double J2, double rE,
                                        const DACE::DA& cont);
    static DACE::vectorDA delaunay2kep(const DACE::vectorDA& delaunay, double mu);
    static DACE::vectorDA analyticJ2propHill(const DACE::vectorDA& x0, double tof, double mu, double rE, double J2,
                                             const DACE::DA& cont);
    static Eigen::Matrix<double, 6, 6> analyticJ2STM(const DACE::vectorDA& x0, double tof, double mu, double rE,
                                                     double J2,
                                                     const DACE::DA& cont);

    using type_array = std::vector<double>;

    class mrplp_J2_analytic_parameters
    {
    public:
        mrplp_J2_analytic_parameters() = delete;
        ~mrplp_J2_analytic_parameters() = default;
        mrplp_J2_analytic_parameters(const type_array& rr1, const type_array& rr2, double tof, const type_array& vv1g,
                                     double mu, double rE, double J2, int order, double tol, double cont,
                                     double dcontMin, double scl, int itermax);

        type_array rr1;
        type_array rr2;
        double tof;
        type_array vv1g;
        double mu;
        double rE;
        double J2;
        int order;
        double tol;
        double cont;
        double dcontMin;
        double scl;
        int itermax;
    };

    class output_mrplp_J2_analytic
    {
    public:
        output_mrplp_J2_analytic() = delete;
        ~output_mrplp_J2_analytic() = default;
        output_mrplp_J2_analytic(const type_array& vv1Sol, const type_array& vv2Sol, const type_array& rr2DA,
                                 const type_array& res, int iter, bool success);

        type_array vv1Sol;
        type_array vv2Sol;
        type_array rr2DA;
        type_array res;
        int iter;
        bool success;
    };

    static Eigen::Matrix<double, 6, 6> stateTransitionMatrix(const std::vector<double>& rr1,
                                                             const std::vector<double>& vv1,
                                                             double tof, const mrplp_J2_analytic_parameters& params);
    static std::tuple<DACE::vectorDA, double, DACE::AlgebraicVector<double>> defectAnalyticJ2(
        const mrplp_J2_analytic_parameters& params);
    static type_array fsolveFromMap(const type_array& dvv1Guess, const DACE::vectorDA& mapD, double dcont);

    struct AdapterClassForHybrd1
    {
        AdapterClassForHybrd1(const DACE::vectorDA& mapD_, double dcont_);

        DACE::vectorDA mapD;
        double dcont{};
    };

    static int adapterForHybrd1(int n, const double* x, double* fvec, int iflag, const double* para = nullptr);
    static output_mrplp_J2_analytic mrplp_J2_analytic(mrplp_J2_analytic_parameters& params);
};

#endif //MRPLLP_J2_ANALYTIC_H
