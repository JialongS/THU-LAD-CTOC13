#ifndef EXPANSION_PERTURBED_LAMBERT_H
#define EXPANSION_PERTURBED_LAMBERT_H

/**
 * expansion_perturbed_lambert.h
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

#include "mrpllp_j2_analytic.h"

#include <tuple>
#include <dace/dace_s.h>
#include <vector>

class ExpansionPerturbedLambert
{
public:
    ExpansionPerturbedLambert() = default;
    ~ExpansionPerturbedLambert() = default;

    static DACE::vectorDA analyticJ2propHill(const DACE::vectorDA& x0, const DACE::DA& tof, double mu, double rE,
                                             double J2, double cont);
    static std::tuple<DACE::vectorDA, DACE::vectorDA> expansionOfPerturbedLambert(
        const std::vector<double>& rr1, const std::vector<double>& vv1, double tof,
        const MultiRevolutionPerturbedLambertSolver::mrplp_J2_analytic_parameters& params);
};

#endif //EXPANSION_PERTURBED_LAMBERT_H
