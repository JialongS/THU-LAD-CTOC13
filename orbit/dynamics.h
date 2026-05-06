#ifndef DYNAMICS_H
#define DYNAMICS_H

/**
 * dynamics.h
 *
 * Copyright (C) Song Jialong
 * Tsinghua University
 * School of Aerospace and Engineering
 * Laboratory of Astrodynamics
 * All rights reserved.
 *
 * This revision: 24-7-28
 */

class DynamicsEarthJ2
{
public:
    DynamicsEarthJ2() = default;
    explicit DynamicsEarthJ2(double homotopyCoefficient_);
    ~DynamicsEarthJ2() = default;

    double homotopyCoefficient{1.0};

    int operator()(double t, const double* rv, double* drv, const double* parameters = nullptr) const;
};

int orbitPropagateJ2(double t0, double tf, const double* rv0, double* rvf, double homotopy_coefficient = 1.0,
                     double atol = 1e-13, double rtol = 1e-13);

void orbitPropagateJ2Approximation(double t0, double tf, const double* rv0, double* rvf);

#endif //DYNAMICS_H
