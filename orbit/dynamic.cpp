# include "dynamics.h"
# include "constant.h"
# include "ODE45.h"
# include "OrbitFun.h"

# include <cmath>
# include <cstring>

DynamicsEarthJ2::DynamicsEarthJ2(const double homotopyCoefficient_): homotopyCoefficient(homotopyCoefficient_)
{
}

int DynamicsEarthJ2::operator()(double t, const double* rv, double* drv, const double* parameters) const
{
    const double r = sqrt(rv[0] * rv[0] + rv[1] * rv[1] + rv[2] * rv[2]);
    const double coe1 = -MU_EARTH / r / r / r * (1. + 3. / 2. * homotopyCoefficient * J2_EARTH * (RADIUS_EARTH / r) * (
        RADIUS_EARTH / r) * (1. - 5.
        * rv[
            2] * rv[2] /
        r / r));
    const double coe2 = -MU_EARTH / r / r / r * (1. + 3. / 2. * homotopyCoefficient * J2_EARTH * (RADIUS_EARTH / r) * (
        RADIUS_EARTH / r) * (3. - 5.
        * rv[
            2] * rv[2] / r / r));

    drv[0] = rv[3];
    drv[1] = rv[4];
    drv[2] = rv[5];
    drv[3] = rv[0] * coe1;
    drv[4] = rv[1] * coe1;
    drv[5] = rv[2] * coe2;

    return 1;
}

int orbitPropagateJ2(const double t0, const double tf, const double* rv0, double* rvf,
                     const double homotopy_coefficient,
                     const double atol,
                     const double rtol)
{
    DynamicsEarthJ2 dynamics(homotopy_coefficient);

    std::memcpy(rvf, rv0, 6 * sizeof(double));
    auto* abstol = new double[6];
    auto* newwork = new double[6 * 10];
    int num;
    for (int i = 0; i < 6; i++) abstol[i] = atol;

    const int flag = ode45(dynamics, rvf, nullptr, t0, tf, 6, num, newwork, abstol, rtol,
                           0, 100000, 0.1);

    delete[] abstol;
    delete[] newwork;

    return flag;
}

void orbitPropagateJ2Approximation(const double t0, const double tf, const double* rv0, double* rvf)
{
    j2rv02rvf_eJ2(rv0, tf - t0, J2_EARTH, rvf);
}

