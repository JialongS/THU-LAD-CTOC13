#include "mrpllp_j2_analytic.h"

#include "constant.h"

#include <algorithm>
#include <limits>
#include<cminpack.h>
#include <omp.h>

DACE::DA MultiRevolutionPerturbedLambertSolver::true2eccAnomaly(const DACE::DA& theta, const DACE::DA& e)
{
    using namespace DACE;
    DA E = 2.0 * atan2(sqrt(1.0 - e) * sin(theta / 2.0),
                       sqrt(1.0 + e) * cos(theta / 2.0));
    return E; // using Named Return Value Optimization, NRVO, not using std::move
}

DACE::DA MultiRevolutionPerturbedLambertSolver::true2meanAnomaly(const DACE::DA& theta, const DACE::DA& e)
{
    using namespace DACE;
    const auto E = true2eccAnomaly(theta, e);
    DA M = E - e * sin(E);
    return M;
}

DACE::DA MultiRevolutionPerturbedLambertSolver::mean2eccAnomaly(const DACE::DA& M, const DACE::DA& e)
{
    using namespace DACE;
    auto E = M;
    for (int i = 0; i < 20; i++) // Simple iteration
        E = M + e * sin(E);
    return E;
}

DACE::DA MultiRevolutionPerturbedLambertSolver::ecc2trueAnomaly(const DACE::DA& E, const DACE::DA& e)
{
    using namespace DACE;
    const auto nu = 2.0 * atan2(sqrt(1.0 + e) * sin(E / 2.0), sqrt(1.0 - e) * cos(E / 2.0));
    return nu;
}

DACE::DA MultiRevolutionPerturbedLambertSolver::mean2trueAnomaly(const DACE::DA& M, const DACE::DA& e)
{
    using namespace DACE;
    const auto E = mean2eccAnomaly(M, e);
    const auto nu = ecc2trueAnomaly(E, e);
    return nu;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::cart2kep(const DACE::vectorDA& x0, double mu)
{
    using namespace DACE;

    vectorDA kep(6);
    DA RAAN, omega, theta;

    vectorDA rr(3), vv(3);
    for (int iter = 0; iter < 3; iter++)
    {
        rr[iter] = x0[iter];
        vv[iter] = x0[iter + 3];
    }

    const auto r = rr.vnorm();
    const auto v = vv.vnorm();
    const auto hh = cross(rr, vv);

    DA sma = mu / (2.0 * (mu / r - v * v / 2.0));
    const auto h1sqr = hh[0] * hh[0];
    const auto h2sqr = hh[1] * hh[1];

    if ((h1sqr + h2sqr).cons() == 0.0)
        RAAN = 0.0;
    else
    {
        const auto sinOMEGA = hh[0] / sqrt(h1sqr + h2sqr);
        const auto cosOMEGA = -1.0 * hh[1] / sqrt(h1sqr + h2sqr);
        if (cosOMEGA.cons() >= 0.0)
        {
            if (sinOMEGA.cons() >= 0.0)
                RAAN = asin(hh[0] / sqrt(h1sqr + h2sqr));
            else
                RAAN = 2.0 * DPI + asin(hh[0] / sqrt(h1sqr + h2sqr));
        }
        else
        {
            if (sinOMEGA.cons() >= 0.0)
                RAAN = acos(-1.0 * hh[1] / sqrt(h1sqr + h2sqr));
            else
                RAAN = 2.0 * DPI - acos(-1.0 * hh[1] / sqrt(h1sqr + h2sqr));
        }
    }

    const auto ee = 1.0 / mu * cross(vv, hh) - rr / r;
    DA e = ee.vnorm();
    DA i = acos(hh[2] / hh.vnorm());

    if (e.cons() <= 1.0e-8 && i.cons() <= 1e-8)
    {
        e = 0.0;
        omega = atan2(rr[1], rr[0]);
        theta = 0.0;
        kep[0] = sma;
        kep[1] = e;
        kep[2] = i;
        kep[3] = RAAN;
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }
    if (e.cons() <= 1.0e-8 && i.cons() > 1e-8)
    {
        omega = 0.0;
        vectorDA P(3);
        vectorDA Q(3);
        vectorDA W(3);
        P[0] = cos(omega) * cos(RAAN) - sin(omega) * sin(i) * sin(RAAN);
        P[1] = -1.0 * sin(omega) * cos(RAAN) - cos(omega) * cos(i) * sin(RAAN);
        P[2] = sin(RAAN) * sin(i);
        Q[0] = cos(omega) * sin(RAAN) + sin(omega) * cos(i) * cos(RAAN);
        Q[1] = -1.0 * sin(omega) * sin(RAAN) + cos(omega) * cos(i) * cos(RAAN);
        Q[2] = -1.0 * cos(RAAN) * sin(i);
        W[0] = sin(omega) * sin(i);
        W[1] = cos(omega) * sin(i);
        W[2] = cos(i);
        auto rrt = P * rr[0] + Q * rr[1] + W * rr[2];
        theta = atan2(rrt[1], rrt[0]);
        // omega = atan2(rr[1], rr[0]);
        // theta = 0.0;
        kep[0] = sma;
        kep[1] = e;
        kep[2] = i;
        kep[3] = RAAN;
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }

    auto dotRxE = dot(rr, ee);
    auto RxE = rr.vnorm() * ee.vnorm();
    if (fabs((dotRxE).cons()) > fabs((RxE).cons()) && fabs((dotRxE).cons()) - fabs((RxE).cons()) < fabs(
        1.0e-6 * (dotRxE).cons()))
        dotRxE = 1.0e-6 * dotRxE;

    theta = acos(dotRxE / RxE);

    if (dot(rr, vv).cons() < 0.0)
        theta = 2.0 * DPI - theta;

    if (i.cons() <= 1.0e-8 && e.cons() >= 1.0e-8)
    {
        i = 0.0;
        omega = atan2(ee[1], ee[0]);
        kep[0] = sma;
        kep[1] = e;
        kep[2] = i;
        kep[3] = RAAN;
        kep[4] = omega;
        kep[5] = theta;
        return kep;
    }

    auto sino = rr[2] / r / sin(i);
    auto coso = (rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r;
    DA argLat;

    if (coso.cons() >= 0.0)
    {
        if (sino.cons() >= 0.0)
            argLat = asin(rr[2] / r / sin(i));
        else
            argLat = 2.0 * DPI + asin(rr[2] / r / sin(i));
    }
    else
    {
        if (coso.cons() >= 0.0) // can't arrvel here, becaues in this case coso must < 0
            argLat = acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
        else
            argLat = 2.0 * DPI - acos((rr[0] * cos(RAAN) + rr[1] * sin(RAAN)) / r);
    }

    omega = argLat - theta;

    if (omega.cons() < 0.0)
        omega += 2.0 * DPI;

    //omega = atan2(ee[1], ee[0]);
    kep[0] = sma;
    kep[1] = e;
    kep[2] = i;
    kep[3] = RAAN;
    kep[4] = omega;
    kep[5] = theta;
    return kep;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::kep2Hill(const DACE::vectorDA& kep, const double mu)
{
    using namespace DACE;

    vectorDA hill(6);
    const auto p = kep[0] * (1.0 - kep[1] * kep[1]);
    const auto f = kep[5];

    hill[4] = sqrt(mu * p);
    hill[0] = p / (1.0 + kep[1] * cos(f));
    hill[1] = f + kep[4];
    hill[2] = kep[3];
    hill[3] = (hill[4] / p) * kep[1] * sin(f);
    hill[5] = hill[4] * cos(kep[2]);

    return hill;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::hill2cart(const DACE::vectorDA& hill, double mu)
{
    using namespace DACE;

    const auto& r = hill[0];
    const auto& th = hill[1];
    const auto& nu = hill[2];
    const auto& R = hill[3];
    const auto& Th = hill[4];
    const auto ci = hill[5] / hill[4];
    const auto si = sqrt(1.0 - ci * ci);

    vectorDA u(3);
    u[0] = cos(th) * cos(nu) - ci * sin(th) * sin(nu);
    u[1] = cos(th) * sin(nu) + ci * sin(th) * cos(nu);
    u[2] = si * sin(th);

    vectorDA cart(6);
    cart[0] = r * u[0];
    cart[1] = r * u[1];
    cart[2] = r * u[2];
    cart[3] = (R * cos(th) - Th * sin(th) / r) * cos(nu) - (
        R * sin(th) + Th * cos(th) / r
    ) * sin(nu) * ci;
    cart[4] = (R * cos(th) - Th * sin(th) / r) * sin(nu) + (
        R * sin(th) + Th * cos(th) / r
    ) * cos(nu) * ci;
    cart[5] = (R * sin(th) + Th * cos(th) / r) * si;

    return cart;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::osculating2meanHill(const DACE::vectorDA& hillOsc, double mu,
                                                                          double J2, double rE, const DACE::DA& cont)
{
    using namespace DACE;

    const auto& r = hillOsc[0];
    const auto& th = hillOsc[1];
    const auto& nu = hillOsc[2];
    const auto& R = hillOsc[3];
    const auto& Th = hillOsc[4];
    const auto& Nu = hillOsc[5];

    const auto ci = Nu / Th;
    const auto si = sqrt(1.0 - ci * ci);
    const auto cs = (-1.0 + pow(Th, 2) / (mu * r)) * cos(th) + (R * Th * sin(th)) / mu;
    const auto ss = -((R * Th * cos(th)) / mu) + (-1.0 + pow(Th, 2) / (mu * r)) * sin(th);
    const auto e = sqrt(cs * cs + ss * ss);
    const auto eta = sqrt(1.0 - e * e);

    const auto beta = 1.0 / (1.0 + eta);
    auto p = Th * Th / mu;
    auto costrue = 1.0 / e * (p / r - 1.0);
    auto f = acos(costrue);

    if (R.cons() < 0.0)
        f = 2.0 * DPI - f;

    auto M = true2meanAnomaly(f, e);
    auto phi = f - M;

    auto rMean = r + (cont) * (
        (DACE::pow(rE, 2) * beta * J2) / (2.0 * r)
        - (3.0 * DACE::pow(rE, 2) * beta * J2 * pow(si, 2)) / (4.0 * r)
        + (DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * r) / pow(Th, 4)
        - (3.0 * DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * r * pow(si, 2))
        / (2.0 * pow(Th, 4))
        + (DACE::pow(rE, 2) * J2 * mu) / (2.0 * pow(Th, 2))
        - (DACE::pow(rE, 2) * beta * J2 * mu) / (2.0 * pow(Th, 2))
        - (3.0 * DACE::pow(rE, 2) * J2 * mu * pow(si, 2)) / (4.0 * pow(Th, 2))
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * mu * pow(si, 2)) / (4.0 * pow(Th, 2))
        - (DACE::pow(rE, 2) * J2 * mu * pow(si, 2) * cos(2 * th)) / (4.0 * pow(Th, 2))
    );

    auto thMean = th + (cont) * (
        (-3.0 * DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * phi) / pow(Th, 4)
        + (15.0 * DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * phi * pow(si, 2))
        / (4.0 * pow(Th, 4))
        - (5.0 * DACE::pow(rE, 2) * J2 * mu * R) / (2.0 * pow(Th, 3))
        - (DACE::pow(rE, 2) * beta * J2 * mu * R) / (2.0 * pow(Th, 3))
        + (3.0 * DACE::pow(rE, 2) * J2 * mu * R * pow(si, 2)) / pow(Th, 3)
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * mu * R * pow(si, 2)) / (4.0 * pow(Th, 3))
        - (DACE::pow(rE, 2) * beta * J2 * R) / (2.0 * r * Th)
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * R * pow(si, 2)) / (4.0 * r * Th)
        + (
            -(DACE::pow(rE, 2) * J2 * mu * R) / (2.0 * pow(Th, 3))
            + (DACE::pow(rE, 2) * J2 * mu * R * pow(si, 2)) / pow(Th, 3)
        )
        * cos(2.0 * th)
        + (
            -(DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2)) / (4.0 * pow(Th, 4))
            + (5.0 * DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * pow(si, 2)) / (8.0 * pow(Th, 4))
            + (DACE::pow(rE, 2) * J2 * mu) / (r * pow(Th, 2))
            - (3.0 * DACE::pow(rE, 2) * J2 * mu * pow(si, 2)) / (2.0 * r * pow(Th, 2))
        )
        * sin(2.0 * th)
    );

    auto nuMean = nu + (cont) * (
        (3.0 * DACE::pow(rE, 2) * ci * J2 * DACE::pow(mu, 2) * phi) / (2.0 * pow(Th, 4))
        + (3.0 * DACE::pow(rE, 2) * ci * J2 * mu * R) / (2.0 * pow(Th, 3))
        + (DACE::pow(rE, 2) * ci * J2 * mu * R * cos(2.0 * th)) / (2.0 * pow(Th, 3))
        + (
            (DACE::pow(rE, 2) * ci * J2 * DACE::pow(mu, 2)) / (4.0 * pow(Th, 4))
            - (DACE::pow(rE, 2) * ci * J2 * mu) / (r * pow(Th, 2))
        )
        * sin(2.0 * th)
    );

    auto RMean = R + (cont) * (
        -(DACE::pow(rE, 2) * beta * J2 * R) / (2.0 * pow(r, 2))
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * R * pow(si, 2)) / (4.0 * pow(r, 2))
        - (DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * R) / (2.0 * pow(Th, 4))
        + (3.0 * DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * R * pow(si, 2))
        / (4.0 * pow(Th, 4))
        + (DACE::pow(rE, 2) * J2 * mu * pow(si, 2) * sin(2.0 * th))
        / (2.0 * pow(r, 2) * Th)
    );

    auto ThMean = Th + (cont) * (
        (
            (DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * pow(si, 2)) / (4.0 * pow(Th, 3))
            - (DACE::pow(rE, 2) * J2 * mu * pow(si, 2)) / (r * Th)
        )
        * cos(2.0 * th)
        - (DACE::pow(rE, 2) * J2 * mu * R * pow(si, 2) * sin(2.0 * th))
        / (2.0 * pow(Th, 2))
    );

    auto NuMean = Nu + 0.0;

    vectorDA hillMean(6);

    hillMean[0] = rMean;
    hillMean[1] = thMean;
    hillMean[2] = nuMean;
    hillMean[3] = RMean;
    hillMean[4] = ThMean;
    hillMean[5] = NuMean;

    return hillMean;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::mean2osculatingHill(const DACE::vectorDA& hillMean, double mu,
                                                                          double J2, double rE, const DACE::DA& cont)
{
    using namespace DACE;
    const auto& r = hillMean[0];
    const auto& th = hillMean[1];
    const auto& nu = hillMean[2];
    const auto& R = hillMean[3];
    const auto& Th = hillMean[4];
    const auto& Nu = hillMean[5];
    auto ci = Nu / Th;
    auto si = sqrt(1.0 - ci * ci);
    auto cs = (-1.0 + pow(Th, 2) / (mu * r)) * cos(th) + (R * Th * sin(th)) / mu;
    auto ss = -((R * Th * cos(th)) / mu) + (-1.0 + pow(Th, 2) / (mu * r)) * sin(th);
    auto e = sqrt(cs * cs + ss * ss);
    auto eta = sqrt(1.0 - e * e);
    auto beta = 1.0 / (1.0 + eta);
    auto p = Th * Th / mu;
    auto costrue = 1.0 / e * (p / r - 1.0);
    auto f = acos(costrue);

    if (R.cons() < 0.0)
        f = D2PI - f;

    auto M = true2meanAnomaly(f, e);
    auto phi = f - M;

    auto rOsc = r - (cont) * (
        (DACE::pow(rE, 2) * beta * J2) / (2.0 * r)
        - (3.0 * DACE::pow(rE, 2) * beta * J2 * pow(si, 2)) / (4.0 * r)
        + (DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * r) / pow(Th, 4)
        - (3.0 * DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * r * pow(si, 2))
        / (2.0 * pow(Th, 4))
        + (DACE::pow(rE, 2) * J2 * mu) / (2.0 * pow(Th, 2))
        - (DACE::pow(rE, 2) * beta * J2 * mu) / (2.0 * pow(Th, 2))
        - (3.0 * DACE::pow(rE, 2) * J2 * mu * pow(si, 2)) / (4.0 * pow(Th, 2))
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * mu * pow(si, 2)) / (4.0 * pow(Th, 2))
        - (DACE::pow(rE, 2) * J2 * mu * pow(si, 2) * cos(2 * th)) / (4.0 * pow(Th, 2))
    );

    auto thOsc = th - (cont) * (
        (-3.0 * DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * phi) / pow(Th, 4)
        + (15.0 * DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * phi * pow(si, 2))
        / (4.0 * pow(Th, 4))
        - (5.0 * DACE::pow(rE, 2) * J2 * mu * R) / (2.0 * pow(Th, 3))
        - (DACE::pow(rE, 2) * beta * J2 * mu * R) / (2.0 * pow(Th, 3))
        + (3.0 * DACE::pow(rE, 2) * J2 * mu * R * pow(si, 2)) / pow(Th, 3)
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * mu * R * pow(si, 2)) / (4.0 * pow(Th, 3))
        - (DACE::pow(rE, 2) * beta * J2 * R) / (2.0 * r * Th)
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * R * pow(si, 2)) / (4.0 * r * Th)
        + (
            -(DACE::pow(rE, 2) * J2 * mu * R) / (2.0 * pow(Th, 3))
            + (DACE::pow(rE, 2) * J2 * mu * R * pow(si, 2)) / pow(Th, 3)
        )
        * cos(2.0 * th)
        + (
            -(DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2)) / (4.0 * pow(Th, 4))
            + (5.0 * DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * pow(si, 2)) / (8.0 * pow(Th, 4))
            + (DACE::pow(rE, 2) * J2 * mu) / (r * pow(Th, 2))
            - (3.0 * DACE::pow(rE, 2) * J2 * mu * pow(si, 2)) / (2.0 * r * pow(Th, 2))
        )
        * sin(2.0 * th)
    );

    auto nuOsc = nu - (cont) * (
        (3.0 * DACE::pow(rE, 2) * ci * J2 * DACE::pow(mu, 2) * phi) / (2.0 * pow(Th, 4))
        + (3.0 * DACE::pow(rE, 2) * ci * J2 * mu * R) / (2.0 * pow(Th, 3))
        + (DACE::pow(rE, 2) * ci * J2 * mu * R * cos(2.0 * th)) / (2.0 * pow(Th, 3))
        + (
            (DACE::pow(rE, 2) * ci * J2 * DACE::pow(mu, 2)) / (4.0 * pow(Th, 4))
            - (DACE::pow(rE, 2) * ci * J2 * mu) / (r * pow(Th, 2))
        )
        * sin(2.0 * th)
    );

    auto ROsc = R - (cont) * (
        -(DACE::pow(rE, 2) * beta * J2 * R) / (2.0 * pow(r, 2))
        + (3.0 * DACE::pow(rE, 2) * beta * J2 * R * pow(si, 2)) / (4.0 * pow(r, 2))
        - (DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * R) / (2.0 * pow(Th, 4))
        + (3.0 * DACE::pow(rE, 2) * eta * J2 * DACE::pow(mu, 2) * R * pow(si, 2))
        / (4.0 * pow(Th, 4))
        + (DACE::pow(rE, 2) * J2 * mu * pow(si, 2) * sin(2.0 * th))
        / (2.0 * pow(r, 2) * Th)
    );

    auto ThOsc = Th - (cont) * (
        (
            (DACE::pow(rE, 2) * J2 * DACE::pow(mu, 2) * pow(si, 2)) / (4.0 * pow(Th, 3))
            - (DACE::pow(rE, 2) * J2 * mu * pow(si, 2)) / (r * Th)
        )
        * cos(2.0 * th)
        - (DACE::pow(rE, 2) * J2 * mu * R * pow(si, 2) * sin(2.0 * th))
        / (2.0 * pow(Th, 2))
    );

    auto NuOsc = Nu + 0.0;

    vectorDA hillOsc(6);

    hillOsc[0] = rOsc;
    hillOsc[1] = thOsc;
    hillOsc[2] = nuOsc;
    hillOsc[3] = ROsc;
    hillOsc[4] = ThOsc;
    hillOsc[5] = NuOsc;

    return hillOsc;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::hill2kep(const DACE::vectorDA& hill, double mu)
{
    using namespace DACE;

    const auto& r = hill[0];
    const auto& th = hill[1];
    const auto& nu = hill[2];
    const auto& R = hill[3];
    const auto& Th = hill[4];
    const auto& Nu = hill[5];

    const auto i = acos(Nu / Th);
    const auto cs = (-1.0 + pow(Th, 2) / (mu * r)) * cos(th) + (R * Th * sin(th)) / mu;
    const auto ss = -((R * Th * cos(th)) / mu) + (-1.0 + pow(Th, 2) / (mu * r)) * sin(th);
    const auto e = sqrt(cs * cs + ss * ss);
    const auto p = Th * Th / mu;
    const auto costrue = 1.0 / e * (p / r - 1.0);
    auto f = acos(costrue);

    if (R.cons() < 0.0)
        f = D2PI - f;

    const auto a = p / (1.0 - e * e);

    vectorDA kep(6);
    kep[0] = a;
    kep[1] = e;
    kep[2] = i;
    kep[3] = nu;
    kep[4] = th - f;
    kep[5] = f;

    return kep;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::kep2hill(const DACE::vectorDA& kep, double mu)
{
    using namespace DACE;

    const auto p = kep[0] * (1.0 - kep[1] * kep[1]);
    const auto& f = kep[5];

    vectorDA hill(6);
    hill[4] = sqrt(mu * p);
    hill[0] = p / (1.0 + kep[1] * cos(f));
    hill[1] = f + kep[4];
    hill[2] = kep[3];
    hill[3] = (hill[4] / p) * kep[1] * sin(f);
    hill[5] = hill[4] * cos(kep[2]);

    return hill;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::kep2delaunay(const DACE::vectorDA& kep, double mu)
{
    using namespace DACE;

    const auto& a = kep[0];
    const auto& e = kep[1];
    const auto& i = kep[2];
    const auto& RAAN = kep[3];
    const auto& omega = kep[4];
    const auto& M = kep[5];

    vectorDA delaunay(6);
    delaunay[0] = M;
    delaunay[1] = omega;
    delaunay[2] = RAAN;
    delaunay[3] = sqrt(mu * a);
    delaunay[4] = sqrt(1.0 - pow(e, 2)) * delaunay[3];
    delaunay[5] = cos(i) * delaunay[4];

    return delaunay;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::averagedJ2rhs(const DACE::vectorDA& xxx, double mu, double J2,
                                                                    double rE, const DACE::DA& cont)
{
    using namespace DACE;

    // auto l = xxx[0];
    // auto g = xxx[1];
    // auto h = xxx[2];
    const auto& L = xxx[3];
    const auto& G = xxx[4];
    const auto& H = xxx[5];
    const auto eta = G / L;
    const auto ci = H / G;
    const auto si = sin(acos(ci));

    const auto dldt = DACE::pow(mu, 2) / pow(L, 3) + (cont) * (
        (3.0 * J2 * DACE::pow(rE, 2) * DACE::pow(mu, 4)) / (2.0 * pow(L, 7) * pow(eta, 3))
        - (9.0 * J2 * pow(si, 2) * DACE::pow(rE, 2) * DACE::pow(mu, 4))
        / (4.0 * pow(L, 7) * pow(eta, 3))
    );

    const auto dgdt = (cont) * (
        (3.0 * J2 * DACE::pow(rE, 2) * DACE::pow(mu, 4)) / (2.0 * pow(L, 7) * pow(eta, 4))
        - (9.0 * J2 * pow(si, 2) * DACE::pow(rE, 2) * DACE::pow(mu, 4))
        / (4.0 * pow(L, 7) * pow(eta, 4))
        + (3.0 * pow(ci, 2) * J2 * DACE::pow(rE, 2) * DACE::pow(mu, 4))
        / (2.0 * G * pow(L, 6) * pow(eta, 3))
    );

    const auto dhdt = (cont) * (
        -(3.0 * pow(ci, 2) * J2 * DACE::pow(rE, 2) * DACE::pow(mu, 4))
        / (2.0 * H * pow(L, 6) * pow(eta, 3))
    );

    vectorDA ff(6);
    ff[0] = dldt;
    ff[1] = dgdt;
    ff[2] = dhdt;
    ff[3] = 0.0;
    ff[4] = 0.0;
    ff[5] = 0.0;

    return ff;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::delaunay2kep(const DACE::vectorDA& delaunay, double mu)
{
    using namespace DACE;

    const auto& l = delaunay[0];
    const auto& g = delaunay[1];
    const auto& h = delaunay[2];
    const auto& L = delaunay[3];
    const auto& G = delaunay[4];
    const auto& H = delaunay[5];

    vectorDA kep(6);
    kep[0] = pow(L, 2) / mu;
    kep[1] = sqrt(1 - pow((G / L), 2));
    kep[2] = acos(H / G);
    kep[3] = h;
    kep[4] = g;
    kep[5] = l;

    return kep;
}

DACE::vectorDA MultiRevolutionPerturbedLambertSolver::analyticJ2propHill(const DACE::vectorDA& x0, double tof,
                                                                         double mu, double rE, double J2,
                                                                         const DACE::DA& cont)
{
    using namespace DACE;

    const auto kep0 = cart2kep(x0, mu);

    // keplerian elements to Hill
    const auto hill0 = kep2Hill(kep0, mu);

    // osculating to mean
    const auto hill0Mean = osculating2meanHill(hill0, mu, J2, rE, cont);

    // hill to kep
    auto kep0Mean = hill2kep(hill0Mean, mu);
    kep0Mean[5] = true2meanAnomaly(kep0Mean[5], kep0Mean[1]);

    const auto del0Mean = kep2delaunay(kep0Mean, mu);
    auto delfMean = averagedJ2rhs(del0Mean, mu, J2, rE, cont);
    delfMean = delfMean * tof + del0Mean;
    auto kepfMean = delaunay2kep(delfMean, mu);
    kepfMean[5] = mean2trueAnomaly(kepfMean[5], kepfMean[1]);
    const auto hillfMean = kep2hill(kepfMean, mu);
    const auto hillf = mean2osculatingHill(hillfMean, mu, J2, rE, cont);

    auto xxf = hill2cart(hillf, mu);

    return xxf;
}

Eigen::Matrix<double, 6, 6> MultiRevolutionPerturbedLambertSolver::analyticJ2STM(const DACE::vectorDA& x0, double tof,
    double mu, double rE, double J2, const DACE::DA& cont)
{
    Eigen::Matrix<double, 6, 6> STM = Eigen::Matrix<double, 6, 6>::Identity();
    auto xfDA = analyticJ2propHill(x0, tof, mu, rE, J2, cont + DACE::DA(4));
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            STM(i, j) = xfDA[i].deriv(j + 1).cons();
        }
    }

    return STM;
}


MultiRevolutionPerturbedLambertSolver::mrplp_J2_analytic_parameters::mrplp_J2_analytic_parameters(const type_array& rr1,
    const type_array& rr2, double tof, const type_array& vv1g, double mu, double rE, double J2, int order, double tol,
    double cont, double dcontMin, double scl, int itermax)
{
    this->rr1 = rr1;
    this->rr2 = rr2;
    this->tof = tof;
    this->vv1g = vv1g;
    this->mu = mu;
    this->rE = rE;
    this->J2 = J2;
    this->order = order;
    this->tol = tol;
    this->cont = cont;
    this->dcontMin = dcontMin;
    this->scl = scl;
    this->itermax = itermax;
}

MultiRevolutionPerturbedLambertSolver::output_mrplp_J2_analytic::output_mrplp_J2_analytic(const type_array& vv1Sol,
    const type_array& vv2Sol, const type_array& rr2DA, const type_array& res, int iter, bool success)
{
    this->vv1Sol = vv1Sol;
    this->vv2Sol = vv2Sol;
    this->rr2DA = rr2DA;
    this->res = res;
    this->iter = iter;
    this->success = success;
}

Eigen::Matrix<double, 6, 6> MultiRevolutionPerturbedLambertSolver::stateTransitionMatrix(const std::vector<double>& rr1,
    const std::vector<double>& vv1, double tof, const mrplp_J2_analytic_parameters& params)
{
    using namespace DACE;

    Eigen::Matrix<double, 6, 6> STM = Eigen::Matrix<double, 6, 6>::Zero();

    auto Lsc = params.rE;
    auto Vsc = sqrt(params.mu / params.rE);
    auto Tsc = Lsc / Vsc;
    auto muSc = params.mu / Lsc / Lsc / Lsc * Tsc * Tsc;
    double scl2 = 1.0;
    auto cont = params.cont;

    vectorDA x0DA(6);
    x0DA[0] = rr1[0] + DA(1);
    x0DA[1] = rr1[1] + DA(2);
    x0DA[2] = rr1[2] + DA(3);
    x0DA[3] = vv1[0] + DA(4);
    x0DA[4] = vv1[1] + DA(5);
    x0DA[5] = vv1[2] + DA(6);

    for (int i = 0; i < 3; i++)
    {
        x0DA[i] /= Lsc;
        x0DA[i + 3] /= Vsc;
    }

    auto xfDA = analyticJ2propHill(x0DA, tof / Tsc, muSc, params.rE / Lsc, params.J2, cont + scl2 * DA(4));

    for (int i = 0; i < 3; i++)
    {
        xfDA[i] *= Lsc;
        xfDA[i + 3] *= Vsc;
    }

    for (int j = 0; j < 6; j++)
    {
        for (int k = 0; k < 6; k++)
        {
            STM(j, k) = xfDA[j].deriv(k + 1).cons();
        }
    }

    return STM;
}

std::tuple<DACE::vectorDA, double, DACE::AlgebraicVector<double>> MultiRevolutionPerturbedLambertSolver::
defectAnalyticJ2(const mrplp_J2_analytic_parameters& params)
{
    using namespace DACE;

    // initialise DA variables --> DA.init(order, num_variables)
    // DA::init(params.order, 4);

    auto rr1 = params.rr1;
    auto vv1 = params.vv1g;
    auto rr2 = params.rr2;
    auto tof = params.tof;

    // scaling units
    const auto Lsc = params.rE;
    const auto Vsc = sqrt(params.mu / params.rE);
    const auto Tsc = Lsc / Vsc;
    const auto muSc = params.mu / Lsc / Lsc / Lsc * Tsc * Tsc;

    const auto tol = params.tol;
    const auto scl = params.scl / Vsc;
    constexpr double scl2 = 1.0;
    const auto cont = params.cont;

    // apply the scaling
    for (int i = 0; i < 3; i++)
    {
        rr1[i] = rr1[i] / Lsc;
        vv1[i] = vv1[i] / Vsc;
        rr2[i] = rr2[i] / Lsc;
    }
    tof = tof / Tsc;

    // Taylor expansion around the initial velocity vector
    vectorDA x0(6);
    x0[0] = rr1[0];
    x0[1] = rr1[1];
    x0[2] = rr1[2];
    x0[3] = vv1[0] + scl * DA(1);
    x0[4] = vv1[1] + scl * DA(2);
    x0[5] = vv1[2] + scl * DA(3);

    // propagate
    auto xfDA = analyticJ2propHill(x0, tof, muSc, params.rE / Lsc, params.J2, cont + scl2 * DA(4));

    vectorDA mapD(4);
    mapD[0] = xfDA[0] - rr2[0];
    mapD[1] = xfDA[1] - rr2[1];
    mapD[2] = xfDA[2] - rr2[2];
    mapD[3] = DA(4);

    // evaluate the map in the zero perturbation to get the convergence radius
    vectorDA dxDA(4);
    dxDA[0] = 0 * DA(1);
    dxDA[1] = 0 * DA(2);
    dxDA[2] = 0 * DA(3);
    dxDA[3] = DA(4);
    dxDA = mapD.eval(dxDA);

    // convergence radius
    vectorDA cr(3);
    cr[0] = dxDA[0].convRadius(tol);
    cr[1] = dxDA[1].convRadius(tol);
    cr[2] = dxDA[2].convRadius(tol);

    // new dcont
    double J2eps = scl2 * std::min(std::min(cr[0].cons(), cr[1].cons()), std::min(cr[0].cons(), cr[2].cons()));

    // scale back
    for (int i = 0; i < 3; i++)
    {
        mapD[i] *= Lsc;
        xfDA[i] *= Lsc;
        xfDA[i + 3] *= Vsc;
    }

    return std::make_tuple(mapD, J2eps, xfDA.cons());
}

MultiRevolutionPerturbedLambertSolver::type_array MultiRevolutionPerturbedLambertSolver::fsolveFromMap(
    const type_array& dvv1Guess, const DACE::vectorDA& mapD, const double dcont)
{
    auto temp = dvv1Guess;
    temp.push_back(dcont);
    auto temp1 = mapD.eval(temp);
    type_array res(temp1.begin(), temp1.begin() + 3);

    return res;
}

MultiRevolutionPerturbedLambertSolver::AdapterClassForHybrd1::AdapterClassForHybrd1(
    const DACE::vectorDA& mapD_, double dcont_): mapD(mapD_), dcont(dcont_)
{
}

int MultiRevolutionPerturbedLambertSolver::adapterForHybrd1(int n, const double* x, double* fvec,
                                                            int iflag, const double* para)
{
    if (iflag == 0)
        return 0;
    const auto* myParams = reinterpret_cast<const AdapterClassForHybrd1*>(para);
    const auto mapD = myParams->mapD;
    const auto dcont = myParams->dcont;
    const type_array dvv1Guess(x, x + n);
    const auto res = fsolveFromMap(dvv1Guess, mapD, dcont);
    for (int i = 0; i < n; i++)
    {
        fvec[i] = res[i];
    }
    return 1;
}

MultiRevolutionPerturbedLambertSolver::output_mrplp_J2_analytic MultiRevolutionPerturbedLambertSolver::
mrplp_J2_analytic(mrplp_J2_analytic_parameters& params)
{
    using namespace DACE;

    auto itermax = params.itermax;
    auto cont = params.cont;

    double residual = 10;
    double errormax = 1e-3;
    int iter = 0;
    int exit = 0;
    type_array epsilon;
    type_array res;
    AlgebraicVector<double> finalState_t;


    while ((residual > errormax) || (exit < 1 && iter < itermax))
    {
        // update the iteration number
        iter += 1;

        // solution did not converge --> try with different initial guess
        if (iter > itermax)
        {
            if (residual > errormax)
            {
                double nan = std::numeric_limits<double>::quiet_NaN();
                type_array vv1Sol(3, nan);
                type_array vv2Sol(3, nan);
                type_array rr2DA(3, nan);
                bool success = false;
                output_mrplp_J2_analytic output(vv1Sol, vv2Sol, rr2DA, res, iter, success);
                return output;
            }
        }

        // compute the maps and the new dcont
        auto [mapD, dcont, finalState] = defectAnalyticJ2(params);
        finalState_t = finalState;
        if (dcont < 1.0e-4 && cont < 1)
        {
            double nan = std::numeric_limits<double>::quiet_NaN();
            type_array vv1Sol(3, nan);
            type_array vv2Sol(3, nan);
            type_array rr2DA(3, nan);
            bool success = false;
            output_mrplp_J2_analytic output(vv1Sol, vv2Sol, rr2DA, res, iter, success);
            return output;
        }
        cont += dcont;
        if (cont >= 1.0)
        {
            cont = 1.0;
            dcont = 0.0;
            exit = exit + 1;
        }
        epsilon.push_back(cont);
        std::vector<double> temp0(0, 4);
        auto temp1 = mapD.eval(temp0);
        residual = std::sqrt(temp1[0] * temp1[0] + temp1[1] * temp1[1] + temp1[2] * temp1[2]);
        res.push_back(residual);

        // solve using cminpack::hybrd1 and the maps
        AdapterClassForHybrd1* myParams;
        myParams = new AdapterClassForHybrd1(mapD, dcont);

        double dvv1Guess[3] = {0.0, 0.0, 0.0};
        double fvec[3] = {0.0, 0.0, 0.0};
        double wa[static_cast<int>((3 * (3 * 3 + 13)) / 2)];
        int info = hybrd1(adapterForHybrd1, 3, dvv1Guess, fvec, reinterpret_cast<double*>(myParams), wa, 1.5e-8, 0,
                          400);

        delete myParams;

        // update the params
        for (int i = 0; i < 3; i++)
        {
            params.vv1g[i] += params.scl * dvv1Guess[i];
        }
        params.cont = cont;
    }

    // extract the solution
    type_array vv1Sol = params.vv1g;
    type_array vv2Sol(3);
    type_array rr2DA(3);
    for (int i = 0; i < 3; i++)
    {
        rr2DA[i] = finalState_t[i];
        vv2Sol[i] = finalState_t[i + 3];
    }
    bool success = true;
    output_mrplp_J2_analytic output(vv1Sol, vv2Sol, rr2DA, res, iter, success);
    return output;
}
