#include "expansion_perturbed_lambert.h"

DACE::vectorDA ExpansionPerturbedLambert::analyticJ2propHill(const DACE::vectorDA& x0, const DACE::DA& tof,
                                                             const double mu, const double rE, const double J2,
                                                             const double cont)
{
    // cartesian to keplerian elements
    const auto kep0 = MultiRevolutionPerturbedLambertSolver::cart2kep(x0, mu);

    // keplerian elements to Hill
    const auto hill0 = MultiRevolutionPerturbedLambertSolver::kep2hill(kep0, mu);

    // osculating to mean
    const auto hill0Mean = MultiRevolutionPerturbedLambertSolver::osculating2meanHill(hill0, mu, J2, rE, cont);

    // hill to kep
    auto kep0Mean = MultiRevolutionPerturbedLambertSolver::hill2kep(hill0Mean, mu);
    kep0Mean[5] = MultiRevolutionPerturbedLambertSolver::true2meanAnomaly(kep0Mean[5], kep0Mean[1]);

    const auto del0Mean = MultiRevolutionPerturbedLambertSolver::kep2delaunay(kep0Mean, mu);
    auto delfMean = MultiRevolutionPerturbedLambertSolver::averagedJ2rhs(del0Mean, mu, J2, rE, cont);
    delfMean = delfMean * tof + del0Mean;

    auto kepfMean = MultiRevolutionPerturbedLambertSolver::delaunay2kep(delfMean, mu);
    kepfMean[5] = MultiRevolutionPerturbedLambertSolver::mean2trueAnomaly(kepfMean[5], kepfMean[1]);
    const auto hillfMean = MultiRevolutionPerturbedLambertSolver::kep2hill(kepfMean, mu);
    const auto hillf = MultiRevolutionPerturbedLambertSolver::mean2osculatingHill(hillfMean, mu, J2, rE, cont);

    auto xxf = MultiRevolutionPerturbedLambertSolver::hill2cart(hillf, mu);

    return xxf;
}

std::tuple<DACE::vectorDA, DACE::vectorDA> ExpansionPerturbedLambert::expansionOfPerturbedLambert(
    const std::vector<double>& rr1, const std::vector<double>& vv1, double tof,
    const MultiRevolutionPerturbedLambertSolver::mrplp_J2_analytic_parameters& params)
{
    using namespace DACE;

    // initialise DA variables --> DA.init(order, num_variables)
    // DA::init(params.order, 7); // (rr1, vv1, tof)

    // scaling
    const auto Lsc = params.rE;
    const auto Vsc = sqrt(params.mu / params.rE);
    const auto Tsc = Lsc / Vsc;
    const auto muSc = params.mu / Lsc / Lsc / Lsc * Tsc * Tsc;
    constexpr auto sclT = 100.0;

    // DA expansion around the initial state and scaling
    vectorDA x0DA(6);
    x0DA[0] = rr1[0] + DA(1);
    x0DA[1] = rr1[1] + DA(2);
    x0DA[2] = rr1[2] + DA(3);
    x0DA[3] = vv1[0] + DA(4);
    x0DA[4] = vv1[1] + DA(5);
    x0DA[5] = vv1[2] + DA(6);
    for (int i = 0; i < 3; i++)
    {
        x0DA[i] = x0DA[i] / Lsc;
        x0DA[i + 3] = x0DA[i + 3] / Vsc;
    }
    const auto tf = (tof + sclT * DA(7)) / Tsc;

    // propagate
    auto xfDA = analyticJ2propHill(x0DA, tf, muSc, params.rE / Lsc, params.J2, 1.0);

    // direct map
    vectorDA mapD(7);
    mapD[0] = DA(1);
    mapD[1] = DA(2);
    mapD[2] = DA(3);
    mapD[3] = (xfDA[0] - xfDA[0].cons()) * Lsc;
    mapD[4] = (xfDA[1] - xfDA[1].cons()) * Lsc;
    mapD[5] = (xfDA[2] - xfDA[2].cons()) * Lsc;
    mapD[6] = DA(7);

    // invert the map and evaluate intial and final state in the inverse map
    auto mapI = mapD.invert();

    // evaluate the intial and final state in the inverse map
    x0DA = x0DA.eval(mapI);
    xfDA = xfDA.eval(mapI);

    // scale back
    for (int i = 0; i < 3; i++)
    {
        x0DA[i] *= Lsc;
        x0DA[i + 3] *= Vsc;
        xfDA[i] *= Lsc;
        xfDA[i + 3] *= Vsc;
    }

    return std::make_tuple(x0DA, xfDA);
}
