#include <dynamics.h>
#include <numbers>

#include "optimization.h"
#include "data_ctoc13.h"
#include "InPlaneOptimization_tym.h"
#include "dace/dace_s.h"
#include <iostream>
#include "io.h"
#include <fstream>
#include <VecMat.h>

#include "MultiImpulseTransferOptimization.h"

#include "toJson.h"
#include "one_leg_beam_search.h"

#include <cmath>
#include <execution>
#include "optmain_hpc.h"

/*void sixShipOld()
{
    load_data_ctoc13();

    DACE::DA::init(3, 7);

    int lineID = 1;

    double t_end = 90.0 * 86400.0;

    std::fstream fout;
    fout.open("../output/result.txt", std::ios::out);
    std::string folder = "../output/20240903/";
    std::fstream log;
    log.open("../output/crossPlane.log", std::ios::out);
    for (int launchId = 1; launchId <= 3; launchId++)
    {
        double dv0[3] = {0, 0, 0};
        std::vector<std::array<int, 2>> planeList;
        std::vector<double> timeList;
        std::vector<double> dvEstList;
        std::ifstream fin;
        std::stringstream name;
        name << folder << "plane_" << launchId << ".txt";
        fin.open(name.str());
        //fin.open("../output/p3.txt");
        if (!fin.is_open())
        {
            std::cout << "Can't open plane file" << name.str() << std::endl;
            break;
        }
        // Ship 1
        int satId, planeNum, satNum;
        fin >> satId >> planeNum >> satNum;
        planeNum++;

        timeList.push_back(0); // t0 = 0;
        int n10, n20;
        fin >> n10 >> n20;
        planeList.push_back({n10, n20});
        double dv0_t = 0, dvf_t;

        for (int i = 0; i < planeNum - 1; i++)
        {
            int temp_i;
            int n1, n2;
            fin >> n1 >> n2;
            planeList.push_back({n1, n2});
            double tp0, tpf;
            fin >> tp0 >> tpf >> dvf_t;
            timeList.push_back(tpf);
            dvEstList.push_back(dvf_t - dv0_t);
            dv0_t = dvf_t;
            if (i != planeNum - 2)
            {
                fin >> temp_i >> temp_i;
            }
        }

        double t0 = 0.0, t2, t1;
        double lambda = 0.0;
        double rv0[6], rv1[6], rv2[6];
        double rv0_launch[6];
        double tof;
        std::array<double, 10> trvdv0{}, trvdv1{}, trvdv2{};

        InPlaneTransfer::IPTResult myResult;
        InPlaneTransfer::IPTResult crossResult;
        myResult.clear();
        crossResult.clear();

        std::cout << "Launch " << launchId << " :" << std::endl;
        std::cout << "Ship " << satId << " :" << std::endl;

        //出发轨道面
        int dir = 1;
        auto plt = getPlane(planeList[0][0], planeList[0][1], t0);
        auto pln = getPlane(planeList[1][0], planeList[1][1], t0);
        dir = sin(pln[2] - plt[2]) > 0 ? 1 : -1;
        bool success = InPlaneTransfer(planeList[0][0], planeList[0][1], t0, rv0_launch, &t2, rv2, myResult,
                                       1);

        printCTOC13(fout, lineID, launchId, satId, 0, 0, 0, 0, 0, rv0_launch, dv0);
        lineID++;

        V_Copy(rv0, rv0_launch, 6);

        std::cout << "Leave plane: " << planeList[0][0] << "-" << planeList[0][1] << std::endl;
        myResult.print();
        myResult.save(fout, lineID, launchId, satId);

        for (int iter = 1; iter < planeNum; iter++)
        {
            t0 = t2;
            V_Copy(rv0, rv2, 6);

            t1 = timeList[iter];
            tof = timeList[iter] - timeList[iter - 1] - 90.0 * 60.0 * static_cast<int>(getPlane(
                planeList[iter - 1][0], planeList[iter - 1][1], 0)[3]);
            if ((t1 - t0) < tof)
                t1 = t0 + tof;


            if (t1 > t_end)
                break;

            myResult.clear();
            crossResult.clear();

            dir = 1;
            if (iter != (planeNum - 1))
            {
                plt = getPlane(planeList[iter][0], planeList[iter][1], t0);
                pln = getPlane(planeList[iter + 1][0], planeList[iter + 1][1], t0);
                dir = sin(pln[2] - plt[2]) > 0 ? 1 : -1;
            }
            success = InPlaneTransfer(planeList[iter][0], planeList[iter][1], t1, rv1, &t2, rv2, myResult, 1);
            std::cout << "InPlaneTransfer" << std::endl;

            MultiImpulseTransferOptimization mto;
            MultiImpulseTransfer::MITResult mit_result;
            mto.maxImpulseNum = 5;
            mto.t0 = t0;
            mto.tf = t1;
            for (int i = 0; i < 6; i++)
            {
                mto.rv0[i] = rv0[i];
                mto.rvf[i] = rv1[i];
            }
            mto.optimization(mit_result);
            crossResult.dv_total = mit_result.dvTotal;
            crossResult.timeList = mit_result.timeList;
            crossResult.rvList = mit_result.rvList;
            crossResult.dvList = mit_result.dvList;
            for (int i = 0; i < crossResult.timeList.size(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    crossResult.rvList[i][j + 3] -= crossResult.dvList[i][j];
                }
            }

            std::cout << "interorbit_trans_final!" << std::endl;

            for (int j = 0; j < crossResult.timeList.size(); j++)
            {
                crossResult.eventList.push_back(static_cast<int>(Event::PULSE));
                crossResult.targetList.push_back({0, 0, 0});
            }
            std::cout << "Transfor form plane " << planeList[iter - 1][0] << "-" << planeList[iter - 1][1] <<
                " to plane "
                << planeList[iter][0] << "-" << planeList[iter][1] << std::endl;

            {
                crossResult.print();
                std::cout << "Leave plane: " << planeList[iter][0] << "-" << planeList[iter][1] << std::endl;
                crossResult.save(fout, lineID, launchId, satId);
            }
            myResult.print();
            myResult.save(fout, lineID, launchId, satId);

            if (!success)
                break;
        }
        printCTOC13(fout, lineID, launchId, satId, 3, t2, 0, 0, 0, rv2, dv0);
        lineID++;
        std::cout << "--------------------------------------" << std::endl;


        // Ship 2
        fin >> satId >> planeNum >> satNum;
        planeNum++;
        timeList.clear();
        planeList.clear();
        dvEstList.clear();

        timeList.push_back(0); // t0 = 0;
        fin >> n10 >> n20;
        planeList.push_back({n10, n20});
        dv0_t = 0;

        for (int i = 0; i < planeNum - 1; i++)
        {
            int temp_i;
            int n1, n2;
            fin >> n1 >> n2;
            planeList.push_back({n1, n2});
            double tp0, tpf;
            fin >> tp0 >> tpf >> dvf_t;
            timeList.push_back(tpf);
            dvEstList.push_back(dvf_t - dv0_t);
            dv0_t = dvf_t;
            if (i != planeNum - 2)
            {
                fin >> temp_i >> temp_i;
            }
        }
        fin.close();

        myResult.clear();
        crossResult.clear();

        std::cout << "Launch " << launchId << " :" << std::endl;
        std::cout << "Ship " << satId << " :" << std::endl;

        printCTOC13(fout, lineID, launchId, satId, 0, 0, 0, 0, 0, rv0_launch, dv0);
        lineID++;
        t0 = 0;
        t2 = t0;
        V_Copy(rv2, rv0_launch, 6);

        for (int iter = 1; iter < planeNum; iter++)
        {
            t0 = t2;
            V_Copy(rv0, rv2, 6);

            t1 = timeList[iter];
            tof = timeList[iter] - timeList[iter - 1] - 90.0 * 60.0 * static_cast<int>(getPlane(
                planeList[iter - 1][0], planeList[iter - 1][1], 0)[3]);
            if ((t1 - t0) < tof)
                t1 = t0 + tof;

            if (t1 > t_end)
                break;

            myResult.clear();
            crossResult.clear();
            dir = 1;
            if (iter != (planeNum - 1))
            {
                plt = getPlane(planeList[iter][0], planeList[iter][1], t0);
                pln = getPlane(planeList[iter + 1][0], planeList[iter + 1][1], t0);
                dir = sin(pln[2] - plt[2]) > 0 ? 1 : -1;
            }
            success = InPlaneTransfer(planeList[iter][0], planeList[iter][1], t1, rv1, &t2, rv2, myResult, 1);
            std::cout << "InPlaneTransfer" << std::endl;

            MultiImpulseTransferOptimization mto;
            MultiImpulseTransfer::MITResult mit_result;
            mto.maxImpulseNum = 5;
            mto.t0 = t0;
            mto.tf = t1;
            for (int i = 0; i < 6; i++)
            {
                mto.rv0[i] = rv0[i];
                mto.rvf[i] = rv1[i];
            }
            mto.optimization(mit_result);
            crossResult.dv_total = mit_result.dvTotal;
            crossResult.timeList = mit_result.timeList;
            crossResult.rvList = mit_result.rvList;
            crossResult.dvList = mit_result.dvList;
            for (int i = 0; i < crossResult.timeList.size(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    crossResult.rvList[i][j + 3] -= crossResult.dvList[i][j];
                }
            }

            for (int j = 0; j < crossResult.timeList.size(); j++)
            {
                crossResult.eventList.push_back(static_cast<int>(Event::PULSE));
                crossResult.targetList.push_back({0, 0, 0});
            }
            std::cout << "interorbit_trans_final!" << std::endl;
            std::cout << "Transfor form plane " << planeList[iter - 1][0] << "-" << planeList[iter - 1][1] <<
                " to plane "
                << planeList[iter][0] << "-" << planeList[iter][1] << std::endl;
            {
                crossResult.print();
                std::cout << "Leave plane: " << planeList[iter][0] << "-" << planeList[iter][1] << std::endl;
                crossResult.save(fout, lineID, launchId, satId);
            }
            myResult.print();
            myResult.save(fout, lineID, launchId, satId);

            if (!success)
                break;
        }
        printCTOC13(fout, lineID, launchId, satId, 3, t2, 0, 0, 0, rv2, dv0);
        lineID++;
        std::cout << "--------------------------------------" << std::endl;
    }

    log.close();
    fout.close();
}

*/
// main
int main()
{
 sixShipParallel_hpc();
}

// one circle mulitSat
// int main()
// {
//     load_data_ctoc13();
//     std::vector<double> x(6);
//     MultiSatData md(17, 1, 2);
//     std::array<double, 6> rv0{};
//     double t0 = 3930.0;
//     get_constellation_RV(17, 1, 2, t0, rv0.data());
//     md.rv0 = rv0;
//     md.t0 = t0;
//
//     std::vector<double> grad;
//     double f_best;
//
//     DE_std(multiSat, &md, x, f_best, 6, 0, 5000, 500);
//
//     std::vector<std::array<double, 6>> rvList;
//     std::vector<std::array<int, 3>> targetList;
//     std::vector<double> tList;
//     multiSat(x, grad, &md, rvList, targetList, tList);
//     std::vector<std::array<double, 3>> dvList;
//     std::vector<int> eventList;
//     for (int i = 0; i < tList.size(); i++)
//     {
//         dvList.push_back({0, 0, 0});
//         eventList.push_back(2);
//     }
//
//     std::fstream fout;
//     fout.open("../output/mult.txt", std::ios::out);
//
//     InPlaneTransfer::IPTResult result;
//     result.targetList = targetList;
//     result.timeList = tList;
//     result.dvList = dvList;
//     result.rvList = rvList;
//     result.eventList = eventList;
//     int lineId = 1;
//     result.save(fout, lineId, 1, 1);
//
//     fout.close();
//
//     double coe[6];
//     int flag;
//     rv2coe(flag, coe, rvList[0].data());
//     std::cout << "coe: " << std::endl;
//     for (const auto i : coe)
//     {
//         std::cout << i << " ";
//     }
//     std::cout << std::endl;
//
//     auto plane17 = getPlane(17, 1, 0);
//     auto plane18 = getPlane(18, 1, 0);
//
//     std::cout << "plane17: ";
//     for (const auto i : plane17)
//     {
//         std::cout << i << " ";
//     }
//     std::cout << std::endl;
//
//
//     std::cout << "plane18: ";
//     for (const auto i : plane18)
//     {
//         std::cout << i << " ";
//     }
//     std::cout << std::endl;
//
//     double n = sqrt(MU_EARTH / coe[0] / coe[0] / coe[0]);
//     std::cout << "n: " << n << std::endl;
//
//     double n17 = sqrt(MU_EARTH / plane17[0] / plane17[0] / plane17[0]);
//     std::cout << "n17: " << n17 << std::endl;
//
//     double n18 = sqrt(MU_EARTH / plane18[0] / plane18[0] / plane18[0]);
//     std::cout << "n18: " << n18 << std::endl;
//
//     double k17 = n / n17;
//     double k18 = n / n18;
//
//     std::cout << "n/n17 = " << k17 << std::endl;
//     std::cout << "n/n18 = " << k18 << std::endl;
//
//
//     return 0;
// }
