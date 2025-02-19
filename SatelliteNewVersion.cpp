#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>

#include "NTN_Deployment.h"
#include "Tables.h"
#include "Matlab_plot.h"
#include "LOS_Probability.h"

int main() {


    auto start = std::chrono::high_resolution_clock::now();

    SatelliteLink UEswithSat(2, 4, 10, 600, 10, 90);
    UEswithSat.generateLinks();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    std::vector<Eigen::Vector3d> users;
    for (const auto& link : UEswithSat.links.getLinks()) {
        users.push_back(link.userPosition);
    }


    // Создание объекта для построения графиков
    MatlabPlot plotter;

    // Построение графиков
    plotter.plotTransformedData(users, UEswithSat.links.getLinks()[0].satellitePosition);
    plotter.plotEarth();

    std::string scenario;
    double f;
    bool isLos;

    std::cout << "\nInput scenario (Dense_Urban/Urban/Suburban/Rural): ";
    std::cin >> scenario;
    std::cout << "\nInput frequency(f) in GHz: ";
    std::cin >> f;


    MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
    MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);

    for (const auto& link : UEswithSat.links.getLinks()) {
        double deg = link.elevationAngle;

        int index = int(AngleForLSP(deg)) / 10;
        std::cout << "User " << ": index: " << index << ", angle: " << AngleForLSP(deg) << "\n";
        isLos = CalculateLOSProbability(index, scenario);

        MatrixXd Table = isLos ? TableLOS : TableNLOS;
        VectorXd Parameters{ Table.rows() };
        Parameters = Table.col(index - 1);

    }


    return 0;
}

// Определение данных из таблиц в зависимости от углов для юзеров, округление к ближайшей колонке


////std::cout << "\n" << Table.rows() << " " << Table.cols();
//std::cout << "\n" << Table << "\n";
//std::vector<std::string> Names;
//// Сценарии Dense_Urban/Urban/Suburban 
//if (scenario != "Rural")
//{
//    if (los)
//    {
//        std::vector<std::string> NamesLOS = {
//            // LSP
//            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ", "K_mu:    ", "K_sg:    ",
//            // Correlation coefficients
//            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ASD_K:   ", "ASA_K:   ", "DS_K:    ", "SF_K:    ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_K:   ", "ZSA_K:   ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
//            // Other Parameters
//            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "DS:      ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     " };
//        Names = NamesLOS;
//    }
//
//    else
//    {
//        std::vector<std::string> NamesNLOS = {
//            // LSP
//            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ",
//            // Correlation coefficients
//            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
//            // Other Parameters
//            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "DS:      ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     ", "CL:      " };
//        Names = NamesNLOS;
//    }
//}
//
////Отдельно сценарий Rural
//if (scenario == "Rural")
//{
//    if (los)
//    {
//        std::vector<std::string> NamesLOS = {
//            // LSP
//            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ", "K_mu:    ", "K_sg:    ",
//            // Correlation coefficients
//            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ASD_K:   ", "ASA_K:   ", "DS_K:    ", "SF_K:    ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_K:   ", "ZSA_K:   ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
//            // Other Parameters
//            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     " };
//        Names = NamesLOS;
//    }
//    else
//    {
//        std::vector<std::string> NamesNLOS = {
//            // LSP
//            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ",
//            // Correlation coefficients
//            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
//            // Other Parameters
//            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     ", "" };
//        Names = NamesNLOS;
//    }
//}
//
//
//for (int i = 0; i < Names.size(); ++i)
//{
//    std::cout << "\033[0m" << Names[i] << "\033[32m" << Parameters[i] << "\n";
//}
//
//std::cout << "\033[0m";