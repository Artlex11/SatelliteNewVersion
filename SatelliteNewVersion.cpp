#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>

// для сортировки векторов
#include <algorithm>

#include "NTN_Deployment.h"
#include "Tables.h"
#include "Matlab_plot.h"
#include "LOS_Probability.h"
#include "LargeScaleParameters.h"
#include "Pathloss.h"


// Функция для вычисления диаграммы направленности
void calculateDishPattern(Eigen::VectorXd& teta_rad, Eigen::VectorXd& AP_dB, double rDish_WL) {
    for (int i = 0; i < teta_rad.size(); ++i) {
        double teta = teta_rad(i);
        if (teta == 0) {
            AP_dB(i) = 1.0;
        }
        else {
            double bessel_arg = 2 * M_PI * rDish_WL * sin(teta);
            double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
            double res = 4 * std::pow(bessel_val / bessel_arg, 2);
            double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;
            AP_dB(i) = 10 * log10(res) + Gain;

        }
    }
}

int main() {
    MatlabPlot plotter;
    //plotter.plotEarth();

    int mode;
    std::cout << "Enter a number 1 or 2: ";
    std::cin >> mode;

    while (true)
    {
        if (mode == 1)
        {
            while (true) {

                // Scenario options
                std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
                std::vector<std::string> frequencyBands = { "S", "Ka" };

                // Scenario selection
                std::cout << "Select a scenario:\n";
                for (size_t i = 0; i < scenarios.size(); ++i) {
                    std::cout << i + 1 << ". " << scenarios[i] << "\n";
                }
                std::cout << "Enter scenario number (1-" << scenarios.size() << ") or 0 to exit: ";

                int scenarioChoice;
                std::cin >> scenarioChoice;

                if (scenarioChoice == 0) {
                    break; // Exit the loop
                }
                else if (scenarioChoice < 1 || scenarioChoice > scenarios.size()) {
                    std::cout << "Invalid choice. Please try again.\n";
                    continue; // Repeat the loop
                }

                std::string scenario = scenarios[scenarioChoice - 1];

                // Frequency band selection
                std::cout << "Select a frequency band:\n";
                for (size_t i = 0; i < frequencyBands.size(); ++i) {
                    std::cout << i + 1 << ". " << frequencyBands[i] << " band\n";
                }
                std::cout << "Enter frequency band number (1-" << frequencyBands.size() << "): ";

                int frequencyChoice;
                std::cin >> frequencyChoice;

                if (frequencyChoice < 1 || frequencyChoice > frequencyBands.size()) {
                    std::cout << "Invalid choice. Please try again.\n";
                    continue; // Repeat the loop
                }

                double f = (frequencyChoice == 1) ? 2.0 : 30.0; // Example: 2 GHz for S and 30 GHz for Ka


                double elTargetDegrees;
                std::cout << "Enter elTargetDegrees  [grad]  (90 - nadir) : ";
                std::cin >> elTargetDegrees;
                double azTargetDegrees;
                std::cout << "Enter azTargetDegrees  [grad]   : ";
                std::cin >> azTargetDegrees;

                double altitude = 600.0;

                auto start = std::chrono::high_resolution_clock::now();
                SatelliteLink UEswithSat(2, 0, 10, altitude, 10.0, elTargetDegrees, azTargetDegrees);
                UEswithSat.generateLinks();

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "Elapsed time: " << elapsed.count() << " s\n";

                // генератор для шума
                std::random_device rd;  // Источник энтропии
                std::mt19937 gen(rd()); // Генератор случайных чисел (Mersenne Twister)
                std::normal_distribution<> dist(0, 0.72); // Нормальное распределение с mean=0 и stddev=0.72

                // параметры для CINR, CNR, CIR
                double BW_MHz = 30;
                double Ptx_dBW = 4 + 10 * log10(BW_MHz);
                double NF_dB = 7;
                double k_dBWHz = -228.6;
                double T = 290;

                double Pn_dBW = k_dBWHz + 10 * log10(T) + 10 * log10(BW_MHz * 1e6) + NF_dB;
                double Pn_lin = pow(10, (Pn_dBW / 10));

                std::vector<Eigen::Vector3d> users;
                for (const auto& link : UEswithSat.links.getLinks()) {
                    users.push_back(link.userPosition);
                }


                plotter.plotTransformedData(users, UEswithSat.links.getLinks()[0].satellitePosition);
                //plotter.plotEarth();

                // Generate tables for LOS and NLOS
                MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
                MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);

                std::vector<Eigen::Vector3d> rRays = UEswithSat.getVectorsToCellCenters();

                MatrixXd AP_lin(UEswithSat.links.getLinks().size(), rRays.size());
                double rDish_WL = 2e9 / 299792458;
                int countLink = 0;

                for (auto& link : UEswithSat.links.getLinks()) {
                    double deg = link.elevationAngle;
                    int index = int(AngleForLSP(deg)) / 10;
                    bool isLos = CalculateLOSProbability(index, scenario);

                    MatrixXd Table = true ? TableLOS : TableNLOS;
                    VectorXd Parameters = Table.col(index - 1);
                    LSP::initializeParameters(isLos, link, Parameters);


                    for (int i = 0; i < static_cast<int> (rRays.size()); ++i) {

                        double theta = std::acos(rRays[i].dot((link.userPosition - link.satellitePosition).normalized()));


                        if (theta == 0) {
                            AP_lin(countLink, i) = pow(10.0, 1.0 / 20.0);
                        }
                        else {
                            double bessel_arg = 2 * M_PI * rDish_WL * sin(theta);
                            double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
                            double res = 4 * std::pow(bessel_val / bessel_arg, 2);
                            double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;
                            AP_lin(countLink, i) = pow(10.0, (10 * log10(res) + Gain) / 20.0);
                        }
                    }

                    double d = CalculateDistance(EARTH_RADIUS * 1000, altitude * 1000, deg * PI / 180);

                    std::cout << "Distance: " << d;

                    // double std = ChooseSTD(isLos, f, deg, scenario);
                    //std::cout << "\nSTD for SF: " << std;

                    // double SF = GenerateSF(std);
                    //std::cout << "\nSF: " << SF;

                    double FSPL = Calculate_FSPL(d, f);
                    std::cout << "\nFSPL: " << FSPL << std::endl;

                    // double CL = ChooseCL(isLos, f, deg, scenario);
                    //std::cout << "\nCL: " << CL;

                    // double PL_b = CalculateBasisPathLoss(FSPL, SF, CL);
                    // std::cout << "\nPL_b: " << PL_b << std::endl;

                    // double PL_dB = PL_b + dist(gen);

                    double PL_dB = FSPL + dist(gen);

                    // Перевод в линейный масштаб
                    double PL_lin_magnitude = pow(10, -1 * PL_dB / 20);

                    double maxVal = AP_lin.row(countLink).maxCoeff(); // Максимальное значение в строке
                    int maxColIndex;
                    AP_lin.row(countLink).maxCoeff(&maxColIndex); // Индекс максимального значения

                    double attSignal_magn = PL_lin_magnitude * maxVal;

                    double attInterference_magn = sqrt(pow(PL_lin_magnitude * maxVal, 2) - attSignal_magn * attSignal_magn);

                    /*double Prx_lin = pow(10, (Ptx_dBW / 10)) * attSignal_magn * attSignal_magn;
                    double Prx_dBW = 10 * log10(Prx_lin);

                    double Pint_lin = pow(10, (Ptx_dBW / 10)) * attInterference_magn * attSignal_magn;
                    double Pint_dBW = 10 * log10(Pint_lin);*/

                    /*double DL_CNR_dB = Prx_dBW - Pn_dBW;
                    double DL_CIR_dB = Prx_dBW - Pint_dBW;
                    double DL_CINR_dB = Prx_dBW - 10 * log10(Pn_lin + Pint_lin);*/

                    double Prx_dBM = 4 + 10 * log10(BW_MHz);
                    double Prd_dBM = k_dBWHz + 10 * log10(290) + 10 * log10(BW_MHz * 160) + NF_dB;
                    double Prx_in = pow(10, Prx_dBM / 10) * attSignal_magn * attSignal_magn;
                    double Pint_in = pow(10, Prx_dBM / 10) * attInterference_magn * attInterference_magn;
                    double DL_CNR_dB = Prx_dBM - Prd_dBM;
                    double DL_CIR_dB = Prx_dBM - 10 * log10(Pint_in);
                    double DL_CINR_dB = Prx_dBM - 10 * log10(Prx_in + Pint_in);

                    std::cout << "Link #" << countLink + 1 << ": Max AP_lin = " << maxVal << ", Ray #" << maxColIndex + 1 << std::endl;
                    std::cout << "CNR, CIR, CINR: " << DL_CNR_dB << ", " << DL_CIR_dB << ", " << DL_CINR_dB << std::endl;

                    countLink += 1;
                }


                for (int i = 0; i < AP_lin.rows(); ++i) {
                    for (int j = 0; j < AP_lin.cols(); ++j) {
                        std::cout << AP_lin(i, j) << " "; // Выводим элемент
                    }
                    std::cout << std::endl; // Переход на новую строку после вывода строки
                }



                for (size_t i = 0; i < rRays.size(); ++i) {
                    std::cout << "Vector " << i + 1 << ": "
                        << rRays[i].transpose() << std::endl;
                }


                //plotter.plotRayPoints(UEswithSat.links.getLinks()[0].satellitePosition, rRays);


                // Вычисление диаграммы направленности
                //calculateDishPattern(elTargetDegrees, AP_dB, rDish_WL);
            }
        }

        // mode 2 - построенеи графиков для LSP

        else
        {
            // 1) разделить на Los и NLOS
            


            // время работы программы
            /*auto start = std::chrono::high_resolution_clock::now();*/

            std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
            std::vector<std::string> frequencyBands = { "S", "Ka" };

            // Scenario selection
            //std::cout << "Select a scenario:\n";
            //for (size_t i = 0; i < scenarios.size(); ++i) {
            //    std::cout << i + 1 << ". " << scenarios[i] << "\n";
            //}
            //std::cout << "Enter scenario number (1-" << scenarios.size() << ") or 0 to exit: ";

            //int scenarioChoice;
            //std::cin >> scenarioChoice;

            //std::string scenario = scenarios[scenarioChoice - 1];

            // Frequency band selection
            std::cout << "Select a frequency band:\n";
            for (size_t i = 0; i < frequencyBands.size(); ++i) {
                std::cout << i + 1 << ". " << frequencyBands[i] << " band\n";
            }
            std::cout << "Enter frequency band number (1-" << frequencyBands.size() << "): ";

            int frequencyChoice;
            std::cin >> frequencyChoice;

            double f = (frequencyChoice == 1) ? 2.0 : 30.0; // Example: 2 GHz for S and 30 GHz for Ka
            double altitude = 600.0;

            int z = 1;


            //std::vector<double> KDU, KU, KS, KR;
            //std::vector<double> DSDU, DSU, DSS, DSR;
            //std::vector<double> ASDDU, ASDU, ASDS, ASDR;
            //std::vector<double> ASADU, ASAU, ASDS, ASDR;
            //std::vector<double> ZSDDU, ZSDU, ZSDS, ZSDR;
            //std::vector<double> ZSADU, ZSAU, ZSAS, ZSAR;

            //std::vector <std::vector<double>> DenseUrban;
            //std::vector <std::vector<double>> Urban;
            //std::vector <std::vector<double>> Suburban;
            //std::vector <std::vector<double>> Rural;

            // генерация спутников для сбора статистики

            //DenseUrban scenario

            std::vector<double> SFDU;
            std::vector<double> KDU;
            std::vector<double> DSDU;
            std::vector<double> ASDDU;
            std::vector<double> ASADU;
            std::vector<double> ZSDDU;
            std::vector<double> ZSADU;

            int level = 2;
            int Nueincell = 10;
            for (double elevation = 10; elevation <= 90; elevation += 10)
            {
                for (int i = 0; i < z; ++i)
                {
                    // double elTargetDegrees = RandomGenerators::generateUniform(10.0, 90.0);
                    double elTargetDegrees = elevation; 
                    double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);

                    SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
                    UEswithSat.generateLinks();

                    MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[0]);
                    MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[0]);

                    for (auto& link : UEswithSat.links.getLinks()) {
                        double deg = link.elevationAngle;
                        int index = int(AngleForLSP(deg)) / 10;
                        // bool isLos = false;
                        bool isLos = false;// CalculateLOSProbability(index, scenarios[0]);

                        MatrixXd Table = isLos ? TableLOS : TableNLOS;
                        VectorXd Parameters = Table.col(index - 1);
                        LSP::initializeParameters(isLos, link, Parameters);

                        if (isLos)
                        {
                            SFDU.push_back(link.SF_db);
                            KDU.push_back(link.K_db);
                            DSDU.push_back(link.DS_sec);
                            ASDDU.push_back(link.ASD_deg);
                            ASADU.push_back(link.ASA_deg);
                            ZSDDU.push_back(link.ZSD_deg);
                            ZSADU.push_back(link.ZSA_deg);
                        }
                        else
                        {
                            SFDU.push_back(link.SF_db);
                            DSDU.push_back(link.DS_sec);
                            ASDDU.push_back(link.ASD_deg);
                            ASADU.push_back(link.ASA_deg);
                            ZSDDU.push_back(link.ZSD_deg);
                            ZSADU.push_back(link.ZSA_deg);
                        }

                    }

                }
                //std::cout << SF.size() << std::endl;
                //std::cout << K.size() << std::endl;
                //std::cout << DS.size() << std::endl;
                //std::cout << ASD.size() << std::endl;
                //std::cout << ASA.size() << std::endl;
                //std::cout << ZSD.size() << std::endl;
                //std::cout << ZSA.size() << std::endl;

                std::sort(SFDU.data(), SFDU.data() + SFDU.size());
                std::sort(KDU.data(), KDU.data() + KDU.size());
                std::sort(DSDU.data(), DSDU.data() + DSDU.size());
                std::sort(ASDDU.data(), ASDDU.data() + ASDDU.size());
                std::sort(ASADU.data(), ASADU.data() + ASADU.size());
                std::sort(ZSDDU.data(), ZSDDU.data() + ZSDDU.size());
                std::sort(ZSADU.data(), ZSADU.data() + ZSADU.size());

                /*std::cout << "SF: \n";
                for (double x : SF)
                    std::cout << x << std::endl;
                std::cout << "K: \n";
                for (double x : K)
                    std::cout << x << std::endl;
                std::cout << "DS: \n";
                for (double x : DS)
                    std::cout << x << std::endl;
                std::cout << "ASD: \n";
                for (double x : ASD)
                    std::cout << x << std::endl;
                std::cout << "ASA: \n";
                for (double x : ASA)
                    std::cout << x << std::endl;
                std::cout << "ZSD: \n";
                for (double x : ZSD)
                    std::cout << x << std::endl;
                std::cout << "ZSA: \n";
                for (double x : ZSA)
                    std::cout << x << std::endl;*/
            }
            

            // Urban scenario

            std::vector<double> SFU;
            std::vector<double> KU;
            std::vector<double> DSU;
            std::vector<double> ASDU;
            std::vector<double> ASAU;
            std::vector<double> ZSDU;
            std::vector<double> ZSAU;

            for (double elevation = 10; elevation <= 90; elevation += 10)
            {
                for (int i = 0; i < z; ++i)
                {
                    // double elTargetDegrees = RandomGenerators::generateUniform(10.0, 90.0);
                    double elTargetDegrees = elevation; // попрпобовать разные углы спутников 
                    double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);

                    SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
                    UEswithSat.generateLinks();

                    MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[1]);
                    MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[1]);

                    for (int i = 0; i < 1; ++i)
                    {
                        for (auto& link : UEswithSat.links.getLinks()) {
                            double deg = link.elevationAngle;
                            int index = int(AngleForLSP(deg)) / 10;
                            //bool isLos = false;
                            bool isLos = false; // CalculateLOSProbability(index, scenarios[1]);

                            MatrixXd Table = isLos ? TableLOS : TableNLOS;
                            VectorXd Parameters = Table.col(index - 1);
                            LSP::initializeParameters(isLos, link, Parameters);

                            if (isLos)
                            {
                                SFU.push_back(link.SF_db);
                                KU.push_back(link.K_db);
                                DSU.push_back(link.DS_sec);
                                ASDU.push_back(link.ASD_deg);
                                ASAU.push_back(link.ASA_deg);
                                ZSDU.push_back(link.ZSD_deg);
                                ZSAU.push_back(link.ZSA_deg);
                            }
                            else
                            {
                                SFU.push_back(link.SF_db);
                                DSU.push_back(link.DS_sec);
                                ASDU.push_back(link.ASD_deg);
                                ASAU.push_back(link.ASA_deg);
                                ZSDU.push_back(link.ZSD_deg);
                                ZSAU.push_back(link.ZSA_deg);
                            }

                        }

                    }
                    //std::cout << SF.size() << std::endl;
                    //std::cout << K.size() << std::endl;
                    //std::cout << DS.size() << std::endl;
                    //std::cout << ASD.size() << std::endl;
                    //std::cout << ASA.size() << std::endl;
                    //std::cout << ZSD.size() << std::endl;
                    //std::cout << ZSA.size() << std::endl;

                    std::sort(SFU.data(), SFU.data() + SFU.size());
                    std::sort(KU.data(), KU.data() + KU.size());
                    std::sort(DSU.data(), DSU.data() + DSU.size());
                    std::sort(ASDU.data(), ASDU.data() + ASDU.size());
                    std::sort(ASAU.data(), ASAU.data() + ASAU.size());
                    std::sort(ZSDU.data(), ZSDU.data() + ZSDU.size());
                    std::sort(ZSAU.data(), ZSAU.data() + ZSAU.size());

                    /*std::cout << "SF: \n";
                    for (double x : SF)
                        std::cout << x << std::endl;
                    std::cout << "K: \n";
                    for (double x : K)
                        std::cout << x << std::endl;
                    std::cout << "DS: \n";
                    for (double x : DS)
                        std::cout << x << std::endl;
                    std::cout << "ASD: \n";
                    for (double x : ASD)
                        std::cout << x << std::endl;
                    std::cout << "ASA: \n";
                    for (double x : ASA)
                        std::cout << x << std::endl;
                    std::cout << "ZSD: \n";
                    for (double x : ZSD)
                        std::cout << x << std::endl;
                    std::cout << "ZSA: \n";
                    for (double x : ZSA)
                        std::cout << x << std::endl;*/
                }
            }
            

            // Suburban scenario

            std::vector<double> SFS;
            std::vector<double> KS;
            std::vector<double> DSS;
            std::vector<double> ASDS;
            std::vector<double> ASAS;
            std::vector<double> ZSDS;
            std::vector<double> ZSAS;


            for (double elevation = 10; elevation <= 90; elevation += 10)
            {
                for (int i = 0; i < z; ++i)
                {
                    // double elTargetDegrees = RandomGenerators::generateUniform(10.0, 90.0);
                    double elTargetDegrees = elevation;
                    double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);

                    SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
                    UEswithSat.generateLinks();

                    MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[2]);
                    MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[2]);

                    for (int i = 0; i < 1; ++i)
                    {
                        for (auto& link : UEswithSat.links.getLinks()) {
                            double deg = link.elevationAngle;
                            int index = int(AngleForLSP(deg)) / 10;
                            //bool isLos = false;
                            bool isLos = false; // CalculateLOSProbability(index, scenarios[2]);

                            MatrixXd Table = isLos ? TableLOS : TableNLOS;
                            VectorXd Parameters = Table.col(index - 1);
                            LSP::initializeParameters(isLos, link, Parameters);

                            if (isLos)
                            {
                                SFS.push_back(link.SF_db);
                                KS.push_back(link.K_db);
                                DSS.push_back(link.DS_sec);
                                ASDS.push_back(link.ASD_deg);
                                ASAS.push_back(link.ASA_deg);
                                ZSDS.push_back(link.ZSD_deg);
                                ZSAS.push_back(link.ZSA_deg);
                            }
                            else
                            {
                                SFS.push_back(link.SF_db);
                                DSS.push_back(link.DS_sec);
                                ASDS.push_back(link.ASD_deg);
                                ASAS.push_back(link.ASA_deg);
                                ZSDS.push_back(link.ZSD_deg);
                                ZSAS.push_back(link.ZSA_deg);
                            }

                        }

                    }
                    //std::cout << SF.size() << std::endl;
                    //std::cout << K.size() << std::endl;
                    //std::cout << DS.size() << std::endl;
                    //std::cout << ASD.size() << std::endl;
                    //std::cout << ASA.size() << std::endl;
                    //std::cout << ZSD.size() << std::endl;
                    //std::cout << ZSA.size() << std::endl;

                    std::sort(SFS.data(), SFS.data() + SFS.size());
                    std::sort(KS.data(), KS.data() + KS.size());
                    std::sort(DSS.data(), DSS.data() + DSS.size());
                    std::sort(ASDS.data(), ASDS.data() + ASDS.size());
                    std::sort(ASAS.data(), ASAS.data() + ASAS.size());
                    std::sort(ZSDS.data(), ZSDS.data() + ZSDS.size());
                    std::sort(ZSAS.data(), ZSAS.data() + ZSAS.size());

                    /*std::cout << "SF: \n";
                    for (double x : SF)
                        std::cout << x << std::endl;
                    std::cout << "K: \n";
                    for (double x : K)
                        std::cout << x << std::endl;
                    std::cout << "DS: \n";
                    for (double x : DS)
                        std::cout << x << std::endl;
                    std::cout << "ASD: \n";
                    for (double x : ASD)
                        std::cout << x << std::endl;
                    std::cout << "ASA: \n";
                    for (double x : ASA)
                        std::cout << x << std::endl;
                    std::cout << "ZSD: \n";
                    for (double x : ZSD)
                        std::cout << x << std::endl;
                    std::cout << "ZSA: \n";
                    for (double x : ZSA)
                        std::cout << x << std::endl;*/
                }
            }
            

            // Rural scenario
            std::vector<double> SFR;
            std::vector<double> KR;
            std::vector<double> DSR;
            std::vector<double> ASDR;
            std::vector<double> ASAR;
            std::vector<double> ZSDR;
            std::vector<double> ZSAR;

            for (double elevation = 10; elevation <= 90; elevation += 10)
            {
                for (int i = 0; i < z; ++i)
                {
                    //double elTargetDegrees = RandomGenerators::generateUniform(10.0, 90.0);
                    double elTargetDegrees = elevation;
                    double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);

                    SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
                    UEswithSat.generateLinks();

                    MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[3]);
                    MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[3]);

                    for (int i = 0; i < 1; ++i)
                    {
                        for (auto& link : UEswithSat.links.getLinks()) {
                            double deg = link.elevationAngle;
                            int index = int(AngleForLSP(deg)) / 10;
                            // bool isLos = false;
                            bool isLos = false; // CalculateLOSProbability(index, scenarios[3]);

                            MatrixXd Table = isLos ? TableLOS : TableNLOS;
                            VectorXd Parameters = Table.col(index - 1);
                            LSP::initializeParameters(isLos, link, Parameters);

                            if (isLos)
                            {
                                SFR.push_back(link.SF_db);
                                KR.push_back(link.K_db);
                                DSR.push_back(link.DS_sec);
                                ASDR.push_back(link.ASD_deg);
                                ASAR.push_back(link.ASA_deg);
                                ZSDR.push_back(link.ZSD_deg);
                                ZSAR.push_back(link.ZSA_deg);
                            }
                            else
                            {
                                SFR.push_back(link.SF_db);
                                DSR.push_back(link.DS_sec);
                                ASDR.push_back(link.ASD_deg);
                                ASAR.push_back(link.ASA_deg);
                                ZSDR.push_back(link.ZSD_deg);
                                ZSAR.push_back(link.ZSA_deg);
                            }

                        }

                    }
                    //std::cout << SF.size() << std::endl;
                    //std::cout << K.size() << std::endl;
                    //std::cout << DS.size() << std::endl;
                    //std::cout << ASD.size() << std::endl;
                    //std::cout << ASA.size() << std::endl;
                    //std::cout << ZSD.size() << std::endl;
                    //std::cout << ZSA.size() << std::endl;

                    std::sort(SFR.data(), SFR.data() + SFR.size());
                    std::sort(KR.data(), KR.data() + KR.size());
                    std::sort(DSR.data(), DSR.data() + DSR.size());
                    std::sort(ASDR.data(), ASDR.data() + ASDR.size());
                    std::sort(ASAR.data(), ASAR.data() + ASAR.size());
                    std::sort(ZSDR.data(), ZSDR.data() + ZSDR.size());
                    std::sort(ZSAR.data(), ZSAR.data() + ZSAR.size());

                    /*std::cout << "SF: \n";
                    for (double x : SF)
                        std::cout << x << std::endl;
                    std::cout << "K: \n";
                    for (double x : K)
                        std::cout << x << std::endl;
                    std::cout << "DS: \n";
                    for (double x : DS)
                        std::cout << x << std::endl;
                    std::cout << "ASD: \n";
                    for (double x : ASD)
                        std::cout << x << std::endl;
                    std::cout << "ASA: \n";
                    for (double x : ASA)
                        std::cout << x << std::endl;
                    std::cout << "ZSD: \n";
                    for (double x : ZSD)
                        std::cout << x << std::endl;
                    std::cout << "ZSA: \n";
                    for (double x : ZSA)
                        std::cout << x << std::endl;*/
                }
            }
            

                MatlabLSPPlot plotter;
                plotter.plotForAllLSP(SFDU, SFU, SFS, SFR, KDU, KU, KS, KR, DSDU, DSU, DSS, DSR, ASDDU, ASDU, ASDS, ASDR, ASADU, ASAU, ASAS, ASAR, ZSDDU, ZSDU, ZSDS, ZSDR, ZSADU, ZSAU, ZSAS, ZSAR, scenarios, frequencyBands[frequencyChoice - 1]);
                // 
                //auto end = std::chrono::high_resolution_clock::now();
                //std::chrono::duration<double> elapsed = end - start;
                //std::cout << "Elapsed time: " << elapsed.count() << " s\n";

            std::cin;
            }
            
        }

}


//#include <iostream>
//#include <Eigen/Dense>
//#include <vector>
//#include <chrono>
//
//#include "NTN_Deployment.h"
//#include "Tables.h"
//#include "Matlab_plot.h"
//#include "LOS_Probability.h"
//#include "Pathloss.h"
//#include "LargeScaleParameters.h"
//
//
//
//int main() 
//{
//    std::string scenario;
//    double f;
//    bool isLos;
//    double altitude;
//
//    std::cout << "\nInput scenario (DenseUrban/Urban/Suburban/Rural): ";
//    std::cin >> scenario;
//    std::cout << "\nInput frequency(f) in GHz: ";
//    std::cin >> f;
//    std::cout << "\nInput altitude in kilometers: ";
//    std::cin >> altitude;
//    auto start = std::chrono::high_resolution_clock::now();
//
//    SatelliteLink UEswithSat(2, 4, 10, altitude, 10, 30); // вместо чисел вставить переменные, потому что где-то надо их вызывать
//    UEswithSat.generateLinks();
//
//    auto end = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = end - start;
//
//    //std::cout << "Elapsed time: " << elapsed.count() << " s\n";
//
//    std::vector<Eigen::Vector3d> users;
//    for (const auto& link : UEswithSat.links.getLinks()) {
//        users.push_back(link.userPosition);
//    }
//
//
//    // Создание объекта для построения графиков
//    MatlabPlot plotter;
//
//    // Построение графиков
//    plotter.plotTransformedData(users, UEswithSat.links.getLinks()[0].satellitePosition);
//    plotter.plotEarth();
//
//    MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
//    MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);
//
//    for (auto& link : UEswithSat.links.getLinks()) {
//        double deg = link.elevationAngle;
//
//        int index = int(AngleForLSP(deg)) / 10;
//        std::cout << "\nUser " << ": index: " << index << ", angle: " << AngleForLSP(deg) << "\n";
//        isLos = CalculateLOSProbability(index, scenario);
//        link.isLOS = isLos;
//        double d = CalculateDistance(EARTH_RADIUS, altitude, deg * PI / 180);
//
//        std::cout << "Distance: " << d;
//
//        double std = ChooseSTD(isLos, f, deg, scenario);
//        //std::cout << "\nSTD for SF: " << std;
//
//        double SF = GenerateSF(std);
//        //std::cout << "\nSF: " << SF;
//
//        double FSPL = Calculate_FSPL(d, f);
//        //std::cout << "\nFSPL: " << FSPL;
//
//        double CL = ChooseCL(isLos, f, deg, scenario);
//        //std::cout << "\nCL: " << CL;
//
//        double PL_b = CalculateBasisPathLoss(FSPL, SF, CL);
//        std::cout << "\nPL_b: " << PL_b << std::endl;
//
//        setlocale(LC_ALL, "RU");
//
//        /*double PL_g = CalculatePathLossInGasses(d, f);
//        std::cout << "PL_g: " << PL_g  << std::endl;
//        std::cout << "PL = PL_b + PL_g = " << PL_b + PL_g << std::endl;*/
//
//        MatrixXd Table = isLos ? TableLOS : TableNLOS;
//        VectorXd Parameters{ Table.rows() };
//        Parameters = Table.col(index - 1);
//        // расстояние в линке 
//
//        LSP::initializeParameters(isLos, link, Parameters);
//    }
//
//    std::cin >> f;
//
//    return 0;
//}
//
//// Определение данных из таблиц в зависимости от углов для юзеров, округление к ближайшей колонке
//
//
//////std::cout << "\n" << Table.rows() << " " << Table.cols();
////std::cout << "\n" << Table << "\n";
////std::vector<std::string> Names;
////// Сценарии DenseUrban/Urban/Suburban 
////if (scenario != "Rural")
////{
////    if (los)
////    {
////        std::vector<std::string> NamesLOS = {
////            // LSP
////            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ", "K_mu:    ", "K_sg:    ",
////            // Correlation coefficients
////            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ASD_K:   ", "ASA_K:   ", "DS_K:    ", "SF_K:    ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_K:   ", "ZSA_K:   ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
////            // Other Parameters
////            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "DS:      ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     " };
////        Names = NamesLOS;
////    }
////
////    else
////    {
////        std::vector<std::string> NamesNLOS = {
////            // LSP
////            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ",
////            // Correlation coefficients
////            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
////            // Other Parameters
////            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "DS:      ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     ", "CL:      " };
////        Names = NamesNLOS;
////    }
////}
////
//////Отдельно сценарий Rural
////if (scenario == "Rural")
////{
////    if (los)
////    {
////        std::vector<std::string> NamesLOS = {
////            // LSP
////            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ", "K_mu:    ", "K_sg:    ",
////            // Correlation coefficients
////            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ASD_K:   ", "ASA_K:   ", "DS_K:    ", "SF_K:    ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_K:   ", "ZSA_K:   ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
////            // Other Parameters
////            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     " };
////        Names = NamesLOS;
////    }
////    else
////    {
////        std::vector<std::string> NamesNLOS = {
////            // LSP
////            "DS_mu:   ", "DS_sg:   ", "ASD_mu:  ", "ASD_sg:  ", "ASA_mu:  ", "ASA_sg:  ", "ZSA_mu:  ", "ZSA_sg:  ", "ZSD_mu:  ", "ZSD_sg:  ", "SF_mu:   ", "SF_sg:   ",
////            // Correlation coefficients
////            "ASD_DS:  ", "ASA_DS:  ", "ASA_SF:  ", "ASD_SF:  ", "DS_SF:   ", "ASD_ASA: ", "ZSD_SF:  ", "ZSA_SF:  ", "ZSD_DS:  ", "ZSA_DS:  ", "ZSD_ASD: ", "ZSA_ASD: ", "ZSD_ASA: ", "ZSA_ASA: ", "ZSD_ZSA: ",
////            // Other Parameters
////            "r_tau:   ", "XPR_mu:  ", "XPR_sg:  ", "N:       ", "M:       ", "ASD:     ", "ASA:     ", "ZSA:     ", "ksi:     ", "" };
////        Names = NamesNLOS;
////    }
////}
////
////
////for (int i = 0; i < Names.size(); ++i)
////{
////    std::cout << "\033[0m" << Names[i] << "\033[32m" << Parameters[i] << "\n";
////}
////
////std::cout << "\033[0m";