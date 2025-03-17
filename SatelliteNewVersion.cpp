#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>
#include <engine.h>
#include <fstream>

#include "NTN_Deployment.h"
#include "Tables.h"
#include "LOS_Probability.h"
#include "Pathloss.h"
#include "LargeScaleParameters.h"
#include "SmallScaleParameters.h"
#include "LinkBudget.h"
#include "Matlab_plot.h"

void plotCDF(const std::vector<double>& data, const std::string& parameterName, const std::string& scenario, Engine* ep, int figureNumber, int frequencyChoice, std::vector<std::string> frequencyBands)
{
    // Передача данных в MATLAB
    mxArray* mxData = mxCreateDoubleMatrix(1, data.size(), mxREAL);
    std::memcpy(mxGetPr(mxData), data.data(), data.size() * sizeof(double));
    engPutVariable(ep, "data", mxData);

    // Создание новой фигуры для каждого параметра
    std::string command =
        "figure(" + std::to_string(figureNumber) + "); "
        "[f, x] = ecdf(data); "
        "plot(x, f, 'DisplayName', '" + scenario + "'); "
        "title('CDF of " + parameterName + " for " + frequencyBands[frequencyChoice-1] + "-band'); "
        "xlabel('" + parameterName + "'); "
        "ylabel('CDF'); "
        "grid on; "
        "hold on; "  // Включаем режим добавления графиков
        "legend ('show', 'location', 'best')";

    engEvalString(ep, command.c_str());

    // Освобождение памяти
    mxDestroyArray(mxData);
}

// Функция для вычисления диаграммы направленности
void calculateDishPattern(Eigen::VectorXd& teta_rad, Eigen::VectorXd& AP_dB, double rDish_WL) {
    double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;
    for (int i = 0; i < teta_rad.size(); ++i) {
        double teta = teta_rad(i);
        if (teta == 0) {
            AP_dB(i) = 1.0;
        }
        else {
            double bessel_arg = 2 * M_PI * rDish_WL * sin(teta);
            double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
            double res = 4 * std::pow(bessel_val / bessel_arg, 2);
            AP_dB(i) = 10 * log10(res) + Gain;
        }
    }
}

// Функция для поиска индексов совпадающих элементов
std::vector<int> findMatchingIndices(const std::vector<Eigen::Vector2d>& uvSet, const std::vector<Eigen::Vector2d>& uvSetCore) {
    std::vector<int> matchingIndices;
    std::unordered_set<Eigen::Vector2d, Vector2dHash> coreSet(uvSetCore.begin(), uvSetCore.end());
    for (int i = 0; i < uvSet.size(); ++i) {
        if (coreSet.find(uvSet[i]) != coreSet.end()) {
            matchingIndices.push_back(i);
        }
    }
    return matchingIndices;
}

// Функция для маски (фильтра)
std::vector<double> filterVector(const std::vector<double>& inputVector, const std::vector<int>& mask) {
    std::vector<double> filteredVector;
    if (inputVector.size() != mask.size()) {
        std::cerr << "Error: inputVector and mask have different sizes!" << std::endl;
        return filteredVector;
    }
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i] == 1) {
            filteredVector.push_back(inputVector[i]);
        }
    }
    return filteredVector;
}

// Функция для фильтрации строк матрицы
MatrixXd filterMatrixRows(const MatrixXd& inputMatrix, const std::vector<int>& mask) {
    if (mask.size() != inputMatrix.rows()) {
        std::cerr << "Error: Mask size does not match the number of rows in the matrix!" << std::endl;
        return MatrixXd();
    }
    int numRowsToKeep = 0;
    for (int value : mask) {
        if (value == 1) {
            numRowsToKeep++;
        }
    }
    MatrixXd filteredMatrix(numRowsToKeep, inputMatrix.cols());
    int filteredRowIndex = 0;
    for (int i = 0; i < mask.size(); ++i) {
        if (mask[i] == 1) {
            filteredMatrix.row(filteredRowIndex) = inputMatrix.row(i);
            filteredRowIndex++;
        }
    }
    return filteredMatrix;
}

int main() {
    Engine* ep = engOpen("");
    if (!ep) {
        std::cerr << "MATLAB Engine no open" << std::endl;
        return 1;
    }
    engSetVisible(ep, false);

    int gear;
    std::cout << "Chose gear: 0 - read file.txt and create Fig , 1 - Create model channel \n";
    std::cin >> gear;

    MatrixXd TableLOS;
    MatrixXd TableNLOS;

    if (gear) {

        std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
        std::vector<std::string> frequencyBands = { "S", "Ka" };

        int nTiersCore = 2, nTiresWrArnd = 2, nUePerCell;
        double satHeightKm = 600.0, elMinDegrees = 10.0, elTargetDegrees = 30.0, azTargetDegrees = 0.0;

        //std::cout << "Enter number of tiers in the core: ";
        //std::cin >> nTiersCore;
        //std::cout << "Enter number of tiers in the wrap-around: ";
        //std::cin >> nTiresWrArnd;
        std::cout << "Enter number of users per cell: ";
        std::cin >> nUePerCell;
        //std::cout << "Enter satellite height in kilometers: ";
        //std::cin >> satHeightKm;
        //std::cout << "Enter minimum elevation angle in degrees: ";
        //std::cin >> elMinDegrees;
        //std::cout << "Enter target elevation angle in degrees (90 - nadir): ";
        //std::cin >> elTargetDegrees;
        //std::cout << "Enter target azimuth angle in degrees: ";
        //std::cin >> azTargetDegrees;

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
            return 1;
        }
        double f = (frequencyChoice == 1) ? 2.0 : 20.0;

        for (size_t i = 0; i < scenarios.size(); ++i) {
            const auto& scenario = scenarios[i];

            auto start = std::chrono::high_resolution_clock::now();

            SatelliteLink UEswithSat(nTiersCore, nTiresWrArnd, nUePerCell, satHeightKm, elMinDegrees, elTargetDegrees, azTargetDegrees);
            UEswithSat.generateLinks();

            //std::vector<Eigen::Vector3d> users;
            //for (const auto& link : UEswithSat.links.getLinks()) {
            //    users.push_back(link.userPosition);
            //}
            //MatlabPlot::plotEarth(ep);
            //MatlabPlot::plotTransformedData( ep ,users, UEswithSat.links.getLinks()[0].satellitePosition);

            std::vector<Eigen::Vector3d> rRays = UEswithSat.getVectorsToCellCenters();
            MatrixXd arrPatt_magn(UEswithSat.links.getLinks().size(), rRays.size());
            std::vector<double> arrPatt;
            std::vector<double> PL_lin_magn;

            double rDish_WL = 2e9 / 299792458;
            int countLink = 0;

            std::vector<int> indRayCorelist = findMatchingIndices(UEswithSat.uvSet, UEswithSat.uvSetCore);
            std::vector<int> RayUElist;
            std::vector<int> coreUElist;

            Eigen::Matrix3d rotMatrix = UEswithSat.rotMatrix;
            std::vector<double> SF_db_values, K_db_values, DS_sec_values, ASA_deg_values, ASD_deg_values, ZSA_deg_values, ZSD_deg_values, ClusterDelay_sec_values, ClusterPower_lin_values;
            std::vector<double> AOD_values, AOA_values, ZOD_values, ZOA_values;

            TableLOS = GenerateMatrix(true, f, scenario);
            TableNLOS = GenerateMatrix(false, f, scenario);

            for (auto& link : UEswithSat.links.getLinks()) {

                int index = int(AngleForLSP(link.elevationAngle)) / 10;
                link.isLos = CalculateLOSProbability(index, scenario);

                MatrixXd Table = link.isLos ? TableLOS : TableNLOS;
                //VectorXd Parameters = Table.col(index - 1);
                VectorXd Parameters = InterpolatedParameters(Table, link.elevationAngle);

                LSP::initializeParameters(link.isLos, link, Parameters);

                SSP::setParameters(Parameters);
                SSP::calculateLosAngles(link, UEswithSat.p1, UEswithSat.p2);
                SSP::generateClusterDelays(link);
                std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);

                //std::cout << "Power before:\n";
                //for (double val : clusterPower.first) {
                //    std::cout << val << " ";
                //}
                //std::cout << "\n";

                SSP::generateArrivalAndDepartureAngles(link, clusterPower.second);
                SSP::generateXPR(link);
                SSP::sortRelativeToFirstVector(clusterPower.first, link.clusterDelays, link.AOA_n_m, link.AOD_n_m, link.ZOA_n_m, link.ZOD_n_m, link.XPR_n_m);

                //std::cout << "Power after:\n";
                //for (double val : clusterPower.first) {
                //    std::cout << val << " ";
                //}
                //std::cout << "\n";

                double d = CalculateDistance(EARTH_RADIUS, satHeightKm, link.elevationAngle * PI / 180);
                double std = ChooseSTD(link.isLos, f, link.elevationAngle, scenario);
                double SF = GenerateSF(std);
                double FSPL = Calculate_FSPL(d * 1e3, f);
                double CL = ChooseCL(link.isLos, f, link.elevationAngle, scenario);
                double PL_dB = CalculateBasisPathLoss(FSPL, SF, CL);

                std::random_device rd;
                std::mt19937 gen(rd());
                std::normal_distribution<> dist(0, 0.72);
                PL_dB = PL_dB + dist(gen);

                double PL_lin_magnitude = pow(10, -1 * PL_dB / 20);

                for (int i = 0; i < static_cast<int>(rRays.size()); ++i) {
                    Eigen::Vector3d rotatedVector = rotMatrix * (link.userPosition - link.satellitePosition).normalized();
                    double theta = std::acos(rRays[i].dot(rotatedVector));
                    double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;

                    if (theta == 0) {
                        arrPatt_magn(countLink, i) = pow(10.0, (1.0 + Gain) / 20.0);
                    }
                    else {
                        double bessel_arg = 2 * M_PI * rDish_WL * sin(theta);
                        double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
                        double res = 4 * std::pow(bessel_val / bessel_arg, 2);
                        arrPatt_magn(countLink, i) = pow(10.0, (10 * log10(res) + Gain) / 20.0);
                    }
                }

                PL_lin_magn.push_back(PL_lin_magnitude);
                arrPatt.push_back(arrPatt_magn.row(countLink).maxCoeff(&RayUElist.emplace_back()));

                ClusterDelay_sec_values.insert(ClusterDelay_sec_values.end(), link.clusterScaledDelays.begin(), link.clusterScaledDelays.end());
                //ClusterPower_lin_values.insert(ClusterPower_lin_values.end(), clusterPower.first.begin(), clusterPower.first.end());

                for (const auto& value : clusterPower.first) {
                    if (value != 0) {
                        ClusterPower_lin_values.push_back(value);
                    }
                }

                SF_db_values.push_back(link.SF_db + CL);
                K_db_values.push_back(link.K_db);
                DS_sec_values.push_back(link.DS_sec);
                ASA_deg_values.push_back(link.ASA_deg);
                ASD_deg_values.push_back(link.ASD_deg);
                ZSA_deg_values.push_back(link.ZSA_deg);
                ZSD_deg_values.push_back(link.ZSD_deg);

                for (int n = 0; n < link.AOA_n_m.rows(); ++n) {
                    if (!link.isLos && n != 0) {
                        if (link.AOA_n_m(n, 0) != INFINITY || link.AOA_n_m(n, 0) != -INFINITY) {
                            AOD_values.insert(AOD_values.end(), link.AOD_n_m.row(n).begin(), link.AOD_n_m.row(n).end());
                            AOA_values.insert(AOA_values.end(), link.AOA_n_m.row(n).begin(), link.AOA_n_m.row(n).end());
                            ZOD_values.insert(ZOD_values.end(), link.ZOD_n_m.row(n).begin(), link.ZOD_n_m.row(n).end());
                            ZOA_values.insert(ZOA_values.end(), link.ZOA_n_m.row(n).begin(), link.ZOA_n_m.row(n).end());
                        }
                    }
                }
                countLink += 1;

            }

            plotCDF(SF_db_values, "SF, db", scenario, ep, 1, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([-40 100])");
            plotCDF(K_db_values, "K, db", scenario, ep, 2, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([-50 100])");
            plotCDF(DS_sec_values, "DS, sec", scenario, ep, 3, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1e-6])");
            plotCDF(ASA_deg_values, "ASA, deg", scenario, ep, 4, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1e3])");
            plotCDF(ASD_deg_values, "ASD, rad", scenario, ep, 5, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1.6])");
            plotCDF(ZSA_deg_values, "ZSA, deg", scenario, ep, 6, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 360])");
            plotCDF(ZSD_deg_values, "ZSD, deg", scenario, ep, 7, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 100])");
            plotCDF(ClusterDelay_sec_values, "Cluster Delays, sec", scenario, ep, 8, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1e-6])");
            plotCDF(ClusterPower_lin_values, "Relative Cluster Powers, linear scale", scenario, ep, 9, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1])");

            plotCDF(AOD_values, "AOD values, degree", scenario, ep, 10, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([-200 200])");
            plotCDF(AOA_values, "AOA values, degree", scenario, ep, 11, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([-200 200])");
            plotCDF(ZOD_values, "ZOD values, degree", scenario, ep, 12, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 180])");
            plotCDF(ZOA_values, "ZOA values, degree", scenario, ep, 13, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 180])");

            for (int rayIndex : RayUElist) {
                bool found = false;
                for (int coreIndex : indRayCorelist) {
                    if (coreIndex == rayIndex) {
                        found = true;
                        break;
                    }
                }
                coreUElist.push_back(found ? 1 : 0);
            }

            std::vector<double> coreArrPatt = filterVector(arrPatt, coreUElist);
            std::vector<double> maxPL_lin_magn = filterVector(PL_lin_magn, coreUElist);
            MatrixXd allArrPatt_magn = filterMatrixRows(arrPatt_magn, coreUElist);

            LinkBudget linkBudget(30, 7, 290);
            linkBudget.calculateMetrics(coreArrPatt, maxPL_lin_magn, allArrPatt_magn);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Elapsed time for scenario " << scenario << ": " << elapsed.count() << " s\n";
        }

        int userInput = 1;
        while (userInput) {
            std::cout << "Enter 0 to close all figures : ";
            std::cin >> userInput;

            if (userInput == 0) {
                // Close all figures in MATLAB
                engEvalString(ep, "close all;");
                engClose(ep);
                std::cout << "All figures closed." << std::endl;
            }
        }
        engClose(ep);
    }
    else {

        std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
        std::vector<std::string> frequencyBands = { "S", "Ka" };


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
            return 1;
        }
        double f = (frequencyChoice == 1) ? 2.0 : 20.0;

        std::ifstream file("setofangles.txt");  // Убедитесь, что путь к файлу правильный
        if (!file.is_open()) {
            std::cerr << "File no open!" << std::endl;
            return 1;
        }

        std::vector<double> data;
        double value;

        while (file >> value) {
            std::vector<double> row;
            data.push_back(value);
            while (file.peek() != '\n' && file >> value) {
                data.push_back(value);
            }
        }

        file.close();

        for (size_t i = 0; i < scenarios.size(); ++i) {
            const auto& scenario = scenarios[i];
            auto start = std::chrono::high_resolution_clock::now();

            std::vector<double> SF_db_values, K_db_values, DS_sec_values, ASA_deg_values, ASD_deg_values, ZSA_deg_values, ZSD_deg_values, ClusterDelay_sec_values, ClusterPower_lin_values;
            std::vector<double> AOD_values, AOA_values, ZOD_values, ZOA_values;

            TableLOS = GenerateMatrix(true, f, scenario);
            TableNLOS = GenerateMatrix(false, f, scenario);

            for (auto& angleEL : data) {
                int index = (AngleForLSP(angleEL)) / 10;
                LinkData link;
                link.elevationAngle = angleEL;
                link.isLos = CalculateLOSProbability(index, scenario);
                MatrixXd Table = link.isLos ? TableLOS : TableNLOS;


                //VectorXd Parameters = Table.col(index - 1);
                VectorXd Parameters = InterpolatedParameters(Table, angleEL);

                LSP::initializeParameters(link.isLos, link, Parameters);
                double CL = ChooseCL(link.isLos, f, angleEL, scenario);

                SSP::setParameters(Parameters);
                SSP::generateClusterDelays(link);
                std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);

                //std::vector<double>& Vector1 = clusterPower.first;
                //std::cout << "clusterPowers:\n";
                //for (size_t i = 0; i < Vector1.size(); ++i) {
                //    std::cout << "# " << i << ": " << Vector1[i] << ",";
                //}
                //std::cout << std::endl;


                ClusterDelay_sec_values.insert(ClusterDelay_sec_values.end(), link.clusterScaledDelays.begin(), link.clusterScaledDelays.end());
                //ClusterPower_lin_values.insert(ClusterPower_lin_values.end(), clusterPower.first.begin(), clusterPower.first.end());

                for (auto& value : clusterPower.first) {
                    if (value != 0.0) {
                        ClusterPower_lin_values.push_back(value);
                    }
                }

                //std::vector<double>& Vector2 = ClusterPower_lin_values;
                //std::cout << "ClusterPower_lin_values:\n";
                //for (size_t i = 0; i < Vector2.size(); ++i) {
                //    std::cout << "# " << i << ": " << Vector2[i] << ",";
                //}
                //std::cout << std::endl;

                // Сбор данных для LSP параметров
                SF_db_values.push_back(link.SF_db + CL);
                K_db_values.push_back(link.K_db);
                DS_sec_values.push_back(link.DS_sec);
                ASA_deg_values.push_back(link.ASA_deg);
                ASD_deg_values.push_back(link.ASD_deg);
                ZSA_deg_values.push_back(link.ZSA_deg);
                ZSD_deg_values.push_back(link.ZSD_deg);

            }
            plotCDF(SF_db_values, "SF, db", scenario, ep, 1, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([-40 100])");
            plotCDF(K_db_values, "K, db", scenario, ep, 2, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([-50 100])");
            plotCDF(DS_sec_values, "DS, sec", scenario, ep, 3, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1e-6])");
            plotCDF(ASA_deg_values, "ASA, deg", scenario, ep, 4, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1000])");
            plotCDF(ASD_deg_values, "ASD, rad", scenario, ep, 5, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1.6])");
            plotCDF(ZSA_deg_values, "ZSA, deg", scenario, ep, 6, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 360])");
            plotCDF(ZSD_deg_values, "ZSD, deg", scenario, ep, 7, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 100])");
            plotCDF(ClusterDelay_sec_values, "Cluster Delays, sec", scenario, ep, 8, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1e-6])");
            plotCDF(ClusterPower_lin_values, "Relative Cluster Powers, linear scale", scenario, ep, 9, frequencyChoice, frequencyBands);
            engEvalString(ep, "xlim([0 1])");

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Elapsed time for scenario " << scenario << ": " << elapsed.count() << " s\n";
        }
        int userInput = 1;
        while (userInput) {
            std::cout << "Enter 0 to close all figures : ";
            std::cin >> userInput;

            if (userInput == 0) {
                // Close all figures in MATLAB
                engEvalString(ep, "close all;");
                engClose(ep);
                std::cout << "All figures closed." << std::endl;
            }
        }
    }
    return 0;
}