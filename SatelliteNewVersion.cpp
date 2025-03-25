#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>
#include <engine.h>
#include <fstream>

#include "satellite_link.h"
#include "Tables.h"
#include "LOS_Probability.h"
#include "Pathloss.h"
#include "LargeScaleParameters.h"
#include "SmallScaleParameters.h"
#include "LinkBudget.h"
#include "matlab_plot.h"
#include "Antennas.h"
#include "ChannelMatrix.h"

void plotAntennaArray(Engine* ep, const Eigen::MatrixXd& antennaArray, const Eigen::Vector3d& UEorSat) {
    // Переносим антенную решетку в точку UEorSat
    Eigen::MatrixXd shiftedAntennaArray = antennaArray;
    for (int i = 0; i < antennaArray.cols(); ++i) {
        shiftedAntennaArray.col(i) += UEorSat;
    }

    // Передача данных в MATLAB
    mxArray* shiftedAntennaArrayMatlab = mxCreateDoubleMatrix(3, shiftedAntennaArray.cols(), mxREAL);
    std::memcpy(mxGetPr(shiftedAntennaArrayMatlab), shiftedAntennaArray.data(), 3 * shiftedAntennaArray.cols() * sizeof(double));
    engPutVariable(ep, "shiftedAntennaArray", shiftedAntennaArrayMatlab);

    mxArray* UEorSatMatlab = mxCreateDoubleMatrix(3, 1, mxREAL);
    std::memcpy(mxGetPr(UEorSatMatlab), UEorSat.data(), 3 * sizeof(double));
    engPutVariable(ep, "UEorSat", UEorSatMatlab);

    // Построение графика в MATLAB
    engEvalString(ep, "hold on;");
    engEvalString(ep, "plot3(shiftedAntennaArray(1,:), shiftedAntennaArray(2,:), shiftedAntennaArray(3,:)); hold on;");
    engEvalString(ep, "plot3(UEorSat(1), UEorSat(2), UEorSat(3), 'ro');");
    engEvalString(ep, "xlabel('X'); ylabel('Y'); zlabel('Z');");
    engEvalString(ep, "grid on; ");

    // Очистка памяти
    mxDestroyArray(shiftedAntennaArrayMatlab);
    mxDestroyArray(UEorSatMatlab);
}

void plotCDF(const std::vector<double>& data, const std::string& parameterName, const std::string& scenario, Engine* ep, int& figureNumber) {
    // Передача данных в MATLAB
    mxArray* mxData = mxCreateDoubleMatrix(1, data.size(), mxREAL);
    std::memcpy(mxGetPr(mxData), data.data(), data.size() * sizeof(double));
    engPutVariable(ep, "data", mxData);

    // Создание новой фигуры для каждого параметра
    std::string command =
        "figure(" + std::to_string(figureNumber) + "); "
        "[f, x] = ecdf(data); "
        "plot(x, f, 'DisplayName', '" + scenario + "'); "
        "title('CDF of " + parameterName + "'); "
        "xlabel('" + parameterName + "'); "
        "ylabel('CDF'); "
        "grid on; "
        "hold on; "  // Включаем режим добавления графиков
        "legend show;";

    engEvalString(ep, command.c_str());

    // Освобождение памяти
    mxDestroyArray(mxData);
    figureNumber++;
}

// Нужно переместить в Бюджет линии //
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
// Нужно переместить в Бюджет линии //

int main() {
    Engine* ep = engOpen("");
    if (!ep) {
        std::cerr << "MATLAB Engine no open" << std::endl;
        return 1;
    }
    engSetVisible(ep, false);

    int gear;
    int figureNumber;
    std::cout << "Chose gear: 0 - read file.txt and create Fig , 1 - Create model channel \n";
    std::cin >> gear;

    MatrixXd TableLOS;
    MatrixXd TableNLOS;

    if (gear) {

        std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
        std::vector<std::string> frequencyBands = { "S", "Ka" };

        int nTiersCore = 2, nTiresWrArnd = 2, nUePerCell;
        double satHeightKm = 600.0, elMinDegrees = 10.0, elTargetDegrees = 50.0, azTargetDegrees = 0.0;

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
            figureNumber = 1;
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



            // Нужно переместить в Бюджет линии //
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
            // Нужно переместить в Бюджет линии //

            std::vector<double> SF_db_values, K_db_values, DS_sec_values, ASA_deg_values, ASD_deg_values, ZSA_deg_values, ZSD_deg_values, ClusterDelay_sec_values, ClusterPower_lin_values;
            std::vector<double> AOD_values, AOA_values, ZOD_values, ZOA_values;

            TableLOS = GenerateMatrix(true, f, scenario);
            TableNLOS = GenerateMatrix(false, f, scenario);

            // for Antennas Sat and UE
            std::vector<double> anglesPol = { -45, 45 };
            Eigen::Vector2d antSatPattern;
            // for Antennas Sat and UE


            //MatlabPlot::plotEarth(ep);

            for (auto& link : UEswithSat.links.getLinks()) {

                //Инициализация параметров канала//
                int index = int(AngleForLSP(link.elevationAngle)) / 10;
                link.isLos = CalculateLOSProbability(index, scenario);
                MatrixXd Table = link.isLos ? TableLOS : TableNLOS;
                //VectorXd Parameters = Table.col(index - 1);
                VectorXd Parameters = InterpolatedParameters(Table, link.elevationAngle);
                //Инициализация параметров канала//

                //ЛСП//
                LSP::initializeParameters(link.isLos, link, Parameters);
                //ЛСП//

                //ССП//
                SSP::setParameters(Parameters);
                SSP::calculateLosAngles(link); // был минус у р2 , стоит перепроверить

                if (countLink == 0) {
                    std::cout << "AOA: " << link.AoA_Los << "," << "AOD: " << link.AoD_Los << "," << "ZOA: " << link.ZoA_Los << "," << "ZOD: " << link.ZoD_Los << "\n";

                }
                SSP::generateClusterDelays(link);
                std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);
                SSP::generateArrivalAndDepartureAngles(link, clusterPower.second);
                SSP::generateXPR(link);
                SSP::sortRelativeToFirstVector(link, clusterPower.first);
                SSP::transformVectors2MatForLink(link, clusterPower.first);
                //ССП//

                // Antenna UE for each link
                if (countLink == 0) {
                    AntennaArray::addAntennaSat(link, 2, 2, "Omni", f, anglesPol, 0.667);
                    //engEvalString(ep, "figure;");
                    //plotAntennaArray(ep, link.antennaSat.antennaArray * 1000, link.satellitePosition); 
                }
                else { link.antennaSat = UEswithSat.links.getLinks()[0].antennaSat; }
                AntennaArray::addAntennaUE(link, 2, 2, "Omni", f, anglesPol, 0.5);


                //plotAntennaArray(ep, link.antennaUE.antennaArray * 1000, link.userPosition);


                Eigen::MatrixXcd rayGainMatrix = ChannelMatrixH::generateRayGain(link);

                //std::cout << rayGainMatrix.size()<< "\n";

                /*auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "Elapsed time for scenario " << scenario << ": " << elapsed.count() << " s\n";*/
















                //Потери в канале//
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
                //Потери в канале//


                // Нужно переместить в Бюджет линии //
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
                // Нужно переместить в Бюджет линии //

                //Сбор данных//
                ClusterDelay_sec_values.insert(ClusterDelay_sec_values.end(), link.clusterScaledDelays.begin(), link.clusterScaledDelays.end());
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
                //Сбор данных//

                countLink += 1;

            }

            //Построение графиков//
            //plotCDF(SF_db_values, "SF, db", scenario, ep, 1);
            //plotCDF(K_db_values, "K, db", scenario, ep, 2);
            //plotCDF(DS_sec_values, "DS, sec", scenario, ep, 3);
            //plotCDF(ASA_deg_values, "ASA, deg", scenario, ep, 4);
            //plotCDF(ASD_deg_values, "ASD, deg", scenario, ep, 5);
            //plotCDF(ZSA_deg_values, "ZSA, deg", scenario, ep, 6);
            //plotCDF(ZSD_deg_values, "ZSD, deg", scenario, ep, 7);
            plotCDF(ClusterDelay_sec_values, "Cluster Delays, sec", scenario, ep, figureNumber);
            plotCDF(ClusterPower_lin_values, "Relative Cluster Powers, linear scale", scenario, ep, figureNumber);
            plotCDF(AOD_values, "AOD_values, degree", scenario, ep, figureNumber);
            plotCDF(AOA_values, "AOA_values, degree", scenario, ep, figureNumber);
            plotCDF(ZOD_values, "ZOD_values, degree", scenario, ep, figureNumber);
            plotCDF(ZOA_values, "ZOA_values, degree", scenario, ep, figureNumber);
            //Построение графиков//

            // Нужно переместить в Бюджет линии //
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
            // Нужно переместить в Бюджет линии //

            // Бюджет линии //
            LinkBudget linkBudget(30, 7, 290);
            linkBudget.calculateMetrics(coreArrPatt, maxPL_lin_magn, allArrPatt_magn);
            // Бюджет линии //

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
            figureNumber = 1;
            const auto& scenario = scenarios[i];
            auto start = std::chrono::high_resolution_clock::now();

            std::vector<double> SF_db_values, K_db_values, DS_sec_values, ASA_deg_values, ASD_deg_values, ZSA_deg_values, ZSD_deg_values, ClusterDelay_sec_values, ClusterPower_lin_values;
            std::vector<double> AOD_values, AOA_values, ZOD_values, ZOA_values;

            TableLOS = GenerateMatrix(true, f, scenario);
            TableNLOS = GenerateMatrix(false, f, scenario);

            for (auto& angleEL : data) {
                LinkData link;
                link.elevationAngle = angleEL;

                int index = (AngleForLSP(link.elevationAngle)) / 10;

                link.isLos = CalculateLOSProbability(index, scenario);
                MatrixXd Table = link.isLos ? TableLOS : TableNLOS;


                //VectorXd Parameters = Table.col(index - 1);
                VectorXd Parameters = InterpolatedParameters(Table, link.elevationAngle);

                LSP::initializeParameters(link.isLos, link, Parameters);
                double CL = ChooseCL(link.isLos, f, link.elevationAngle, scenario);

                SSP::setParameters(Parameters);
                SSP::generateClusterDelays(link);
                std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);

                //std::vector<double>& Vector1 = clusterPower.first;
                //std::cout << "clusterPowers:\n";
                //for (size_t i = 0; i < Vector1.size(); ++i) {
                //    std::cout << "# " << i << ": " << Vector1[i] << ",";
                //}
                //std::cout << std::endl;
                if (true) {
                    ClusterDelay_sec_values.insert(ClusterDelay_sec_values.end(), link.clusterScaledDelays.begin(), link.clusterScaledDelays.end());
                    //ClusterPower_lin_values.insert(ClusterPower_lin_values.end(), clusterPower.first.begin(), clusterPower.first.end());

                    for (auto& value : clusterPower.first) {
                        if (value != 0.0) {
                            ClusterPower_lin_values.push_back((value));
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

            }

            std::cout << SF_db_values.size() << "\n";
            plotCDF(SF_db_values, "SF, db", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([-40 100])");
            plotCDF(K_db_values, "K, db", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([-50 100])");
            plotCDF(DS_sec_values, "DS, sec", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([0 1e-6])");
            plotCDF(ASA_deg_values, "ASA, deg", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([0 1000])");
            plotCDF(ASD_deg_values, "ASD, deg", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([0 1.6])");
            plotCDF(ZSA_deg_values, "ZSA, deg", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([0 350])");
            plotCDF(ZSD_deg_values, "ZSD, deg", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([0 100])");
            plotCDF(ClusterDelay_sec_values, "Cluster Delays, sec", scenario, ep, figureNumber);
            engEvalString(ep, "xlim([0 1e-6])");
            plotCDF(ClusterPower_lin_values, "Relative Cluster Powers, linear scale", scenario, ep, figureNumber);

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