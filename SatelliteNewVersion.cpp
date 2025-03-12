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
#include "matlab_plot.h"


void plotCDF(const std::vector<double>& data, const std::string& parameterName, const std::string& scenario, Engine* ep, int figureNumber) {
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


            MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
            MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);

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


            TableLOS = GenerateMatrix(true, f, scenario);
            TableNLOS = GenerateMatrix(false, f, scenario);

            for (auto& link : UEswithSat.links.getLinks()) {

                int index = int(AngleForLSP(link.elevationAngle)) / 10;
                link.isLos = CalculateLOSProbability(index, scenario);

                MatrixXd Table = link.isLos ? TableLOS : TableNLOS;
                VectorXd Parameters = Table.col(index - 1);

                LSP::initializeParameters(link.isLos, link, Parameters);

                SSP::setParameters(Parameters);
                SSP::calculateLosAngles(link, UEswithSat.p1, UEswithSat.p2);
                SSP::generateClusterDelays(link);
                std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);
                SSP::generateArrivalAndDepartureAngles(link, clusterPower.second);



                SSP::generateXPR(link);

                std::vector<double>& secondVector = clusterPower.first;
                std::cout << "clusterPowers:\n";
                for (size_t i = 0; i < secondVector.size(); ++i) {
                    std::cout << "# " << i << ": " << secondVector[i] << ",";
                }
                std::cout << std::endl;

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



                ClusterDelay_sec_values.insert(ClusterDelay_sec_values.end(), link.clusterDelays.begin(), link.clusterDelays.end());
                //ClusterPower_lin_values.insert(ClusterPower_lin_values.end(), clusterPower.first.begin(), clusterPower.first.end());

                for (const auto& value : clusterPower.first) {
                    if (value != 0.0) {
                        ClusterPower_lin_values.push_back(value);
                    }
                }

                SF_db_values.push_back(link.SF_db + CL);
                if (link.isLos) { K_db_values.push_back(link.K_db); }
                else { K_db_values.push_back(-INFINITY); }
                DS_sec_values.push_back(link.DS_sec);
                ASA_deg_values.push_back(link.ASA_deg);
                ASD_deg_values.push_back(link.ASD_deg);
                ZSA_deg_values.push_back(link.ZSA_deg);
                ZSD_deg_values.push_back(link.ZSD_deg);

                countLink += 1;

            }


            plotCDF(SF_db_values, "SF, db", scenario, ep, 1);
            plotCDF(K_db_values, "K, db", scenario, ep, 2);
            plotCDF(DS_sec_values, "DS, sec", scenario, ep, 3);
            plotCDF(ASA_deg_values, "ASA, deg", scenario, ep, 4);
            plotCDF(ASD_deg_values, "ASD, deg", scenario, ep, 5);
            plotCDF(ZSA_deg_values, "ZSA, deg", scenario, ep, 6);
            plotCDF(ZSD_deg_values, "ZSD, deg", scenario, ep, 7);
            plotCDF(ClusterDelay_sec_values, "Cluster Delays, sec", scenario, ep, 8);
            plotCDF(ClusterPower_lin_values, "Relative Cluster Powers, linear scale", scenario, ep, 9);

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

            TableLOS = GenerateMatrix(true, f, scenario);
            TableNLOS = GenerateMatrix(false, f, scenario);

            for (auto& angleEL : data) {
                int index = int(AngleForLSP(angleEL)) / 10;
                LinkData link;
                double isLos = CalculateLOSProbability(index, scenario);

                MatrixXd Table = isLos ? TableLOS : TableNLOS;
                VectorXd Parameters = Table.col(index - 1);

                LSP::initializeParameters(isLos, link, Parameters);
                double CL = ChooseCL(isLos, f, angleEL, scenario);

                SSP::setParameters(Parameters);
                SSP::generateClusterDelays(link);
                std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);



                ClusterDelay_sec_values.insert(ClusterDelay_sec_values.end(), link.clusterDelays.begin(), link.clusterDelays.end());
                //ClusterPower_lin_values.insert(ClusterPower_lin_values.end(), clusterPower.first.begin(), clusterPower.first.end());
                for (const auto& value : clusterPower.first) {
                    if (value != 0.0) {
                        ClusterPower_lin_values.push_back(value);
                    }
                }
                // Сбор данных для LSP параметров
                SF_db_values.push_back(link.SF_db + CL);
                if (isLos) { K_db_values.push_back(link.K_db); }
                else { K_db_values.push_back(-INFINITY); }
                DS_sec_values.push_back(link.DS_sec);
                ASA_deg_values.push_back(link.ASA_deg);
                ASD_deg_values.push_back(link.ASD_deg);
                ZSA_deg_values.push_back(link.ZSA_deg);
                ZSD_deg_values.push_back(link.ZSD_deg);

            }
            plotCDF(SF_db_values, "SF, db", scenario, ep, 1);
            plotCDF(K_db_values, "K, db", scenario, ep, 2);
            plotCDF(DS_sec_values, "DS, sec", scenario, ep, 3);
            plotCDF(ASA_deg_values, "ASA, deg", scenario, ep, 4);
            plotCDF(ASD_deg_values, "ASD, deg", scenario, ep, 5);
            plotCDF(ZSA_deg_values, "ZSA, deg", scenario, ep, 6);
            plotCDF(ZSD_deg_values, "ZSD, deg", scenario, ep, 7);
            plotCDF(ClusterDelay_sec_values, "Cluster Delays, sec", scenario, ep, 8);
            plotCDF(ClusterPower_lin_values, "Relative Cluster Powers, linear scale", scenario, ep, 9);

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

//#include <iostream>
//#include <Eigen/Dense>
//#include <vector>
//#include <chrono>
//
//// для сортировки векторов
//#include <algorithm>
//
//#include "NTN_Deployment.h"
//#include "Tables.h"
//#include "Matlab_plot.h"
//#include "LOS_Probability.h"
//#include "LargeScaleParameters.h"
//#include "Pathloss.h"
//#include "SSP.h"
//#include "Budget.h"
//
//
//// Функция для вычисления диаграммы направленности
//void calculateDishPattern(Eigen::VectorXd& teta_rad, Eigen::VectorXd& AP_dB, double rDish_WL) {
//    for (int i = 0; i < teta_rad.size(); ++i) {
//        double teta = teta_rad(i);
//        if (teta == 0) {
//            AP_dB(i) = 1.0;
//        }
//        else {
//            double bessel_arg = 2 * M_PI * rDish_WL * sin(teta);
//            double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
//            double res = 4 * std::pow(bessel_val / bessel_arg, 2);
//            double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;
//            AP_dB(i) = 10 * log10(res) + Gain;
//
//        }
//    }
//}
//
//// Функция для поиска индексов совпадающих элементов
//std::vector<int> findMatchingIndices(const std::vector<Eigen::Vector2d>& uvSet, const std::vector<Eigen::Vector2d>& uvSetCore) {
//    std::vector<int> matchingIndices;
//
//    // Создаем unordered_set из элементов uvSetCore для быстрого поиска
//    std::unordered_set<Eigen::Vector2d, Vector2dHash> coreSet(uvSetCore.begin(), uvSetCore.end());
//
//    for (int i = 0; i < uvSet.size(); ++i) {
//        if (coreSet.find(uvSet[i]) != coreSet.end()) {
//            matchingIndices.push_back(i);
//        }
//    }
//
//    return matchingIndices;
//}
//
//// функция для маски ( фильтра)
//std::vector<double> filterVector(const std::vector<double>& inputVector, const std::vector<int>& mask) {
//    std::vector<double> filteredVector;
//
//    // Проверяем, что размеры векторов совпадают
//    if (inputVector.size() != mask.size()) {
//        std::cerr << "Error: inputVector and mask have different sizes!" << std::endl;
//        return filteredVector; // Возвращаем пустой вектор
//    }
//
//    // Фильтруем элементы
//    for (size_t i = 0; i < mask.size(); ++i) {
//        if (mask[i] == 1) {
//            filteredVector.push_back(inputVector[i]);
//        }
//    }
//
//    return filteredVector;
//}
//
//// Функция для фильтрации строк матрицы
//MatrixXd filterMatrixRows(const MatrixXd& inputMatrix, const std::vector<int>& mask) {
//    // Проверяем, что размер маски соответствует количеству строк в матрице
//    if (mask.size() != inputMatrix.rows()) {
//        std::cerr << "Error: Mask size does not match the number of rows in the matrix!" << std::endl;
//        return MatrixXd(); // Возвращаем пустую матрицу
//    }
//
//    // Считаем количество строк, которые нужно сохранить
//    int numRowsToKeep = 0;
//    for (int value : mask) {
//        if (value == 1) {
//            numRowsToKeep++;
//        }
//    }
//
//    // Создаем новую матрицу для хранения отфильтрованных строк
//    MatrixXd filteredMatrix(numRowsToKeep, inputMatrix.cols());
//
//    // Копируем строки, соответствующие маске
//    int filteredRowIndex = 0;
//    for (int i = 0; i < mask.size(); ++i) {
//        if (mask[i] == 1) {
//            filteredMatrix.row(filteredRowIndex) = inputMatrix.row(i);
//            filteredRowIndex++;
//        }
//    }
//    return filteredMatrix;
//}
//
//
//int main() {
//    MatlabPlot plotter;
//    //plotter.plotEarth();
//
//    int mode;
//    std::cout << "Enter a number 1 or 2: ";
//    std::cin >> mode;
//
//    while (true)
//    {
//        if (mode == 1)
//        {
//            double altitude = 600.0;
//
//            while (true)
//            {
//
//                std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
//                std::vector<std::string> frequencyBands = { "S", "Ka" };
//
//                std::cout << "Select a scenario:\n";
//                for (size_t i = 0; i < scenarios.size(); ++i) {
//                    std::cout << i + 1 << ". " << scenarios[i] << "\n";
//                }
//                std::cout << "Enter scenario number (1-" << scenarios.size() << ") or 0 to exit: ";
//
//                int scenarioChoice;
//                std::cin >> scenarioChoice;
//
//                if (scenarioChoice == 0) {
//                    break;
//                }
//                else if (scenarioChoice < 1 || scenarioChoice > scenarios.size()) {
//                    std::cout << "Invalid choice. Please try again.\n";
//                    continue;
//                }
//
//                std::string scenario = scenarios[scenarioChoice - 1];
//
//                // Frequency band selection
//                std::cout << "Select a frequency band:\n";
//                for (size_t i = 0; i < frequencyBands.size(); ++i) {
//                    std::cout << i + 1 << ". " << frequencyBands[i] << " band\n";
//                }
//                std::cout << "Enter frequency band number (1-" << frequencyBands.size() << "): ";
//
//                int frequencyChoice;
//                std::cin >> frequencyChoice;
//
//                if (frequencyChoice < 1 || frequencyChoice > frequencyBands.size()) {
//                    std::cout << "Invalid choice. Please try again.\n";
//                    continue; // Repeat the loop
//                }
//
//                double f = (frequencyChoice == 1) ? 2.0 : 20.0; // Example: 2 GHz for S and 20 GHz for Ka
//
//                double elTargetDegrees;
//                std::cout << "Enter elTargetDegrees  [grad]  (90 - nadir) : ";
//                std::cin >> elTargetDegrees;
//                double azTargetDegrees;
//                std::cout << "Enter azTargetDegrees  [grad]   : ";
//                std::cin >> azTargetDegrees;
//
//
//
//                auto start = std::chrono::high_resolution_clock::now();
//
//                //for (int i = 0; i < 10; ++i) {
//                SatelliteLink UEswithSat(2, 2, 20, altitude, 10.0, elTargetDegrees, azTargetDegrees);
//                UEswithSat.generateLinks();
//                //}
//
//                auto end = std::chrono::high_resolution_clock::now();
//                std::chrono::duration<double> elapsed = end - start;
//                std::cout << "Elapsed time: " << elapsed.count() << " s\n";
//
//
//                // Пользователи на сфере
//                //std::vector<Eigen::Vector3d> users;
//                //for (const auto& link : UEswithSat.links.getLinks()) {
//                //    //std::cout << link.userPosition;
//                //    users.push_back(link.userPosition);
//                //}
//
//                //plotter.plotTransformedData(users, UEswithSat.links.getLinks()[0].satellitePosition);
//
//                // Generate tables for LOS and NLOS
//                MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
//                MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);
//
//
//                std::vector<Eigen::Vector3d> rRays = UEswithSat.getVectorsToCellCenters(); // Лучи к центрам ячеек 
//                MatrixXd arrPatt_magn(UEswithSat.links.getLinks().size(), rRays.size()); // Массив всех амплитуд  для каждого пользователя от любого центра ячейки
//                std::vector<double> arrPatt; // ДН максимальные по всем лучам
//                std::vector<double> PL_lin_magn; // PL по всем лучам
//
//                double rDish_WL = 2e9 / 299792458;
//                int countLink = 0;
//
//                std::vector<int> indRayCorelist = findMatchingIndices(UEswithSat.uvSet, UEswithSat.uvSetCore);
//                std::cout << "\n Count indUvCorelist: " << indRayCorelist.size() << std::endl;
//                std::vector<int>  RayUElist;
//                std::vector<int> coreUElist;
//
//
//                Eigen::Matrix3d rotMatrix = UEswithSat.rotMatrix;
//
//                for (auto& link : UEswithSat.links.getLinks())
//                {
//
//                    double deg = link.elevationAngle;
//
//                    int index = int(AngleForLSP(deg)) / 10;
//
//                    link.isLos = CalculateLOSProbability(index, scenario);
//
//
//                    MatrixXd Table = link.isLos ? TableLOS : TableNLOS;
//                    VectorXd Parameters = Table.col(index - 1);
//
//                    LSP::initializeParameters(link.isLos, link, Parameters);
//                    std::cout << link.isLos << "\n";
//                    SSP::setParameters(Parameters);
//                    SSP::calculateLosAngles(link, UEswithSat.p1, UEswithSat.p2);
//                    SSP::generateClusterDelays(link);
//                    std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);
//
//                    //Вывод мощностей ( без / с масштабирующим фактором) 
//                    //std::cout << "clusterPower: " << clusterPower.first[0] << ", " << clusterPower.first[1] << ", " << clusterPower.first[2] << "\n";
//                    //std::cout << "clusterPowersWithScalingFactors: " << clusterPower.second[0] << ", " << clusterPower.second[1] << ", " << clusterPower.second[2] << "\n";
//
//
//                    /////////////////////////////////////////Нужно позже оформить ...  ниже почти всё для расчёта бюджета линии////////////////////////////////////
//
//                    // для PL 
//                    double d = CalculateDistance(EARTH_RADIUS, altitude, deg * PI / 180);
//                    double std = ChooseSTD(link.isLos, f, deg, scenario);
//                    double SF = GenerateSF(std);
//                    double FSPL = Calculate_FSPL(d * 1e3, f);
//                    double CL = ChooseCL(link.isLos, f, deg, scenario);
//                    double PL_dB = CalculateBasisPathLoss(FSPL, SF, CL);
//
//
//                    // генератор для шума
//                    std::random_device rd;  // Источник энтропии
//                    std::mt19937 gen(rd()); // Генератор случайных чисел (Mersenne Twister)
//                    std::normal_distribution<> dist(0, 0.72); // Нормальное распределение с mean=0 и stddev=0.72
//                    PL_dB = PL_dB + dist(gen); // добавляем SF
//
//                    // Перевод в линейный масштаб
//                    double PL_lin_magnitude = pow(10, -1 * PL_dB / 20);
//                    //std::cout << "\nPL_lin_magnitude: " << PL_lin_magnitude << std::endl;
//
//
//
//                    for (int i = 0; i < static_cast<int> (rRays.size()); ++i)
//                    {
//
//                        Eigen::Vector3d rotatedVector = rotMatrix * (link.userPosition - link.satellitePosition).normalized();
//
//                        double theta = std::acos(rRays[i].dot(rotatedVector));
//                        double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;
//
//                        if (theta == 0) {
//                            arrPatt_magn(countLink, i) = pow(10.0, (1.0 + Gain) / 20.0);
//                        }
//                        else {
//                            double bessel_arg = 2 * M_PI * rDish_WL * sin(theta);
//                            double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
//                            double res = 4 * std::pow(bessel_val / bessel_arg, 2);
//                            arrPatt_magn(countLink, i) = pow(10.0, (10 * log10(res) + Gain) / 20.0);
//
//                        }
//                    }
//
//                    PL_lin_magn.push_back(PL_lin_magnitude);
//                    arrPatt.push_back(arrPatt_magn.row(countLink).maxCoeff(&RayUElist.emplace_back()));
//
//                    countLink += 1;
//                }
//
//                // Вывод задержек во времени
//                for (auto& link : UEswithSat.links.getLinks()) {
//                    std::cout << "Cluster Delay: " << link.clusterDelays[0] << ", " << link.clusterDelays[1] << ", " << link.clusterDelays[2] << ", ";
//                    if (!link.isLos) { std::cout << link.clusterDelays[3]; }
//                    std::cout << " \n";
//                }
//
//
//                // Заполняем coreUElist
//
//                for (int rayIndex : RayUElist) {
//                    bool found = false;
//                    for (int coreIndex : indRayCorelist) {
//                        if (coreIndex == rayIndex) {
//                            found = true;
//                            break;
//                        }
//                    }
//                    coreUElist.push_back(found ? 1 : 0); // Добавляем 1, если элемент найден, иначе 0
//                }
//
//                std::vector<double> coreArrPatt = filterVector(arrPatt, coreUElist);
//                std::vector<double> maxPL_lin_magn = filterVector(PL_lin_magn, coreUElist);
//                MatrixXd allArrPatt_magn = filterMatrixRows(arrPatt_magn, coreUElist);
//
//                LinkBudget linkBudget(30, 7, 290); //  BW_MHz = 30 МГц, NF_dB = 7 дБ, T = 290 К
//                linkBudget.calculateMetrics(coreArrPatt, maxPL_lin_magn, allArrPatt_magn);
//                plotter.plotCDF(linkBudget.getDL_CNR_dB(), linkBudget.getDL_CIR_dB(), linkBudget.getDL_CINR_dB(), scenario);
//
//            }
//            return 0;
//
//            //plotter.plotRayPoints(UEswithSat.links.getLinks()[0].satellitePosition, rRays);
//
//
//            // Вычисление диаграммы направленности
//            //calculateDishPattern(elTargetDegrees, AP_dB, rDish_WL);
//        }
//
//
//        // mode 2 - построенеи графиков для LSP
//
// else
// {
//     // 1) разделить на Los и NLOS
//     // 2) занулить ASD, ZSD
//     std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
//     std::vector<std::string> frequencyBands = { "S", "Ka" };
//
//     // Frequency band selection
//     std::cout << "Select a frequency band:\n";
//     for (size_t i = 0; i < frequencyBands.size(); ++i) {
//         std::cout << i + 1 << ". " << frequencyBands[i] << " band\n";
//     }
//     std::cout << "Enter frequency band number (1-" << frequencyBands.size() << "): ";
//
//     int frequencyChoice;
//     std::cin >> frequencyChoice;
//
//     double f = (frequencyChoice == 1) ? 2.0 : 30.0; // Example: 2 GHz for S and 30 GHz for Ka
//
//     double altitude = 600.0;
//     int z = 5;
//     int level = 2;
//     int Nueincell = 20;
//
//     // генерация спутников для сбора статистики
//
//     //DenseUrban scenario
//
//     std::vector<double> SFDU;
//     std::vector<double> KDU;
//     std::vector<double> DSDU;
//     std::vector<double> ASDDU;
//     std::vector<double> ASADU;
//     std::vector<double> ZSDDU;
//     std::vector<double> ZSADU;
//
//     for (int i = 0; i < z; ++i)
//     {
//         double elTargetDegrees = 90;
//         double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);
//         int numLos = 0;
//         SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
//         UEswithSat.generateLinks();
//
//         MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[0]);
//         MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[0]);
//
//         for (auto& link : UEswithSat.links.getLinks()) {
//             double deg = link.elevationAngle;
//             int index = int(AngleForLSP(deg)) / 10;
//             bool isLos = CalculateLOSProbability(index, scenarios[0]);
//             
//             MatrixXd Table = isLos ? TableLOS : TableNLOS;
//             VectorXd Parameters = Table.col(index - 1);
//             LSP::initializeParameters(isLos, link, Parameters);
//
//             if (isLos)
//             {
//                 SFDU.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[0]));
//                 KDU.push_back(link.K_db);
//                 DSDU.push_back(pow(10, link.DS_sec));
//                 ASDDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                 ASADU.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                 ZSDDU.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                 ZSADU.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                 numLos += 1;
//             }
//             else
//             {
//                 SFDU.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[0]));
//                 DSDU.push_back(pow(10, link.DS_sec));
//                 ASDDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                 ASADU.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                 ZSDDU.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                 ZSADU.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//             }
//
//         }
//         
//
//         std::sort(SFDU.data(), SFDU.data() + SFDU.size());
//         std::sort(KDU.data(), KDU.data() + KDU.size());
//         std::sort(DSDU.data(), DSDU.data() + DSDU.size());
//         std::sort(ASDDU.data(), ASDDU.data() + ASDDU.size());
//         std::sort(ASADU.data(), ASADU.data() + ASADU.size());
//         std::sort(ZSDDU.data(), ZSDDU.data() + ZSDDU.size());
//         std::sort(ZSADU.data(), ZSADU.data() + ZSADU.size());
//
//         std::cout << numLos << std::endl;
//     }
//
//
//     // Urban scenario
//
//     std::vector<double> SFU;
//     std::vector<double> KU;
//     std::vector<double> DSU;
//     std::vector<double> ASDU;
//     std::vector<double> ASAU;
//     std::vector<double> ZSDU;
//     std::vector<double> ZSAU;
//
//    
//     for (int i = 0; i < z; ++i)
//     {
//         double elTargetDegrees = 90; // попрпобовать разные углы спутников 
//         double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);
//         int numLos = 0;
//         SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
//         UEswithSat.generateLinks();
//
//         MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[1]);
//         MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[1]);
//
//         for (int i = 0; i < 1; ++i)
//         {
//             for (auto& link : UEswithSat.links.getLinks()) {
//                 double deg = link.elevationAngle;
//                 int index = int(AngleForLSP(deg)) / 10;
//                 bool isLos = CalculateLOSProbability(index, scenarios[1]);
//
//                 MatrixXd Table = isLos ? TableLOS : TableNLOS;
//                 VectorXd Parameters = Table.col(index - 1);
//                 LSP::initializeParameters(isLos, link, Parameters);
//
//                 if (isLos)
//                 {
//                     SFDU.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[1]));
//                     KU.push_back(link.K_db);
//                     DSU.push_back(pow(10, link.DS_sec));
//                     ASDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                     ASAU.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                     ZSDU.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                     ZSAU.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                     numLos += 1;
//                 }
//                 else
//                 {
//                     SFU.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[1]));
//                     DSU.push_back(pow(10, link.DS_sec));
//                     ASDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                     ASAU.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                     ZSDU.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                     ZSAU.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                 }
//
//             }
//
//         }
//
//         std::sort(SFU.data(), SFU.data() + SFU.size());
//         std::sort(KU.data(), KU.data() + KU.size());
//         std::sort(DSU.data(), DSU.data() + DSU.size());
//         std::sort(ASDU.data(), ASDU.data() + ASDU.size());
//         std::sort(ASAU.data(), ASAU.data() + ASAU.size());
//         std::sort(ZSDU.data(), ZSDU.data() + ZSDU.size());
//         std::sort(ZSAU.data(), ZSAU.data() + ZSAU.size());
//
//         std::cout << numLos << std::endl;
//
//     }
//     
//
//
//     // Suburban scenario
//
//     std::vector<double> SFS;
//     std::vector<double> KS;
//     std::vector<double> DSS;
//     std::vector<double> ASDS;
//     std::vector<double> ASAS;
//     std::vector<double> ZSDS;
//     std::vector<double> ZSAS;
//
//
//     for (int i = 0; i < z; ++i)
//     {
//         double elTargetDegrees = 90;
//         double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);
//         int numLos = 0;
//         SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
//         UEswithSat.generateLinks();
//
//         MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[2]);
//         MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[2]);
//
//         for (int i = 0; i < 1; ++i)
//         {
//             for (auto& link : UEswithSat.links.getLinks()) {
//                 double deg = link.elevationAngle;
//                 int index = int(AngleForLSP(deg)) / 10;
//                 bool isLos = CalculateLOSProbability(index, scenarios[2]);
//
//                 MatrixXd Table = isLos ? TableLOS : TableNLOS;
//                 VectorXd Parameters = Table.col(index - 1);
//                 LSP::initializeParameters(isLos, link, Parameters);
//
//                 if (isLos)
//                 {
//                     SFDU.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[2]));
//                     KS.push_back(link.K_db);
//                     DSS.push_back(pow(10, link.DS_sec));
//                     ASDS.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                     ASAS.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                     ZSDS.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                     ZSAS.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                     numLos += 1;
//                 }
//                 else
//                 {
//                     SFS.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[2]));
//                     DSS.push_back(pow(10, link.DS_sec));
//                     ASDS.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                     ASAS.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                     ZSDS.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                     ZSAS.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                 }
//
//             }
//
//         }
//
//         std::sort(SFS.data(), SFS.data() + SFS.size());
//         std::sort(KS.data(), KS.data() + KS.size());
//         std::sort(DSS.data(), DSS.data() + DSS.size());
//         std::sort(ASDS.data(), ASDS.data() + ASDS.size());
//         std::sort(ASAS.data(), ASAS.data() + ASAS.size());
//         std::sort(ZSDS.data(), ZSDS.data() + ZSDS.size());
//         std::sort(ZSAS.data(), ZSAS.data() + ZSAS.size());
//
//         std::cout << numLos << std::endl;
//     }
//
//
//     // Rural scenario
//     std::vector<double> SFR;
//     std::vector<double> KR;
//     std::vector<double> DSR;
//     std::vector<double> ASDR;
//     std::vector<double> ASAR;
//     std::vector<double> ZSDR;
//     std::vector<double> ZSAR;
//
//     for (int i = 0; i < z; ++i)
//     {
//         double elTargetDegrees = 90;
//         double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);
//         int numLos = 0;
//         SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
//         UEswithSat.generateLinks();
//
//         MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[3]);
//         MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[3]);
//
//         for (int i = 0; i < 1; ++i)
//         {
//             for (auto& link : UEswithSat.links.getLinks()) {
//                 double deg = link.elevationAngle;
//                 int index = int(AngleForLSP(deg)) / 10;
//                 bool isLos = CalculateLOSProbability(index, scenarios[3]);
//
//                 MatrixXd Table = isLos ? TableLOS : TableNLOS;
//                 VectorXd Parameters = Table.col(index - 1);
//                 LSP::initializeParameters(isLos, link, Parameters);
//
//                 if (isLos)
//                 {
//                     SFDU.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[3]));
//                     KR.push_back(link.K_db);
//                     DSR.push_back(pow(10, link.DS_sec));
//                     ASDR.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                     ASAR.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                     ZSDR.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                     ZSAR.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                     numLos += 1;
//                 }
//                 else
//                 {
//                     SFR.push_back(link.SF_db + ChooseCL(isLos, f, link.elevationAngle, scenarios[3]));
//                     DSR.push_back(pow(10, link.DS_sec));
//                     ASDR.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
//                     ASAR.push_back(M_PI / 180 * pow(10, std::min(link.ASA_deg, log10(90))));
//                     ZSDR.push_back(M_PI / 180 * pow(10, std::min(link.ZSD_deg, log10(360))));
//                     ZSAR.push_back(M_PI / 180 * pow(10, std::min(link.ZSA_deg, log10(360))));
//                 }
//
//             }
//
//         }
//
//         std::sort(SFR.data(), SFR.data() + SFR.size());
//         std::sort(KR.data(), KR.data() + KR.size());
//         std::sort(DSR.data(), DSR.data() + DSR.size());
//         std::sort(ASDR.data(), ASDR.data() + ASDR.size());
//         std::sort(ASAR.data(), ASAR.data() + ASAR.size());
//         std::sort(ZSDR.data(), ZSDR.data() + ZSDR.size());
//         std::sort(ZSAR.data(), ZSAR.data() + ZSAR.size());
//
//         std::cout << numLos << std::endl;
//     }
//
//     MatlabLSPPlot plotter;
//     plotter.plotForAllLSP(SFDU, SFU, SFS, SFR, KDU, KU, KS, KR, DSDU, DSU, DSS, DSR, ASDDU, ASDU, ASDS, ASDR, ASADU, ASAU, ASAS, ASAR, ZSDDU, ZSDU, ZSDS, ZSDR, ZSADU, ZSAU, ZSAS, ZSAR, scenarios, frequencyBands[frequencyChoice - 1]);
//
//
//     std::cin;
//    }
//
//    }
//
//}
