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
#include "SSP.h"
#include "Budget.h"


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

// Функция для поиска индексов совпадающих элементов
std::vector<int> findMatchingIndices(const std::vector<Eigen::Vector2d>& uvSet, const std::vector<Eigen::Vector2d>& uvSetCore) {
    std::vector<int> matchingIndices;

    // Создаем unordered_set из элементов uvSetCore для быстрого поиска
    std::unordered_set<Eigen::Vector2d, Vector2dHash> coreSet(uvSetCore.begin(), uvSetCore.end());

    for (int i = 0; i < uvSet.size(); ++i) {
        if (coreSet.find(uvSet[i]) != coreSet.end()) {
            matchingIndices.push_back(i);
        }
    }

    return matchingIndices;
}

// функция для маски ( фильтра)
std::vector<double> filterVector(const std::vector<double>& inputVector, const std::vector<int>& mask) {
    std::vector<double> filteredVector;

    // Проверяем, что размеры векторов совпадают
    if (inputVector.size() != mask.size()) {
        std::cerr << "Error: inputVector and mask have different sizes!" << std::endl;
        return filteredVector; // Возвращаем пустой вектор
    }

    // Фильтруем элементы
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i] == 1) {
            filteredVector.push_back(inputVector[i]);
        }
    }

    return filteredVector;
}

// Функция для фильтрации строк матрицы
MatrixXd filterMatrixRows(const MatrixXd& inputMatrix, const std::vector<int>& mask) {
    // Проверяем, что размер маски соответствует количеству строк в матрице
    if (mask.size() != inputMatrix.rows()) {
        std::cerr << "Error: Mask size does not match the number of rows in the matrix!" << std::endl;
        return MatrixXd(); // Возвращаем пустую матрицу
    }

    // Считаем количество строк, которые нужно сохранить
    int numRowsToKeep = 0;
    for (int value : mask) {
        if (value == 1) {
            numRowsToKeep++;
        }
    }

    // Создаем новую матрицу для хранения отфильтрованных строк
    MatrixXd filteredMatrix(numRowsToKeep, inputMatrix.cols());

    // Копируем строки, соответствующие маске
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
    MatlabPlot plotter;
    //plotter.plotEarth();

    int mode;
    std::cout << "Enter a number 1 or 2: ";
    std::cin >> mode;

    while (true)
    {
        if (mode == 1)
        {
            double altitude = 600.0;

            while (true)
            {

                std::vector<std::string> scenarios = { "DenseUrban", "Urban", "Suburban", "Rural" };
                std::vector<std::string> frequencyBands = { "S", "Ka" };

                std::cout << "Select a scenario:\n";
                for (size_t i = 0; i < scenarios.size(); ++i) {
                    std::cout << i + 1 << ". " << scenarios[i] << "\n";
                }
                std::cout << "Enter scenario number (1-" << scenarios.size() << ") or 0 to exit: ";

                int scenarioChoice;
                std::cin >> scenarioChoice;

                if (scenarioChoice == 0) {
                    break;
                }
                else if (scenarioChoice < 1 || scenarioChoice > scenarios.size()) {
                    std::cout << "Invalid choice. Please try again.\n";
                    continue;
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

                double f = (frequencyChoice == 1) ? 2.0 : 20.0; // Example: 2 GHz for S and 20 GHz for Ka

                double elTargetDegrees;
                std::cout << "Enter elTargetDegrees  [grad]  (90 - nadir) : ";
                std::cin >> elTargetDegrees;
                double azTargetDegrees;
                std::cout << "Enter azTargetDegrees  [grad]   : ";
                std::cin >> azTargetDegrees;



                auto start = std::chrono::high_resolution_clock::now();

                //for (int i = 0; i < 10; ++i) {
                SatelliteLink UEswithSat(2, 2, 20, altitude, 10.0, elTargetDegrees, azTargetDegrees);
                UEswithSat.generateLinks();
                //}

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                std::cout << "Elapsed time: " << elapsed.count() << " s\n";


                // Пользователи на сфере
                //std::vector<Eigen::Vector3d> users;
                //for (const auto& link : UEswithSat.links.getLinks()) {
                //    //std::cout << link.userPosition;
                //    users.push_back(link.userPosition);
                //}

                //plotter.plotTransformedData(users, UEswithSat.links.getLinks()[0].satellitePosition);

                // Generate tables for LOS and NLOS
                MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
                MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);


                std::vector<Eigen::Vector3d> rRays = UEswithSat.getVectorsToCellCenters(); // Лучи к центрам ячеек 
                MatrixXd arrPatt_magn(UEswithSat.links.getLinks().size(), rRays.size()); // Массив всех амплитуд  для каждого пользователя от любого центра ячейки
                std::vector<double> arrPatt; // ДН максимальные по всем лучам
                std::vector<double> PL_lin_magn; // PL по всем лучам

                double rDish_WL = 2e9 / 299792458;
                int countLink = 0;

                std::vector<int> indRayCorelist = findMatchingIndices(UEswithSat.uvSet, UEswithSat.uvSetCore);
                std::cout << "\n Count indUvCorelist: " << indRayCorelist.size() << std::endl;
                std::vector<int>  RayUElist;
                std::vector<int> coreUElist;


                Eigen::Matrix3d rotMatrix = UEswithSat.rotMatrix;

                for (auto& link : UEswithSat.links.getLinks())
                {

                    double deg = link.elevationAngle;

                    int index = int(AngleForLSP(deg)) / 10;

                    link.isLos = CalculateLOSProbability(index, scenario);


                    MatrixXd Table = link.isLos ? TableLOS : TableNLOS;
                    VectorXd Parameters = Table.col(index - 1);

                    LSP::initializeParameters(link.isLos, link, Parameters);
                    std::cout << link.isLos << "\n";
                    SSP::setParameters(Parameters);
                    SSP::calculateLosAngles(link, UEswithSat.p1, UEswithSat.p2);
                    SSP::generateClusterDelays(link);
                    std::pair<std::vector<double>, std::vector<double>> clusterPower = SSP::generateClusterPowers(link);

                    //Вывод мощностей ( без / с масштабирующим фактором) 
                    //std::cout << "clusterPower: " << clusterPower.first[0] << ", " << clusterPower.first[1] << ", " << clusterPower.first[2] << "\n";
                    //std::cout << "clusterPowersWithScalingFactors: " << clusterPower.second[0] << ", " << clusterPower.second[1] << ", " << clusterPower.second[2] << "\n";


                    /////////////////////////////////////////Нужно позже оформить ...  ниже почти всё для расчёта бюджета линии////////////////////////////////////

                    // для PL 
                    double d = CalculateDistance(EARTH_RADIUS, altitude, deg * PI / 180);
                    double std = ChooseSTD(link.isLos, f, deg, scenario);
                    double SF = GenerateSF(std);
                    double FSPL = Calculate_FSPL(d * 1e3, f);
                    double CL = ChooseCL(link.isLos, f, deg, scenario);
                    double PL_dB = CalculateBasisPathLoss(FSPL, SF, CL);


                    // генератор для шума
                    std::random_device rd;  // Источник энтропии
                    std::mt19937 gen(rd()); // Генератор случайных чисел (Mersenne Twister)
                    std::normal_distribution<> dist(0, 0.72); // Нормальное распределение с mean=0 и stddev=0.72
                    PL_dB = PL_dB + dist(gen); // добавляем SF

                    // Перевод в линейный масштаб
                    double PL_lin_magnitude = pow(10, -1 * PL_dB / 20);
                    //std::cout << "\nPL_lin_magnitude: " << PL_lin_magnitude << std::endl;



                    for (int i = 0; i < static_cast<int> (rRays.size()); ++i)
                    {

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

                    countLink += 1;
                }

                // Вывод задержек во времени
                for (auto& link : UEswithSat.links.getLinks()) {
                    std::cout << "Cluster Delay: " << link.clusterDelays[0] << ", " << link.clusterDelays[1] << ", " << link.clusterDelays[2] << ", ";
                    if (!link.isLos) { std::cout << link.clusterDelays[3]; }
                    std::cout << " \n";
                }


                // Заполняем coreUElist

                for (int rayIndex : RayUElist) {
                    bool found = false;
                    for (int coreIndex : indRayCorelist) {
                        if (coreIndex == rayIndex) {
                            found = true;
                            break;
                        }
                    }
                    coreUElist.push_back(found ? 1 : 0); // Добавляем 1, если элемент найден, иначе 0
                }

                std::vector<double> coreArrPatt = filterVector(arrPatt, coreUElist);
                std::vector<double> maxPL_lin_magn = filterVector(PL_lin_magn, coreUElist);
                MatrixXd allArrPatt_magn = filterMatrixRows(arrPatt_magn, coreUElist);

                LinkBudget linkBudget(30, 7, 290); //  BW_MHz = 30 МГц, NF_dB = 7 дБ, T = 290 К
                linkBudget.calculateMetrics(coreArrPatt, maxPL_lin_magn, allArrPatt_magn);
                plotter.plotCDF(linkBudget.getDL_CNR_dB(), linkBudget.getDL_CIR_dB(), linkBudget.getDL_CINR_dB(), scenario);

            }
            return 0;

            //plotter.plotRayPoints(UEswithSat.links.getLinks()[0].satellitePosition, rRays);


            // Вычисление диаграммы направленности
            //calculateDishPattern(elTargetDegrees, AP_dB, rDish_WL);
        }


        // mode 2 - построенеи графиков для LSP

 else
 {
     // 1) разделить на Los и NLOS

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

     double f = (frequencyChoice == 1) ? 2.0 : 30.0; // Example: 2 GHz for S and 30 GHz for Ka
     double altitude = 600.0;

     int z = 1;

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
             double elTargetDegrees = elevation;
             double azTargetDegrees = RandomGenerators::generateUniform(0.0, 360.0);

             SatelliteLink UEswithSat(2, level, Nueincell, altitude, 10.0, elTargetDegrees, azTargetDegrees);
             UEswithSat.generateLinks();

             MatrixXd TableLOS = GenerateMatrix(true, f, scenarios[0]);
             MatrixXd TableNLOS = GenerateMatrix(false, f, scenarios[0]);

             for (auto& link : UEswithSat.links.getLinks()) {
                 double deg = link.elevationAngle;
                 int index = int(AngleForLSP(deg)) / 10;
                 bool isLos = CalculateLOSProbability(index, scenarios[0]);

                 MatrixXd Table = isLos ? TableLOS : TableNLOS;
                 VectorXd Parameters = Table.col(index - 1);
                 LSP::initializeParameters(isLos, link, Parameters);

                 if (isLos)
                 {
                     SFDU.push_back(link.SF_db);
                     KDU.push_back(link.K_db);
                     DSDU.push_back(pow(10, link.DS_sec));
                     ASDDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                     ASADU.push_back(pow(10, link.ASA_deg));
                     ZSDDU.push_back(pow(10, link.ZSD_deg));
                     ZSADU.push_back(pow(10, link.ZSA_deg));
                 }
                 else
                 {
                     SFDU.push_back(link.SF_db);
                     DSDU.push_back(pow(10, link.DS_sec));
                     ASDDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                     ASADU.push_back(pow(10, link.ASA_deg));
                     ZSDDU.push_back(pow(10, link.ZSD_deg));
                     ZSADU.push_back(pow(10, link.ZSA_deg));
                 }

             }

         }

         std::sort(SFDU.data(), SFDU.data() + SFDU.size());
         std::sort(KDU.data(), KDU.data() + KDU.size());
         std::sort(DSDU.data(), DSDU.data() + DSDU.size());
         std::sort(ASDDU.data(), ASDDU.data() + ASDDU.size());
         std::sort(ASADU.data(), ASADU.data() + ASADU.size());
         std::sort(ZSDDU.data(), ZSDDU.data() + ZSDDU.size());
         std::sort(ZSADU.data(), ZSADU.data() + ZSADU.size());

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
                     bool isLos = CalculateLOSProbability(index, scenarios[1]);

                     MatrixXd Table = isLos ? TableLOS : TableNLOS;
                     VectorXd Parameters = Table.col(index - 1);
                     LSP::initializeParameters(isLos, link, Parameters);

                     if (isLos)
                     {
                         SFU.push_back(link.SF_db);
                         KU.push_back(link.K_db);
                         DSU.push_back(pow(10, link.DS_sec));
                         ASDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                         ASAU.push_back(pow(10, link.ASA_deg));
                         ZSDU.push_back(pow(10, link.ZSD_deg));
                         ZSAU.push_back(pow(10, link.ZSA_deg));
                     }
                     else
                     {
                         SFU.push_back(link.SF_db);
                         DSU.push_back(pow(10, link.DS_sec));
                         ASDU.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                         ASAU.push_back(pow(10, link.ASA_deg));
                         ZSDU.push_back(pow(10, link.ZSD_deg));
                         ZSAU.push_back(pow(10, link.ZSA_deg));
                     }

                 }

             }

             std::sort(SFU.data(), SFU.data() + SFU.size());
             std::sort(KU.data(), KU.data() + KU.size());
             std::sort(DSU.data(), DSU.data() + DSU.size());
             std::sort(ASDU.data(), ASDU.data() + ASDU.size());
             std::sort(ASAU.data(), ASAU.data() + ASAU.size());
             std::sort(ZSDU.data(), ZSDU.data() + ZSDU.size());
             std::sort(ZSAU.data(), ZSAU.data() + ZSAU.size());

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
                     bool isLos = CalculateLOSProbability(index, scenarios[2]);

                     MatrixXd Table = isLos ? TableLOS : TableNLOS;
                     VectorXd Parameters = Table.col(index - 1);
                     LSP::initializeParameters(isLos, link, Parameters);

                     if (isLos)
                     {
                         SFS.push_back(link.SF_db);
                         KS.push_back(link.K_db);
                         DSS.push_back(pow(10, link.DS_sec));
                         ASDS.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                         ASAS.push_back(pow(10, link.ASA_deg));
                         ZSDS.push_back(pow(10, link.ZSD_deg));
                         ZSAS.push_back(pow(10, link.ZSA_deg));
                     }
                     else
                     {
                         SFS.push_back(link.SF_db);
                         DSS.push_back(pow(10, link.DS_sec));
                         ASDS.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                         ASAS.push_back(pow(10, link.ASA_deg));
                         ZSDS.push_back(pow(10, link.ZSD_deg));
                         ZSAS.push_back(pow(10, link.ZSA_deg));
                     }

                 }

             }

             std::sort(SFS.data(), SFS.data() + SFS.size());
             std::sort(KS.data(), KS.data() + KS.size());
             std::sort(DSS.data(), DSS.data() + DSS.size());
             std::sort(ASDS.data(), ASDS.data() + ASDS.size());
             std::sort(ASAS.data(), ASAS.data() + ASAS.size());
             std::sort(ZSDS.data(), ZSDS.data() + ZSDS.size());
             std::sort(ZSAS.data(), ZSAS.data() + ZSAS.size());
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
                     bool isLos = CalculateLOSProbability(index, scenarios[3]);

                     MatrixXd Table = isLos ? TableLOS : TableNLOS;
                     VectorXd Parameters = Table.col(index - 1);
                     LSP::initializeParameters(isLos, link, Parameters);

                     if (isLos)
                     {
                         SFR.push_back(link.SF_db);
                         KR.push_back(link.K_db);
                         DSR.push_back(pow(10, link.DS_sec));
                         ASDR.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                         ASAR.push_back(pow(10, link.ASA_deg));
                         ZSDR.push_back(pow(10, link.ZSD_deg));
                         ZSAR.push_back(pow(10, link.ZSA_deg));
                     }
                     else
                     {
                         SFR.push_back(link.SF_db);
                         DSR.push_back(pow(10, link.DS_sec));
                         ASDR.push_back(M_PI / 180 * pow(10, std::min(link.ASD_deg, log10(90))));
                         ASAR.push_back(pow(10, link.ASA_deg));
                         ZSDR.push_back(pow(10, link.ZSD_deg));
                         ZSAR.push_back(pow(10, link.ZSA_deg));
                     }

                 }

             }

             std::sort(SFR.data(), SFR.data() + SFR.size());
             std::sort(KR.data(), KR.data() + KR.size());
             std::sort(DSR.data(), DSR.data() + DSR.size());
             std::sort(ASDR.data(), ASDR.data() + ASDR.size());
             std::sort(ASAR.data(), ASAR.data() + ASAR.size());
             std::sort(ZSDR.data(), ZSDR.data() + ZSDR.size());
             std::sort(ZSAR.data(), ZSAR.data() + ZSAR.size());

         }
     }

     MatlabLSPPlot plotter;
     plotter.plotForAllLSP(SFDU, SFU, SFS, SFR, KDU, KU, KS, KR, DSDU, DSU, DSS, DSR, ASDDU, ASDU, ASDS, ASDR, ASADU, ASAU, ASAS, ASAR, ZSDDU, ZSDU, ZSDS, ZSDR, ZSADU, ZSAU, ZSAS, ZSAR, scenarios, frequencyBands[frequencyChoice - 1]);


     std::cin;
    }

    }

}
