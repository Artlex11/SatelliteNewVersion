#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <chrono>

#include "NTN_Deployment.h"
#include "Tables.h"
#include "Matlab_plot.h"
#include "LOS_Probability.h"
#include "Pathloss.h"

int main() 
{

    std::string scenario;
    double f;
    bool isLos;
    double altitude;

    std::cout << "\nInput scenario (Dense_Urban/Urban/Suburban/Rural): ";
    std::cin >> scenario;
    std::cout << "\nInput frequency(f) in GHz: ";
    std::cin >> f;
    std::cout << "\nInput altitude in kilometers: ";
    std::cin >> altitude;
    auto start = std::chrono::high_resolution_clock::now();

    SatelliteLink UEswithSat(2, 4, 10, altitude, 10, 90);
    UEswithSat.generateLinks();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    //std::cout << "Elapsed time: " << elapsed.count() << " s\n";

    std::vector<Eigen::Vector3d> users;
    for (const auto& link : UEswithSat.links.getLinks()) {
        users.push_back(link.userPosition);
    }


    // Создание объекта для построения графиков
    MatlabPlot plotter;

    // Построение графиков
    plotter.plotTransformedData(users, UEswithSat.links.getLinks()[0].satellitePosition);
    plotter.plotEarth();

    MatrixXd TableLOS = GenerateMatrix(true, f, scenario);
    MatrixXd TableNLOS = GenerateMatrix(false, f, scenario);

    for (auto& link : UEswithSat.links.getLinks()) {
        double deg = link.elevationAngle;

        int index = int(AngleForLSP(deg)) / 10;
        std::cout << "\nUser " << ": index: " << index << ", angle: " << AngleForLSP(deg) << "\n";
        isLos = CalculateLOSProbability(index, scenario);
        link.isLOS = isLos;
        double d = CalculateDistance(EARTH_RADIUS, altitude, deg * PI / 180);

        std::cout << "Distance: " << d;

        double std = ChooseSTD(isLos, f, deg, scenario);
        //std::cout << "\nSTD for SF: " << std;

        double SF = GenerateSF(std);
        //std::cout << "\nSF: " << SF;

        double FSPL = Calculate_FSPL(d, f);
        //std::cout << "\nFSPL: " << FSPL;

        double CL = ChooseCL(isLos, f, deg, scenario);
        //std::cout << "\nCL: " << CL;

        double PL_b = CalculateBasisPathLoss(FSPL, SF, CL);
        std::cout << "\nPL_b: " << PL_b << std::endl;

        setlocale(LC_ALL, "RU");

        // Расчёт газовых потерь 
        Engine* ep;
        MATFile* pmat;
        mxArray* pa_xData;
        mxArray* pa_yData;
        double* xData;
        double* yData;
        int numElements;

        // Запускаем MATLAB Engine
        ep = engOpen(NULL);
        if (ep == NULL) {
            std::cerr << "Не удалось запустить MATLAB Engine!" << std::endl;
            return 1;
        }

        // Открываем .mat файл
        pmat = matOpen("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 2 (СПУТНИК)\\ChannelModel_NTN\\LSP\\gasAttTable.mat", "r");
        if (pmat == NULL) {
            std::cerr << "Не удалось открыть .mat файл!" << std::endl;
            engClose(ep);
            return 1;
        }

        // Читаем xData из .mat файла
        pa_xData = matGetVariable(pmat, "xData");
        if (pa_xData == NULL) {
            std::cerr << "Не удалось прочитать xData из .mat файла!" << std::endl;
            matClose(pmat);
            engClose(ep);
            return 1;
        }

        // Читаем yData из .mat файла
        pa_yData = matGetVariable(pmat, "yData");
        if (pa_yData == NULL) {
            std::cerr << "Не удалось прочитать yData из .mat файла!" << std::endl;
            mxDestroyArray(pa_xData); // Освобождаем память, выделенную для xData
            matClose(pmat);
            engClose(ep);
            return 1;
        }

        // Получаем данные из mxArray для xData
        numElements = mxGetNumberOfElements(pa_xData); // Должно быть 1000
        xData = mxGetPr(pa_xData);

        // Получаем данные из mxArray для yData
        int numElementsY = mxGetNumberOfElements(pa_yData);  // Должно быть 1000
        yData = mxGetPr(pa_yData);

        // Проверка, что оба массива имеют одинаковую длину (опционально)
        if (numElements != numElementsY) {
            std::cerr << "Внимание: массивы xData и yData имеют разную длину!" << std::endl;
        }

        // Выводим несколько элементов для проверки

        /*std::cout << "Первые 5 элементов xData:" << std::endl;
        for (int i = 0; i < 5; ++i) {
            std::cout << xData[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Первые 5 элементов yData:" << std::endl;
        for (int i = 0; i < 5; ++i) {
            std::cout << yData[i] << " ";
        }
        std::cout << std::endl;*/

        double logFC = log10(f);
        //std::cout << "logFC: " << logFC << std::endl;
        int ind; // Переменная, которая ищет элемент из mat таблицы

        for (int i = 0; i < 1000; ++i)
        {
            if (logFC > xData[i])
            {
                continue;
            }
            else
            {
                ind = i;
                break;
            }

        }

        //std::cout << "Index of element: " << ind << std::endl;
        //std::cout << "xData[index]: " << xData[ind] << std::endl;
        //std::cout << "yData[index]: " << yData[ind] << std::endl;

        double PL_g;

        PL_g = pow(10, yData[ind]) * d;

        std::cout << "PL_g: " << PL_g << std::endl;
        std::cout << "PL = PL_b + PL_g = " << PL_b + PL_g << std::endl;

        // Закрываем .mat файл и MATLAB Engine
        mxDestroyArray(pa_xData);
        mxDestroyArray(pa_yData);
        matClose(pmat);
        engClose(ep);

        MatrixXd Table = isLos ? TableLOS : TableNLOS;
        VectorXd Parameters{ Table.rows() };
        Parameters = Table.col(index - 1);
        // расстояние в линке 
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