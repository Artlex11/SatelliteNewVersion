#include "Matlab_plot.h"
#include <iostream>

MatlabPlot::MatlabPlot() {
    if (!(ep = engOpen(""))) {
        std::cerr << "Не удалось открыть MATLAB Engine" << std::endl;
        exit(1);
    }
}

MatlabPlot::~MatlabPlot() {
    engClose(ep);
}

void MatlabPlot::plotTransformedData(const std::vector<Eigen::Vector3d>& users, const Eigen::Vector3d& satellite) {
    // Передача данных пользователей в MATLAB
    mxArray* usersX = mxCreateDoubleMatrix(1, users.size(), mxREAL);
    mxArray* usersY = mxCreateDoubleMatrix(1, users.size(), mxREAL);
    mxArray* usersZ = mxCreateDoubleMatrix(1, users.size(), mxREAL);

    double* px = mxGetPr(usersX);
    double* py = mxGetPr(usersY);
    double* pz = mxGetPr(usersZ);

    for (size_t i = 0; i < users.size(); ++i) {
        px[i] = users[i].x();
        py[i] = users[i].y();
        pz[i] = users[i].z();
    }

    engPutVariable(ep, "usersX", usersX);
    engPutVariable(ep, "usersY", usersY);
    engPutVariable(ep, "usersZ", usersZ);

    // Передача данных спутника в MATLAB
    mxArray* satX = mxCreateDoubleScalar(satellite.x());
    mxArray* satY = mxCreateDoubleScalar(satellite.y());
    mxArray* satZ = mxCreateDoubleScalar(satellite.z());

    engPutVariable(ep, "satX", satX);
    engPutVariable(ep, "satY", satY);
    engPutVariable(ep, "satZ", satZ);

    // Построение графиков в MATLAB
    engEvalString(ep, "hold on; grid on;"); // Добавляем новые данные на текущий график
    engEvalString(ep, "plot3(usersX, usersY, usersZ, '.', 'DisplayName', 'Координаты пользователей');");
    engEvalString(ep, "plot3(satX, satY, satZ, '.', 'MarkerSize', 10, 'DisplayName', 'Спутник');");

    engEvalString(ep, "plot3(0,0,0,'r.'); axis equal;xlabel('x, km'); ylabel('y, km'); zlabel('z, km'); ");

    // Освобождение памяти
    mxDestroyArray(usersX);
    mxDestroyArray(usersY);
    mxDestroyArray(usersZ);
    mxDestroyArray(satX);
    mxDestroyArray(satY);
    mxDestroyArray(satZ);

}

void MatlabPlot::plotRayPoints(const Eigen::Vector3d& satellitePosition, const std::vector<Eigen::Vector3d>& rays) {
    // Проверка данных
    std::cout << "Satellite position: " << satellitePosition.transpose() << std::endl;
    for (size_t i = 0; i < rays.size(); ++i) {
        std::cout << "Ray " << i + 1 << ": " << rays[i].transpose() << std::endl;
    }

    // Проверка на наличие лучей
    if (rays.empty()) {
        std::cerr << "Ошибка: Массив лучей пуст." << std::endl;
        return;
    }

    // Вычисление конечных точек лучей
    std::vector<Eigen::Vector3d> rayPoints;
    for (const auto& ray : rays) {
        rayPoints.push_back(satellitePosition + ray);
    }

    // Передача данных в MATLAB
    mxArray* satPosArray = mxCreateDoubleMatrix(1, 3, mxREAL);
    std::memcpy(mxGetPr(satPosArray), satellitePosition.data(), 3 * sizeof(double));
    engPutVariable(ep, "satPos", satPosArray);

    // Создание массива для конечных точек лучей
    mxArray* rayPointsArray = mxCreateDoubleMatrix(rayPoints.size(), 3, mxREAL);
    double* rayPointsData = mxGetPr(rayPointsArray);

    // Заполнение массива данными
    for (size_t i = 0; i < rayPoints.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rayPointsData[i + j * rayPoints.size()] = rayPoints[i](j);  // Транспонирование данных
        }
    }
    engPutVariable(ep, "rayPoints", rayPointsArray);

    // Отрисовка спутника
    engEvalString(ep, " hold on;");

    // Отрисовка конечных точек лучей
    engEvalString(ep, "plot3(satPos(1), satPos(2), satPos(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Спутник');");

    engEvalString(ep, "plot3(rayPoints(:, 1), rayPoints(:, 2), rayPoints(:, 3),'.','DisplayName', 'Точки лучей');");

    // Настройка графика
    engEvalString(ep, "grid on;");
    engEvalString(ep, "xlabel('X'); ylabel('Y'); zlabel('Z');");
    engEvalString(ep, "title('Точки лучей от спутника'); legend show;");
    engEvalString(ep, "hold off;");

    // Освобождение ресурсов
    mxDestroyArray(satPosArray);
    mxDestroyArray(rayPointsArray);
}


void MatlabPlot::plotEarth() 
{
    engEvalString(ep, "hold on;");
    engEvalString(ep, "[phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));");
    engEvalString(ep, "R = 6371;");
    engEvalString(ep, "x_sphere = R * cos(phi) .* cos(theta);");
    engEvalString(ep, "y_sphere = R * sin(phi) .* cos(theta);");
    engEvalString(ep, "z_sphere = R * sin(theta);");
    engEvalString(ep, "surf(x_sphere, y_sphere, z_sphere, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Сфера Земли');");
    engEvalString(ep, "axis equal;");
    engEvalString(ep, "xlabel('x, km'); ylabel('y, km'); zlabel('z, km');");
    engEvalString(ep, "grid on");
    engEvalString(ep, "legend show;");
}

void MatlabPlot::plotCDF(const std::vector<double>& DL_CNR_dB_vec, const std::vector<double>& DL_CIR_dB_vec, const std::vector<double>& DL_CINR_dB_vec, const std::string& scenario) {
    // Проверка на пустые векторы
    if (DL_CNR_dB_vec.empty() || DL_CIR_dB_vec.empty() || DL_CINR_dB_vec.empty()) {
        std::cerr << "Error: One or more vectors are empty." << std::endl;
        return;
    }

    // Передача данных CNR в MATLAB
    mxArray* cnrArray = mxCreateDoubleMatrix(1, DL_CNR_dB_vec.size(), mxREAL);
    std::memcpy(mxGetPr(cnrArray), DL_CNR_dB_vec.data(), DL_CNR_dB_vec.size() * sizeof(double));
    engPutVariable(ep, "cnr", cnrArray);

    // Передача данных CIR в MATLAB
    mxArray* cirArray = mxCreateDoubleMatrix(1, DL_CIR_dB_vec.size(), mxREAL);
    std::memcpy(mxGetPr(cirArray), DL_CIR_dB_vec.data(), DL_CIR_dB_vec.size() * sizeof(double));
    engPutVariable(ep, "cir", cirArray);

    // Передача данных CINR в MATLAB
    mxArray* cinrArray = mxCreateDoubleMatrix(1, DL_CINR_dB_vec.size(), mxREAL);
    std::memcpy(mxGetPr(cinrArray), DL_CINR_dB_vec.data(), DL_CINR_dB_vec.size() * sizeof(double));
    engPutVariable(ep, "cinr", cinrArray);

    // Передача сценария в MATLAB
    mxArray* scenarioArray = mxCreateString(scenario.c_str());
    engPutVariable(ep, "scenario", scenarioArray);

    // Построение CDF в MATLAB
    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_cnr, x_cnr] = ecdf(cnr); plot(x_cnr, f_cnr, 'DisplayName', 'CNR');");
    engEvalString(ep, "[f_cir, x_cir] = ecdf(cir); plot(x_cir, f_cir, 'DisplayName', 'CIR');");
    engEvalString(ep, "[f_cinr, x_cinr] = ecdf(cinr); plot(x_cinr, f_cinr, 'DisplayName', 'CINR');");
    engEvalString(ep, "xlabel('SNR, дБ'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF для CNR, CIR и CINR (Сценарий: ', scenario, ')']);");
    engEvalString(ep, "legend show;");

    // Освобождение памяти
    mxDestroyArray(cnrArray);
    mxDestroyArray(cirArray);
    mxDestroyArray(cinrArray);
    mxDestroyArray(scenarioArray);
}

// LSP Graphics

MatlabLSPPlot::MatlabLSPPlot() {
    if (!(ep = engOpen(""))) {
        std::cerr << "Не удалось открыть MATLAB Engine" << std::endl;
        exit(1);
    }
}

MatlabLSPPlot::~MatlabLSPPlot() {
    engClose(ep);
}

void MatlabLSPPlot::plotForAllLSP(std::vector<double> SFDU, std::vector<double> SFU, std::vector<double> SFS, std::vector<double> SFR, std::vector<double> KDU, std::vector<double> KU, std::vector<double> KS, std::vector<double> KR, std::vector<double> DSDU, std::vector<double> DSU, std::vector<double> DSS, std::vector<double> DSR, std::vector<double> ASDDU, std::vector<double> ASDU, std::vector<double> ASDS, std::vector<double> ASDR, std::vector<double> ASADU, std::vector<double> ASAU, std::vector<double> ASAS, std::vector<double> ASAR, std::vector<double> ZSDDU, std::vector<double> ZSDU, std::vector<double> ZSDS, std::vector<double> ZSDR, std::vector<double> ZSADU, std::vector<double> ZSAU, std::vector<double> ZSAS, std::vector<double> ZSAR, std::vector<std::string> scenarios, std::string frequencyBand)
{
    // Передача сценария в MATLAB
    mxArray* scenarioArrayDU = mxCreateString(scenarios[0].c_str());
    mxArray* scenarioArrayU = mxCreateString(scenarios[1].c_str());
    mxArray* scenarioArrayS = mxCreateString(scenarios[2].c_str());
    mxArray* scenarioArrayR = mxCreateString(scenarios[3].c_str());
    mxArray* BandsArray = mxCreateString(frequencyBand.c_str());
    engPutVariable(ep, "scenarioDU", scenarioArrayDU);
    engPutVariable(ep, "scenarioU", scenarioArrayU);
    engPutVariable(ep, "scenarioS", scenarioArrayS);
    engPutVariable(ep, "scenarioR", scenarioArrayR);
    engPutVariable(ep, "frequencyBand", BandsArray);

    // Передача данных SF в MATLAB

    mxArray* SFArrayDU = mxCreateDoubleMatrix(1, SFDU.size(), mxREAL);
    std::memcpy(mxGetPr(SFArrayDU), SFDU.data(), SFDU.size() * sizeof(double));
    engPutVariable(ep, "SFDU", SFArrayDU);

    mxArray* SFArrayU = mxCreateDoubleMatrix(1, SFU.size(), mxREAL);
    std::memcpy(mxGetPr(SFArrayU), SFU.data(), SFU.size() * sizeof(double));
    engPutVariable(ep, "SFU", SFArrayU);

    mxArray* SFArrayS = mxCreateDoubleMatrix(1, SFS.size(), mxREAL);
    std::memcpy(mxGetPr(SFArrayS), SFS.data(), SFS.size() * sizeof(double));
    engPutVariable(ep, "SFS", SFArrayS);

    mxArray* SFArrayR = mxCreateDoubleMatrix(1, SFR.size(), mxREAL);
    std::memcpy(mxGetPr(SFArrayR), SFR.data(), SFR.size() * sizeof(double));
    engPutVariable(ep, "SFR", SFArrayR);

    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_SF, x_SF] = ecdf(SFDU); plot(x_SF, f_SF, 'DisplayName', scenarioDU);");   
    engEvalString(ep, "[f_SF, x_SF] = ecdf(SFU); plot(x_SF, f_SF, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_SF, x_SF] = ecdf(SFS); plot(x_SF, f_SF, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_SF, x_SF] = ecdf(SFR); plot(x_SF, f_SF, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('SF, дБ'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF SF for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");

    //// Передача данных K в MATLAB
    mxArray* KArrayDU = mxCreateDoubleMatrix(1, KDU.size(), mxREAL);
    std::memcpy(mxGetPr(KArrayDU), KDU.data(), KDU.size() * sizeof(double));
    engPutVariable(ep, "KDU", KArrayDU);
    
    mxArray* KArrayU = mxCreateDoubleMatrix(1, KU.size(), mxREAL);
    std::memcpy(mxGetPr(KArrayU), KU.data(), KU.size() * sizeof(double));
    engPutVariable(ep, "KU", KArrayU);

    mxArray* KArrayS = mxCreateDoubleMatrix(1, KS.size(), mxREAL);
    std::memcpy(mxGetPr(KArrayS), KS.data(), KS.size() * sizeof(double));
    engPutVariable(ep, "KS", KArrayS);

    mxArray* KArrayR = mxCreateDoubleMatrix(1, KR.size(), mxREAL);
    std::memcpy(mxGetPr(KArrayR), KR.data(), KR.size() * sizeof(double));
    engPutVariable(ep, "KR", KArrayR);

    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_K, x_K] = ecdf(KDU); plot(x_K, f_K, 'DisplayName', scenarioDU);");
    engEvalString(ep, "[f_K, x_K] = ecdf(KU); plot(x_K, f_K, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_K, x_K] = ecdf(KS); plot(x_K, f_K, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_K, x_K] = ecdf(KR); plot(x_K, f_K, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('K, дБ'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF K for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");

    // Передача данных DS в MATLAB
    mxArray* DSArrayDU = mxCreateDoubleMatrix(1, DSDU.size(), mxREAL);
    std::memcpy(mxGetPr(DSArrayDU), DSDU.data(), DSDU.size() * sizeof(double));
    engPutVariable(ep, "DSDU", DSArrayDU);

    mxArray* DSArrayU = mxCreateDoubleMatrix(1, DSU.size(), mxREAL);
    std::memcpy(mxGetPr(DSArrayU), DSU.data(), DSU.size() * sizeof(double));
    engPutVariable(ep, "DSU", DSArrayU);

    mxArray* DSArrayS = mxCreateDoubleMatrix(1, DSS.size(), mxREAL);
    std::memcpy(mxGetPr(DSArrayS), DSS.data(), DSS.size() * sizeof(double));
    engPutVariable(ep, "DSS", DSArrayS);

    mxArray* DSArrayR = mxCreateDoubleMatrix(1, DSR.size(), mxREAL);
    std::memcpy(mxGetPr(DSArrayR), DSR.data(), DSR.size() * sizeof(double));
    engPutVariable(ep, "DSR", DSArrayR);


    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_DS, x_DS] = ecdf(DSDU); plot(x_DS, f_DS, 'DisplayName', scenarioDU);");
    engEvalString(ep, "[f_DS, x_DS] = ecdf(DSU); plot(x_DS, f_DS, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_DS, x_DS] = ecdf(DSS); plot(x_DS, f_DS, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_DS, x_DS] = ecdf(DSR); plot(x_DS, f_DS, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('DS, сек'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF DS for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");

    //// Передача данных ASD в MATLAB
    mxArray* ASDArrayDU = mxCreateDoubleMatrix(1, ASDDU.size(), mxREAL);
    std::memcpy(mxGetPr(ASDArrayDU), ASDDU.data(), ASDDU.size() * sizeof(double));
    engPutVariable(ep, "ASDDU", ASDArrayDU);

    mxArray* ASDArrayU = mxCreateDoubleMatrix(1, ASDU.size(), mxREAL);
    std::memcpy(mxGetPr(ASDArrayU), ASDU.data(), ASDU.size() * sizeof(double));
    engPutVariable(ep, "ASDU", ASDArrayU);

    mxArray* ASDArrayS = mxCreateDoubleMatrix(1, ASDS.size(), mxREAL);
    std::memcpy(mxGetPr(ASDArrayS), ASDS.data(), ASDS.size() * sizeof(double));
    engPutVariable(ep, "ASDS", ASDArrayS);

    mxArray* ASDArrayR = mxCreateDoubleMatrix(1, ASDR.size(), mxREAL);
    std::memcpy(mxGetPr(ASDArrayR), ASDR.data(), ASDR.size() * sizeof(double));
    engPutVariable(ep, "ASDR", ASDArrayR);


    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_ASD, x_ASD] = ecdf(ASDDU); plot(x_ASD, f_ASD, 'DisplayName', scenarioDU);");
    engEvalString(ep, "[f_ASD, x_ASD] = ecdf(ASDU); plot(x_ASD, f_ASD, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_ASD, x_ASD] = ecdf(ASDS); plot(x_ASD, f_ASD, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_ASD, x_ASD] = ecdf(ASDR); plot(x_ASD, f_ASD, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('ASD, rad'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF ASD for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");

    //// Передача данных ASA в MATLAB
    mxArray* ASAArrayDU = mxCreateDoubleMatrix(1, ASADU.size(), mxREAL);
    std::memcpy(mxGetPr(ASAArrayDU), ASADU.data(), ASADU.size() * sizeof(double));
    engPutVariable(ep, "ASADU", ASAArrayDU);

    mxArray* ASAArrayU = mxCreateDoubleMatrix(1, ASAU.size(), mxREAL);
    std::memcpy(mxGetPr(ASAArrayU), ASAU.data(), ASAU.size() * sizeof(double));
    engPutVariable(ep, "ASAU", ASAArrayU);

    mxArray* ASAArrayS = mxCreateDoubleMatrix(1, ASAS.size(), mxREAL);
    std::memcpy(mxGetPr(ASAArrayS), ASAS.data(), ASAS.size() * sizeof(double));
    engPutVariable(ep, "ASAS", ASAArrayS);

    mxArray* ASAArrayR = mxCreateDoubleMatrix(1, ASAR.size(), mxREAL);
    std::memcpy(mxGetPr(ASAArrayR), ASAR.data(), ASAR.size() * sizeof(double));
    engPutVariable(ep, "ASAR", ASAArrayR);

    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_ASA, x_ASA] = ecdf(ASADU); plot(x_ASA, f_ASA, 'DisplayName', scenarioDU);");
    engEvalString(ep, "[f_ASA, x_ASA] = ecdf(ASAU); plot(x_ASA, f_ASA, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_ASA, x_ASA] = ecdf(ASAS); plot(x_ASA, f_ASA, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_ASA, x_ASA] = ecdf(ASAR); plot(x_ASA, f_ASA, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('ASA, rad'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF ASA for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");

    //// Передача данных ZSD в MATLAB
    mxArray* ZSDArrayDU = mxCreateDoubleMatrix(1, ZSDDU.size(), mxREAL);
    std::memcpy(mxGetPr(ZSDArrayDU), ZSDDU.data(), ZSDDU.size() * sizeof(double));
    engPutVariable(ep, "ZSDDU", ZSDArrayDU);

    mxArray* ZSDArrayU = mxCreateDoubleMatrix(1, ZSDU.size(), mxREAL);
    std::memcpy(mxGetPr(ZSDArrayU), ZSDU.data(), ZSDU.size() * sizeof(double));
    engPutVariable(ep, "ZSDU", ZSDArrayU);

    mxArray* ZSDArrayS = mxCreateDoubleMatrix(1, ZSDS.size(), mxREAL);
    std::memcpy(mxGetPr(ZSDArrayS), ZSDS.data(), ZSDS.size() * sizeof(double));
    engPutVariable(ep, "ZSDS", ZSDArrayS);

    mxArray* ZSDArrayR = mxCreateDoubleMatrix(1, ZSDR.size(), mxREAL);
    std::memcpy(mxGetPr(ZSDArrayR), ZSDR.data(), ZSDR.size() * sizeof(double));
    engPutVariable(ep, "ZSDR", ZSDArrayR);

    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_ZSD, x_ZSD] = ecdf(ZSDDU); plot(x_ZSD, f_ZSD, 'DisplayName', scenarioDU);");
    engEvalString(ep, "[f_ZSD, x_ZSD] = ecdf(ZSDU); plot(x_ZSD, f_ZSD, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_ZSD, x_ZSD] = ecdf(ZSDS); plot(x_ZSD, f_ZSD, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_ZSD, x_ZSD] = ecdf(ZSDR); plot(x_ZSD, f_ZSD, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('ZSD, rad'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF ZSD for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");

    //// Передача данных ZSA в MATLAB
    mxArray* ZSAArrayDU = mxCreateDoubleMatrix(1, ZSADU.size(), mxREAL);
    std::memcpy(mxGetPr(ZSAArrayDU), ZSADU.data(), ZSADU.size() * sizeof(double));
    engPutVariable(ep, "ZSADU", ZSAArrayDU);

    mxArray* ZSAArrayU = mxCreateDoubleMatrix(1, ZSAU.size(), mxREAL);
    std::memcpy(mxGetPr(ZSAArrayU), ZSAU.data(), ZSAU.size() * sizeof(double));
    engPutVariable(ep, "ZSAU", ZSAArrayU);

    mxArray* ZSAArrayS = mxCreateDoubleMatrix(1, ZSAS.size(), mxREAL);
    std::memcpy(mxGetPr(ZSAArrayS), ZSAS.data(), ZSAS.size() * sizeof(double));
    engPutVariable(ep, "ZSAS", ZSAArrayS);

    mxArray* ZSAArrayR = mxCreateDoubleMatrix(1, ZSAR.size(), mxREAL);
    std::memcpy(mxGetPr(ZSAArrayR), ZSAR.data(), ZSAR.size() * sizeof(double));
    engPutVariable(ep, "ZSAR", ZSAArrayR);


    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[f_ZSA, x_ZSA] = ecdf(ZSADU); plot(x_ZSA, f_ZSA, 'DisplayName', scenarioDU);");
    engEvalString(ep, "[f_ZSA, x_ZSA] = ecdf(ZSAU); plot(x_ZSA, f_ZSA, 'DisplayName', scenarioU);");
    engEvalString(ep, "[f_ZSA, x_ZSA] = ecdf(ZSAS); plot(x_ZSA, f_ZSA, 'DisplayName', scenarioS);");
    engEvalString(ep, "[f_ZSA, x_ZSA] = ecdf(ZSAR); plot(x_ZSA, f_ZSA, 'DisplayName', scenarioR);");
    engEvalString(ep, "xlabel('ZSA, rad'); ylabel('CDF');");
    engEvalString(ep, "title(['CDF ZSA for ', frequencyBand,'-band']);");
    engEvalString(ep, "legend show;");



    // Освобождение памяти
    mxDestroyArray(SFArrayDU);
    mxDestroyArray(SFArrayU);
    mxDestroyArray(SFArrayS);
    mxDestroyArray(SFArrayR);
    mxDestroyArray(KArrayDU);
    mxDestroyArray(KArrayU);
    mxDestroyArray(KArrayS);
    mxDestroyArray(KArrayR);
    mxDestroyArray(DSArrayDU);
    mxDestroyArray(DSArrayU);
    mxDestroyArray(DSArrayS);
    mxDestroyArray(DSArrayR);
    mxDestroyArray(ASDArrayDU);
    mxDestroyArray(ASDArrayU);
    mxDestroyArray(ASDArrayS);
    mxDestroyArray(ASDArrayR);
    mxDestroyArray(ASAArrayDU);
    mxDestroyArray(ASAArrayU);
    mxDestroyArray(ASAArrayS);
    mxDestroyArray(ASAArrayR);
    mxDestroyArray(ZSDArrayDU);
    mxDestroyArray(ZSDArrayU);
    mxDestroyArray(ZSDArrayS);
    mxDestroyArray(ZSDArrayR);
    mxDestroyArray(ZSAArrayDU);
    mxDestroyArray(ZSAArrayU);
    mxDestroyArray(ZSAArrayS);
    mxDestroyArray(ZSAArrayR);
    mxDestroyArray(scenarioArrayDU);
    mxDestroyArray(scenarioArrayU);
    mxDestroyArray(scenarioArrayS);
    mxDestroyArray(scenarioArrayR);
    mxDestroyArray(BandsArray);
}
