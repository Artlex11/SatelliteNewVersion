#include "Matlab_plot.h"
#include <iostream>
//
//MatlabPlot::MatlabPlot() {
//    if (!(ep = engOpen(""))) {
//        std::cerr << "MATLAB Engine no open" << std::endl;
//        exit(1);
//    }
//}
//
//MatlabPlot::~MatlabPlot() {
//    engClose(ep);
//}

void MatlabPlot::plotTransformedData(Engine* ep, const std::vector<Eigen::Vector3d>& users, const Eigen::Vector3d& satellite) {
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
    engEvalString(ep, " hold on; grid on;");
    // Отрисовка линии из центра Земли (0, 0, 0) в спутник
    engEvalString(ep, "line([0, satX], [0, satY], [0, satZ], 'Color', 'b', 'LineWidth', 1);");
    engEvalString(ep, "plot3(usersX, usersY, usersZ, '.', 'DisplayName', 'Координаты пользователей');");
    engEvalString(ep, "plot3(satX, satY, satZ, '.', 'MarkerSize', 10, 'DisplayName', 'Спутник');");



    // Освобождение памяти
    mxDestroyArray(usersX);
    mxDestroyArray(usersY);
    mxDestroyArray(usersZ);
    mxDestroyArray(satX);
    mxDestroyArray(satY);
    mxDestroyArray(satZ);

}

void MatlabPlot::plotEarth(Engine* ep) {

    engEvalString(ep, "figure; hold on; grid on;");
    engEvalString(ep, "[phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));");
    engEvalString(ep, "R = 6371;");
    engEvalString(ep, "x_sphere = R * cos(phi) .* cos(theta);");
    engEvalString(ep, "y_sphere = R * sin(phi) .* cos(theta);");
    engEvalString(ep, "z_sphere = R * sin(theta);");
    engEvalString(ep, "plot3(0,0,0,'r.'); ");
    engEvalString(ep, "surf(x_sphere, y_sphere, z_sphere, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Сфера Земли');");
    engEvalString(ep, "axis equal;");
    engEvalString(ep, "xlabel('x, km'); ylabel('y, km'); zlabel('z, km');");
    engEvalString(ep, "grid on");
    engEvalString(ep, "legend show;");

}

void MatlabPlot::plotCDF(Engine* ep, const std::vector<double>& DL_CNR_dB_vec, const std::vector<double>& DL_CIR_dB_vec, const std::vector<double>& DL_CINR_dB_vec, const std::string& scenario) {
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