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





void MatlabPlot::plotEarth() {

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