#include "Matlab_plot.h"
#include <iostream>

MatlabPlot::MatlabPlot() {
    if (!(ep = engOpen(""))) {
        std::cerr << "�� ������� ������� MATLAB Engine" << std::endl;
        exit(1);
    }
}

MatlabPlot::~MatlabPlot() {
    engClose(ep);
}

void MatlabPlot::plotTransformedData(const std::vector<Eigen::Vector3d>& users, const Eigen::Vector3d& satellite) {
    // �������� ������ ������������� � MATLAB
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

    // �������� ������ �������� � MATLAB
    mxArray* satX = mxCreateDoubleScalar(satellite.x());
    mxArray* satY = mxCreateDoubleScalar(satellite.y());
    mxArray* satZ = mxCreateDoubleScalar(satellite.z());

    engPutVariable(ep, "satX", satX);
    engPutVariable(ep, "satY", satY);
    engPutVariable(ep, "satZ", satZ);

    // ���������� �������� � MATLAB
    engEvalString(ep, "hold on; grid on;"); // ��������� ����� ������ �� ������� ������
    engEvalString(ep, "plot3(usersX, usersY, usersZ, '.', 'DisplayName', '���������� �������������');");
    engEvalString(ep, "plot3(satX, satY, satZ, '.', 'MarkerSize', 10, 'DisplayName', '�������');");

    engEvalString(ep, "plot3(0,0,0,'r.'); axis equal;xlabel('x, km'); ylabel('y, km'); zlabel('z, km'); ");

    // ������������ ������
    mxDestroyArray(usersX);
    mxDestroyArray(usersY);
    mxDestroyArray(usersZ);
    mxDestroyArray(satX);
    mxDestroyArray(satY);
    mxDestroyArray(satZ);

}

void MatlabPlot::plotRayPoints(const Eigen::Vector3d& satellitePosition, const std::vector<Eigen::Vector3d>& rays) {
    // �������� ������
    std::cout << "Satellite position: " << satellitePosition.transpose() << std::endl;
    for (size_t i = 0; i < rays.size(); ++i) {
        std::cout << "Ray " << i + 1 << ": " << rays[i].transpose() << std::endl;
    }

    // �������� �� ������� �����
    if (rays.empty()) {
        std::cerr << "������: ������ ����� ����." << std::endl;
        return;
    }

    // ���������� �������� ����� �����
    std::vector<Eigen::Vector3d> rayPoints;
    for (const auto& ray : rays) {
        rayPoints.push_back(satellitePosition + ray);
    }

    // �������� ������ � MATLAB
    mxArray* satPosArray = mxCreateDoubleMatrix(1, 3, mxREAL);
    std::memcpy(mxGetPr(satPosArray), satellitePosition.data(), 3 * sizeof(double));
    engPutVariable(ep, "satPos", satPosArray);

    // �������� ������� ��� �������� ����� �����
    mxArray* rayPointsArray = mxCreateDoubleMatrix(rayPoints.size(), 3, mxREAL);
    double* rayPointsData = mxGetPr(rayPointsArray);

    // ���������� ������� �������
    for (size_t i = 0; i < rayPoints.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            rayPointsData[i + j * rayPoints.size()] = rayPoints[i](j);  // ���������������� ������
        }
    }
    engPutVariable(ep, "rayPoints", rayPointsArray);

    // ��������� ��������
    engEvalString(ep, " hold on;");

    // ��������� �������� ����� �����
    engEvalString(ep, "plot3(satPos(1), satPos(2), satPos(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '�������');");

    engEvalString(ep, "plot3(rayPoints(:, 1), rayPoints(:, 2), rayPoints(:, 3),'.','DisplayName', '����� �����');");

    // ��������� �������
    engEvalString(ep, "grid on;");
    engEvalString(ep, "xlabel('X'); ylabel('Y'); zlabel('Z');");
    engEvalString(ep, "title('����� ����� �� ��������'); legend show;");
    engEvalString(ep, "hold off;");

    // ������������ ��������
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
    engEvalString(ep, "surf(x_sphere, y_sphere, z_sphere, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', '����� �����');");
    engEvalString(ep, "axis equal;");
    engEvalString(ep, "xlabel('x, km'); ylabel('y, km'); zlabel('z, km');");
    engEvalString(ep, "grid on");
    engEvalString(ep, "legend show;");


}