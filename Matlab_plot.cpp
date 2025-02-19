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
    engEvalString(ep, "figure;");
    engEvalString(ep, "plot3(usersX, usersY, usersZ, '.', 'DisplayName', '��������������� ���������� �������������'); hold on;");
    engEvalString(ep, "plot3(satX, satY, satZ, 'b*', 'MarkerSize', 10, 'DisplayName', '��������������� �������');");

    // ������������ ������
    mxDestroyArray(usersX);
    mxDestroyArray(usersY);
    mxDestroyArray(usersZ);
    mxDestroyArray(satX);
    mxDestroyArray(satY);
    mxDestroyArray(satZ);
}

void MatlabPlot::plotEarth() {
    // ���������� ����� ����� � MATLAB
    engEvalString(ep, "[phi, theta] = meshgrid(linspace(0, 2*pi, 30), linspace(-pi/2, pi/2, 30));");
    engEvalString(ep, "R = 6371;");
    engEvalString(ep, "x_sphere = R * cos(phi) .* cos(theta);");
    engEvalString(ep, "y_sphere = R * sin(phi) .* cos(theta);");
    engEvalString(ep, "z_sphere = R * sin(theta);");
    engEvalString(ep, "surf(x_sphere, y_sphere, z_sphere, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', '����� �����');");
    engEvalString(ep, "axis equal;");
    engEvalString(ep, "xlabel('x, km'); ylabel('y, km'); zlabel('z, km');");
    engEvalString(ep, "grid on");
    engEvalString(ep, "legend show; hold off;");
}