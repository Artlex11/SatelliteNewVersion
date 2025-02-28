#ifndef MATLAB_PLOT_H
#define MATLAB_PLOT_H

#include <engine.h>
#include <Eigen/Dense>
#include <vector>

class MatlabPlot {
public:
    MatlabPlot();
    ~MatlabPlot();

    void plotTransformedData(const std::vector<Eigen::Vector3d>& users, const Eigen::Vector3d& satellite);
    void plotEarth();
    void plotRayPoints(const Eigen::Vector3d& satellitePosition, const std::vector<Eigen::Vector3d>& rays);

private:
    Engine* ep;
};

#endif // MATLAB_PLOT_H#pragma once
