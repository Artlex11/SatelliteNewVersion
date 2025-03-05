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
    void plotCDF(const std::vector<double>& DL_CNR_dB_vec, const std::vector<double>& DL_CIR_dB_vec, const std::vector<double>& DL_CINR_dB_vec, const std::string& scenario);

private:
    Engine* ep;
};


class MatlabLSPPlot
{
public:
    MatlabLSPPlot();
    ~MatlabLSPPlot();

    void plotForAllLSP(std::vector<double> SFDU, std::vector<double> SFU, std::vector<double> SFS, std::vector<double> SFR, std::vector<double> KDU, std::vector<double> KU, std::vector<double> KS, std::vector<double> KR, std::vector<double> DSDU, std::vector<double> DSU, std::vector<double> DSS, std::vector<double> DSR, std::vector<double> ASDDU, std::vector<double> ASDU, std::vector<double> ASDS, std::vector<double> ASDR, std::vector<double> ASADU, std::vector<double> ASAU, std::vector<double> ASAS, std::vector<double> ASAR, std::vector<double> ZSDDU, std::vector<double> ZSDU, std::vector<double> ZSDS, std::vector<double> ZSDR, std::vector<double> ZSADU, std::vector<double> ZSAU, std::vector<double> ZSAS, std::vector<double> ZSAR, std::vector<std::string> scenarios, std::string frequencyBand);

private:
    Engine* ep;

};

#endif // MATLAB_PLOT_H#pragma once
