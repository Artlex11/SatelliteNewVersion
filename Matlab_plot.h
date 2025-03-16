#ifndef MATLAB_PLOT_H
#define MATLAB_PLOT_H

#include <engine.h>
#include <Eigen/Dense>
#include <vector>

namespace MatlabPlot {

    //MatlabPlot();
    //~MatlabPlot();

    void plotTransformedData(Engine* ep, const std::vector<Eigen::Vector3d>& users, const Eigen::Vector3d& satellite);
    void plotEarth(Engine* ep);
    void plotCDF(Engine* ep, const std::vector<double>& DL_CNR_dB_vec, const std::vector<double>& DL_CIR_dB_vec, const std::vector<double>& DL_CINR_dB_vec, const std::string& scenario);


};

#endif // MATLAB_PLOT_H