#ifndef LINKBUDGET_H
#define LINKBUDGET_H

#include <vector>
#include <Eigen/Dense>

class LinkBudget {
public:
    LinkBudget(double bandwidth_MHz, double noiseFigure_dB, double temperature_K);

    void calculateMetrics(const std::vector<double>& coreArrPatt, const std::vector<double>& maxPL_lin_magn, const Eigen::MatrixXd& allArrPatt_magn);

    const std::vector<double>& getDL_CNR_dB() const { return DL_CNR_dB_vec; }
    const std::vector<double>& getDL_CIR_dB() const { return DL_CIR_dB_vec; }
    const std::vector<double>& getDL_CINR_dB() const { return DL_CINR_dB_vec; }

private:
    double BW_MHz;
    double NF_dB;
    double T;
    double k_dBWHz;

    std::vector<double> DL_CNR_dB_vec;
    std::vector<double> DL_CIR_dB_vec;
    std::vector<double> DL_CINR_dB_vec;


};

#endif // LINKBUDGET_H