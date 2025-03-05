#include "Budget.h"
#include <cmath>

LinkBudget::LinkBudget(double bandwidth_MHz, double noiseFigure_dB, double temperature_K)
    : BW_MHz(bandwidth_MHz), NF_dB(noiseFigure_dB), T(temperature_K), k_dBWHz(-228.6) {
}

void LinkBudget::calculateMetrics(const std::vector<double>& coreArrPatt, const std::vector<double>& maxPL_lin_magn, const Eigen::MatrixXd& allArrPatt_magn) {
    double transmitPower_dBW = 4 + 10 * log10(BW_MHz);
    double Pn_dBW = k_dBWHz + 10 * std::log10(T) + 10 * std::log10(BW_MHz * 1e6) + NF_dB;
    double Pn_lin = std::pow(10, Pn_dBW / 10);

    DL_CNR_dB_vec.clear();
    DL_CIR_dB_vec.clear();
    DL_CINR_dB_vec.clear();

    for (size_t i = 0; i < coreArrPatt.size(); ++i) {
        double attSignal_magn = maxPL_lin_magn[i] * coreArrPatt[i];
        double attInterference_magn = std::sqrt((maxPL_lin_magn[i] * allArrPatt_magn.row(i)).array().square().sum() - (attSignal_magn * attSignal_magn));

        double Prx_lin = std::pow(10, transmitPower_dBW / 10) * attSignal_magn * attSignal_magn;
        double Prx_dBW = 10 * std::log10(Prx_lin);
        double Pint_lin = std::pow(10, transmitPower_dBW / 10) * attInterference_magn * attInterference_magn;

        double DL_CNR_dB = Prx_dBW - Pn_dBW;
        double DL_CIR_dB = Prx_dBW - 10 * std::log10(Pint_lin);
        double DL_CINR_dB = Prx_dBW - 10 * std::log10(Pn_lin + Pint_lin);

        DL_CNR_dB_vec.push_back(DL_CNR_dB);
        DL_CIR_dB_vec.push_back(DL_CIR_dB);
        DL_CINR_dB_vec.push_back(DL_CINR_dB);
    }
}

