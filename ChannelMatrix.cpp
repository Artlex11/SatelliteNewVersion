#include "ChannelMatrix.h"

namespace ChannelMatrixH {


    Eigen::MatrixXcd generateRayGain(LinkData& link) {
        int u = link.antennaUE.antennaArray.size() * link.antennaUE.anglesPol.size(); // numbers ant elements in UE AntArray 
        int s = link.antennaSat.antennaArray.size() * link.antennaSat.anglesPol.size(); // numbers ant elements in Sat AntArray 
        int numbersUSonePol = u * s / 4;
        int numClusts = link.powerInRays.rows();
        int numRaysInClust = link.powerInRays.cols();

        Eigen::MatrixXcd rGainMatrix(u * s, numClusts * numRaysInClust);

        double antPatternUE = AntennaPattern::AntennaPattern_Omni();
        double antPatternSat = AntennaPattern::AntennaPattern_Omni();

        Eigen::MatrixXd XPR_n_m = link.XPR_n_m;

        // Angles -> LSC
        Eigen::MatrixXd AOD_n_m = link.AOD_n_m.array() + (link.AoD_Los - link.antennaSat.bearingAntArray);
        Eigen::MatrixXd AOA_n_m = link.AOA_n_m.array() + (link.AoA_Los - link.antennaUE.bearingAntArray);
        Eigen::MatrixXd ZOD_n_m = link.ZOD_n_m.array() + (link.ZoD_Los - 90.0 - link.antennaSat.elevationAntArray);
        Eigen::MatrixXd ZOA_n_m = link.ZOA_n_m.array() + (link.ZoA_Los - 90.0 - link.antennaUE.elevationAntArray);




        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> randPhase(-M_PI, M_PI);

        // ѕредварительное вычисление дл€ каждой пары (n, m)
        std::vector<Eigen::Matrix2cd> XPR_and_InitialRandomPhases(numClusts * numRaysInClust);
        std::vector<Eigen::Vector2d> F_theta_F_phi_TXPol1(numClusts * numRaysInClust);
        std::vector<Eigen::Vector2d> F_theta_F_phi_TXPol2(numClusts * numRaysInClust);
        std::vector<Eigen::Vector2d> F_theta_F_phi_RXPol1(numClusts * numRaysInClust);
        std::vector<Eigen::Vector2d> F_theta_F_phi_RXPol2(numClusts * numRaysInClust);

        double anglesPolTX1 = link.antennaSat.anglesPol[0];
        double anglesPolTX2 = link.antennaSat.anglesPol[1];
        double anglesPolRX1 = link.antennaUE.anglesPol[0];
        double anglesPolRX2 = link.antennaUE.anglesPol[1];





        for (int n = 0; n < numClusts; ++n) {
            for (int m = 0; m < numRaysInClust; ++m) {
                int index = n * numRaysInClust + m;
                double xpr_inv = sqrt(1. / XPR_n_m(n, m));

                XPR_and_InitialRandomPhases[index] <<
                    exp(std::complex<double>(0.0, randPhase(gen))),
                    xpr_inv* exp(std::complex<double>(0.0, randPhase(gen))),
                    xpr_inv* exp(std::complex<double>(0.0, randPhase(gen))),
                    exp(std::complex<double>(0.0, randPhase(gen)));





                F_theta_F_phi_TXPol1[index] = AntennaArray::FieldPatternWithPolarization(antPatternSat, anglesPolTX1); // -45 TX
                F_theta_F_phi_TXPol2[index] = AntennaArray::FieldPatternWithPolarization(antPatternSat, anglesPolTX2); //  45 TX
                F_theta_F_phi_RXPol1[index] = AntennaArray::FieldPatternWithPolarization(antPatternUE, anglesPolRX1);  // -45 RX
                F_theta_F_phi_RXPol2[index] = AntennaArray::FieldPatternWithPolarization(antPatternUE, anglesPolRX2);  //  45 RX
            }
        }




        // ѕараллелизаци€ внешнего цикла
#pragma omp parallel for
        for (int us = 0; us < numbersUSonePol; ++us) {
            for (int n = 0; n < numClusts; ++n) {
                for (int m = 0; m < numRaysInClust; ++m) {
                    int index = n * numRaysInClust + m;

                    auto T1 = F_theta_F_phi_RXPol1[index].transpose() * XPR_and_InitialRandomPhases[index];
                    auto T2 = F_theta_F_phi_RXPol2[index].transpose() * XPR_and_InitialRandomPhases[index];

                    rGainMatrix(us, index) = (T1 * F_theta_F_phi_TXPol1[index])(0, 0);
                    rGainMatrix(us + numbersUSonePol, index) = (T2 * F_theta_F_phi_TXPol1[index])(0, 0);
                    rGainMatrix(us + 2 * numbersUSonePol, index) = (T1 * F_theta_F_phi_TXPol2[index])(0, 0);
                    rGainMatrix(us + 3 * numbersUSonePol, index) = (T2 * F_theta_F_phi_TXPol2[index])(0, 0);
                }
            }
        }

        return rGainMatrix;
    }

}


