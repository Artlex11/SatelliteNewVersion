#include "AntennasPatterns.h"

namespace AntennaPattern {

    double AntennaPattern_NTN1_3D(double az_deg, double el_deg)
    {
        double AP_gain = 4; // dBi
        double AP_max = 25; // dB
        double phi_3dB = 65; //degrees

        auto AntennaPattern_NTN1_V = [&](double el_deg) -> double {
            el_deg = el_deg > 180 ? el_deg - 360 : el_deg;
            el_deg = el_deg < -180 ? el_deg + 360 : el_deg;
            return -std::min(12 * pow((el_deg - 90) / phi_3dB, 2), AP_max);
            };

        auto AntennaPattern_NTN1_H = [&](double az_deg) -> double {
            az_deg = az_deg > 180 ? az_deg - 360 : az_deg;
            az_deg = az_deg < -180 ? az_deg + 360 : az_deg;
            return -std::min(12 * pow(az_deg / phi_3dB, 2), AP_max);
            };


        return AP_gain - std::min(-(AntennaPattern_NTN1_V(el_deg) + AntennaPattern_NTN1_H(az_deg)), AP_max);
    }

    double AntennaPattern_Omni()
    {
        return 0.0; // dB
    }

    double AntennaPattern_38901_3D(double az_deg, double el_deg)
    {
        double AP_gain = 8; // dBi
        double AP_max = 30; // dB
        double phi_3dB = 65; //degrees

        auto AntennaPattern_38901_V = [&](double el_deg) -> double {
            el_deg = el_deg > 180 ? el_deg - 360 : el_deg;
            el_deg = el_deg < -180 ? el_deg + 360 : el_deg;
            return -std::min(12 * ((el_deg - 90) / phi_3dB) * ((el_deg - 90) / phi_3dB), AP_max);
            };

        auto AntennaPattern_38901_H = [&](double az_deg) -> double {
            az_deg = az_deg > 180 ? az_deg - 360 : az_deg;
            az_deg = az_deg < -180 ? az_deg + 360 : az_deg;
            return -std::min(12 * (az_deg / phi_3dB) * (az_deg / phi_3dB), AP_max);
            };


        return AP_gain - std::min(-(AntennaPattern_38901_V(el_deg) + AntennaPattern_38901_H(az_deg)), AP_max);
    }

    double AntennaPattern_37885_PedCellUE_3D(double az_deg, double el_deg)
    {
        double AP_gain = 5; // dBi
        double AP_max = 25; // dB
        double phi_3dB = 90; //degrees

        auto AntennaPattern_37885_PedCellUE_V = [&](double el_deg) -> double {
            el_deg = el_deg > 180 ? el_deg - 360 : el_deg;
            el_deg = el_deg < -180 ? el_deg + 360 : el_deg;
            return -std::min(12 * ((el_deg - 90) / phi_3dB) * ((el_deg - 90) / phi_3dB), AP_max);
            };

        auto AntennaPattern_37885_PedCellUE_H = [&](double az_deg) -> double {
            az_deg = az_deg > 180 ? az_deg - 360 : az_deg;
            az_deg = az_deg < -180 ? az_deg + 360 : az_deg;
            return -std::min(12 * (az_deg / phi_3dB) * (az_deg / phi_3dB), AP_max);
            };


        return AP_gain - std::min(-(AntennaPattern_37885_PedCellUE_V(el_deg) + AntennaPattern_37885_PedCellUE_H(az_deg)), AP_max);
    }



}

namespace AntennaArray {
    Eigen::MatrixXd createAntennaArray_Mg_Ng_P_M_N(int Mg, int Ng, int M, int N, double frequency, double spacing)
    {
        double lambda = 3e8 / (frequency * 1e9);
        double elementSpacing = lambda * spacing;

        Eigen::MatrixXd antennaMatrix(Mg * M * Ng * N, 3);
        int index = 0;


        double maxY = (M * Mg - 1) * elementSpacing;
        double maxZ = (N * Ng - 1) * elementSpacing;


        double centerY = maxY / 2.0;
        double centerZ = maxZ / 2.0;


        for (int i = 0; i < Mg; ++i) {
            for (int m = 0; m < M; ++m) {
                for (int j = 0; j < Ng; ++j) {
                    for (int n = 0; n < N; ++n) {
                        double elementY = (i * M + m) * elementSpacing - centerY;
                        double elementZ = (j * N + n) * elementSpacing - centerZ;

                        antennaMatrix(index, 0) = 0.0;
                        antennaMatrix(index, 1) = elementY;
                        antennaMatrix(index, 2) = elementZ;
                        index++;
                    }
                }
            }
        }
        std::cout << antennaMatrix.transpose() << "\n";
        return antennaMatrix.transpose();

    }

    Eigen::MatrixXd rotateAntennaArray(double alpha, double beta, const Eigen::MatrixXd& antennaArray) {

        double alpha_rad = alpha * M_PI / 180.0;
        double beta_rad = beta * M_PI / 180.0;

        // Исходная нормаль 
        Eigen::Vector3d panNormal = { 1.0, 0.0, 0.0 };

        // Шаг 1: Поворот вокруг оси Z на азимутальный угол (alpha)
        Eigen::Vector3d axisZ(0.0, 0.0, 1.0);
        Eigen::Matrix3d rotationZ = Eigen::AngleAxisd(alpha_rad, axisZ).toRotationMatrix();
        panNormal = rotationZ * panNormal;

        // Шаг 2: Находим новую ось для поворота на угол места (beta)
        Eigen::Vector3d ax = panNormal.cross(axisZ).normalized();

        // Шаг 3: Поворот вокруг новой оси на угол места (beta)
        Eigen::Matrix3d rotationX = Eigen::AngleAxisd(beta_rad, ax).toRotationMatrix();

        // Комбинируем повороты
        Eigen::Matrix3d totalRotation = rotationX * rotationZ;


        Eigen::MatrixXd rotatedAntennaArray = totalRotation * antennaArray;

        // В общем случае необходим поворот по ещё одной оси , но пока не нужен  такой функционал

        return rotatedAntennaArray;
    }

    void addAntennaUE(LinkData& link, int numElemensInCol, int numElemensInRow, std::string pattern, double frequency, std::vector<double> anglesPol, double spacing)
    {
        Eigen::MatrixXd antennaArray = createAntennaArray_Mg_Ng_P_M_N(1, 1, numElemensInRow, numElemensInCol, frequency, spacing);

        double bearing = link.AoA_Los + 180.0;
        double elevation = link.ZoA_Los - 90.0;
        bearing = bearing > 180 ? bearing - 360 : bearing;
        antennaArray = rotateAntennaArray(bearing, elevation, antennaArray);

        link.antennaUE = { antennaArray, "Omni", anglesPol ,bearing, elevation };
    }

    void addAntennaSat(LinkData& link, int numElemensInCol, int numElemensInRow, std::string pattern, double frequency, std::vector<double> anglesPol, double spacing)
    {
        Eigen::MatrixXd antennaArray = createAntennaArray_Mg_Ng_P_M_N(1, 1, numElemensInRow, numElemensInCol, frequency, spacing);

        Eigen::Vector3d xyzSat = link.satellitePosition;
        double bearing = (atan2(xyzSat(1), xyzSat(0)) * 180.0 / M_PI) + 180.0;
        double elevation = (atan2(xyzSat(2), sqrt(xyzSat(0) * xyzSat(0) + xyzSat(1) * xyzSat(1))) * 180. / M_PI) - 90.0;
        antennaArray = rotateAntennaArray(bearing, elevation, antennaArray);

        link.antennaSat = { antennaArray, "Omni", anglesPol , bearing , elevation };
    }

    Eigen::Vector2d FieldPatternWithPolarization(double& antennaPattern3D_dB, double anglePol_deg)
    {
        double fieldPattern3D_lin = pow(10., antennaPattern3D_dB / 20.);
        double anglePol_rad = anglePol_deg * M_PI / 180.;
        return { fieldPattern3D_lin * std::cos(anglePol_rad), fieldPattern3D_lin * std::sin(anglePol_rad) }; // F_theta (vertical) , F_phi (horizontal)
    }
}