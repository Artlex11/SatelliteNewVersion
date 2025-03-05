#include "SSP.h"

// Определение переменных
double SSP::r_tau = 0.0;
double SSP::numberClusters = 0.0;
double SSP::ksi = 3.0;

void SSP::setParameters(Eigen::VectorXd Parameters) {
    //std::cout << "Parameters.size(): " << Parameters.size() << '\n';
    r_tau = Parameters[Parameters.size() - 11];
    numberClusters = Parameters[Parameters.size() - 8];
    ksi = Parameters[Parameters.size() - 2]; //вроде как для всех сценариев 3 ...
    //std::cout << "ksi: " << ksi << '\n';
    //std::cout << "r_tau: " << r_tau << '\n';
    //std::cout << "numberClusters: " << numberClusters << '\n';
}

//Исправленный метод 
void SSP::generateClusterDelays(LinkData& link) {

    double riceanK = link.K_db;
    double delaySpread = link.DS_sec;

    std::vector<double> delays_tau;
    double delay_tau_n;

    for (int n = 0; n < numberClusters; ++n) {
        double Xn = RandomGenerators::generateUniform(0.0, 1.0);
        delay_tau_n = -1 * r_tau * log(Xn) * delaySpread;
        delays_tau.push_back(delay_tau_n);
    }

    // Нормализация и сортировка задержек
    double minDelay_tau_n = *std::min_element(delays_tau.begin(), delays_tau.end());
    for (auto& delay_tau_n : delays_tau) {
        delay_tau_n -= minDelay_tau_n;
    }
    std::sort(delays_tau.begin(), delays_tau.end());

    link.clusterDelays = delays_tau;

    // Масштабирование для LOS
    if (link.isLos) {
        double scalingFactor = (0.000017 * riceanK * riceanK * riceanK) + (0.0002 * riceanK * riceanK) - (0.0433 * riceanK) + 0.7705; //(TR38.901v17.0.0 Rel_17, formula 7.5 - 3)
        for (int n = 0; n < delays_tau.size(); ++n) {
            delays_tau[n] /= scalingFactor; //(TR38.901v17.0.0 Rel_17, formula 7.5-4)
        }
        link.clusterScaledDelays = delays_tau;
    }


    //std::cout << "ksi: " << ksi << '\n';
    //std::cout << "r_tau: " << r_tau << '\n';
    //std::cout << "numberClusters: " << numberClusters << '\n';
}



std::pair<std::vector<double>, std::vector<double>> SSP::generateClusterPowers(LinkData& link) {

    double delaySpread = link.DS_sec;
    double riceanK = link.K_db;
    std::vector<double> clusterDelays = link.clusterDelays;
    std::vector<double> clusterPowers(clusterDelays.size());

    double power = 0.0;
    double sumclusterPowers = 0.0;

    for (size_t n = 0; n < clusterDelays.size(); ++n) {
        double shadowing = RandomGenerators::generateGauss(0.0, ksi);
        power = exp(((-1.0) * clusterDelays[n] * (r_tau - 1.0)) / r_tau / delaySpread) * pow(10, (-0.1 * shadowing));
        clusterPowers[n] = power;
    }


    // Нормализация мощностей кластеров
    for (auto& n : clusterPowers) sumclusterPowers += n;

    for (size_t n = 0; n < clusterPowers.size(); ++n) {
        clusterPowers[n] = clusterPowers[n] / sumclusterPowers;
    }

    // Добавление дополнительных коэффициентов для LOS
    std::vector<double> clusterPowersWithScalingFactors = clusterPowers;
    if (link.isLos) {
        for (size_t n = 0; n < clusterPowersWithScalingFactors.size(); ++n) {
            clusterPowersWithScalingFactors[n] *= (1 / (pow(10, riceanK / 10) + 1));
        }
        clusterPowersWithScalingFactors[0] += pow(10, riceanK / 10) / (pow(10, riceanK / 10) + 1);
    }

    // Нормализация по максимальной мощности
    double threshold = *std::max_element(clusterPowers.begin(), clusterPowers.end()) * 0.00316228; // -25 дБ
    for (size_t n = 0; n < clusterPowers.size(); ++n) {
        if (clusterPowers[n] < threshold) {
            clusterPowers[n] = 0.0;
        }
    }

    //std::vector<double> filteredDelays; // Вектор для хранения отфильтрованных задержек
    //double maxPowerWithScalingFactors = *std::max_element(clusterPowersWithScalingFactors.begin(), clusterPowersWithScalingFactors.end());
    double thresholdWithScalingFactors = *std::max_element(clusterPowersWithScalingFactors.begin(), clusterPowersWithScalingFactors.end()) * 0.00316228; // -25 дБ
    for (size_t n = 0; n < clusterPowersWithScalingFactors.size(); ++n) {
        if (clusterPowersWithScalingFactors[n] < thresholdWithScalingFactors) {
            clusterPowersWithScalingFactors[n] = 0.0;

            //filteredDelays.emplace_back(clusterDelays[n]); // Сохраняем соответствующую задержку
        }
    }
    //clusterDelays = filteredDelays; // Обновляем массив задержек

    return { clusterPowers, clusterPowersWithScalingFactors };
}

//Исправленный метод 
void SSP::calculateLosAngles(LinkData& link, Eigen::Vector3d y_local, Eigen::Vector3d z_local) {

    // Вектор "спутник-пользователь"
    Eigen::Vector3d v = link.userPosition - link.satellitePosition;
    double v_norm = v.norm();

    Eigen::Vector3d x_local = link.satellitePosition.normalized();

    // Проекции вектора на локальные оси
    double v_x_local = v.dot(x_local);
    double v_y_local = v.dot(y_local);
    double v_z_local = v.dot(z_local);

    link.AoD_Los = std::atan2(v_y_local, v_x_local) * 180.0 / M_PI;
    link.AoA_Los = std::atan2(-v_y_local, -v_x_local) * 180.0 / M_PI;
    link.ZoD_Los = std::acos(v_z_local / v_norm) * 180.0 / M_PI;
    link.ZoA_Los = std::acos(-v_z_local / v_norm) * 180.0 / M_PI;
}



//
//// Function to generate AOA or AOD based on given parameters
//Eigen::MatrixXd SSP::generateAOAandAOD_n_m(bool los, const std::vector<double>& clusterPowers,
//    double ASAorASD, double riceanK, double AOAorAOD, int AOA_0_or_AOD_1) {
//
//
//    // 
//    AOAorAOD = AOAorAOD * 180 / M_PI; // Convert to degrees
//    ASAorASD = std::pow(10, ASAorASD); // Convert to linear scale
//
//
//
//    double maxPower = *std::max_element(clusterPowers.begin(), clusterPowers.end());
//    double C_phi = 1.273 * ((0.0001 * riceanK * riceanK * riceanK) - (0.002 * riceanK * riceanK) - (0.028 * riceanK) + 1.1035);
//
//    Eigen::MatrixXd AOAorAOD_n_m(clusterPowers.size() + (los ? 1 : 0), 20);
//    Eigen::VectorXd AOAorAOD_n(clusterPowers.size());
//
//    for (size_t n = 0; n < clusterPowers.size(); ++n) {
//        double Xn = RandomGenerators::generateUniformFromVector(std::vector<int>{ -1, 1 });
//        double Yn = RandomGenerators::generateGauss(0.0, (ASAorASD / 7.0));
//
//        // Calculate AOA or AOD
//        AOAorAOD_n(n) = 2.0 * (ASAorASD / 1.4) * std::sqrt(-std::log(clusterPowers[n] / maxPower)) / C_phi;
//
//        if (los && n == 0) {
//            AOAorAOD_n_m(n, 0) = AOAorAOD; // Initialize first value
//        }
//
//        // Update the value based on previous results
//        AOAorAOD_n(n) = AOAorAOD_n(n) * Xn + Yn + AOAorAOD;
//
//        for (int m = 0; m < 20; ++m) {
//            AOAorAOD_n_m(n + (los ? 1 : 0), m) = AOAorAOD_n(n) + (los ? los_C[AOA_0_or_AOD_1] : nlos_C[AOA_0_or_AOD_1]) * am[m];
//            // Normalize angles to [-180, 180]
//            while (AOAorAOD_n_m(n + (los ? 1 : 0), m) < -180) AOAorAOD_n_m(n + (los ? 1 : 0), m) += 360;
//            while (AOAorAOD_n_m(n + (los ? 1 : 0), m) > 180) AOAorAOD_n_m(n + (los ? 1 : 0), m) -= 360;
//        }
//    }
//    return AOAorAOD_n_m;
//}
//
//// Function to generate ZOA or ZOD based on given parameters
//Eigen::MatrixXd SSP::generateZOAandZOD_n_m(bool los, const std::vector<double>& clusterPowers,
//    double ZSAorZSD, double riceanK, double ZOAorZOD, int ZOA_2_or_ZOD_3) {
//    ZOAorZOD = ZOAorZOD * 180 / M_PI; // Convert to degrees
//    ZSAorZSD = std::pow(10, ZSAorZSD); // Convert to linear scale
//
//    double maxPower = *std::max_element(clusterPowers.begin(), clusterPowers.end());
//    double C_theta = 1.184 * ((0.0002 * riceanK * riceanK * riceanK) - (0.0077 * riceanK * riceanK) + (0.0339 * riceanK) + 1.3086);
//
//    Eigen::MatrixXd ZOAorZOD_n_m(clusterPowers.size() + (los ? 1 : 0), 20);
//    Eigen::VectorXd ZOAorZOD_n(clusterPowers.size());
//
//    for (size_t n = 0; n < clusterPowers.size(); ++n) {
//        double Xn = RandomGenerators::generateUniformFromVector(std::vector<int>{ -1, 1 });
//        double Yn = RandomGenerators::generateGauss(0.0, (ZSAorZSD / 7.0));
//
//        // Calculate ZOA or ZOD
//        ZOAorZOD_n(n) = -1 * ZSAorZSD * std::log(clusterPowers[n] / maxPower) / C_theta;
//
//        if (los && n == 0) {
//            ZOAorZOD_n_m(n, 0) = ZOAorZOD; // Initialize first value
//        }
//
//        // Update the value based on previous results
//        ZOAorZOD_n(n) = ZOAorZOD_n(n) * Xn + Yn + ZOAorZOD;
//
//        for (int m = 0; m < 20; ++m) {
//            ZOAorZOD_n_m(n + (los ? 1 : 0), m) = ZOAorZOD_n(n) + (los ? los_C[ZOA_2_or_ZOD_3] : nlos_C[ZOA_2_or_ZOD_3]) * am[m];
//            // Normalize angles to [0, 360]
//            if (ZOAorZOD_n_m(n + (los ? 1 : 0), m) >= 180) {
//                ZOAorZOD_n_m(n + (los ? 1 : 0), m) = 360 - ZOAorZOD_n_m(n + (los ? 1 : 0), m);
//            }
//        }
//    }
//    return ZOAorZOD_n_m;
//}
