#include "SmallScaleParameters.h"

// Определение переменных
double r_tau = 0.0;
int numberClusters = 0.0;
int numberRaysPerCluster = 20;
double ksi = 3.0;

double clustASD = 0.0;
double clustASA = 0.0;
double clustZSD = 0.0;
double clustZSA = 0.0;

double C_phi = 0.0;
double C_theta = 0.0;

double maxPower;

double meanXPR = 0.0;
double standardDeviationXPR = 0.0;

double clucterDS_ns = 0.0;


void SSP::setParameters(Eigen::VectorXd Parameters) 
{
    //std::cout << "Parameters.size(): " << Parameters.size() << '\n';
    r_tau = Parameters[Parameters.size() - 11];

    numberClusters = Parameters[Parameters.size() - 8];
    //ksi = Parameters[Parameters.size() - 2]; //вроде как для всех сценариев 3 ...

    clustASD = Parameters[Parameters.size() - 5];
    clustASA = Parameters[Parameters.size() - 4];
    clustZSA = Parameters[Parameters.size() - 3];
    clustZSD = 0.375 * pow(10, Parameters[8]); // Как по формуле  3GPP_38901 (7.5-20) 

    C_phi = tableScalingFactorsAOAandAOD[numberClusters - 1];
    C_theta = tableScalingFactorsZOAandZOD[numberClusters - 1];

    meanXPR = Parameters[Parameters.size() - 10];
    standardDeviationXPR = Parameters[Parameters.size() - 9];

    clucterDS_ns = Parameters[Parameters.size() - 6];

    //std::cout << "ksi: " << ksi << '\n';
    //std::cout << "r_tau: " << r_tau << '\n';
    //std::cout << "numberClusters: " << numberClusters << '\n';
    //std::cout << "clustASD: " << clustASD << '\n';
    //std::cout << "clustASA: " << clustASA << '\n';
    //std::cout << "clustZSD: " << clustZSD << '\n';
    //std::cout << "clustZSA: " << clustZSA << '\n';
    //std::cout << "C_Nlos_phi: " << C_Nlos_phi << '\n';
    //std::cout << "C_Nlos_theta: " << C_Nlos_theta << '\n';
}

//Исправленный метод 
void SSP::generateClusterDelays(LinkData& link) 
{
    double riceanK = link.K_db;
    double delaySpread = link.DS_sec;

    std::vector<double> delays_tau;
    double delay_tau_n;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int n = 0; n < numberClusters; ++n) 
    {
        double Xn = dis(gen);// RandomGenerators::generateUniform(0.0, 1.0);
        delay_tau_n = (-1.0) * r_tau * log(Xn) * delaySpread;
        delays_tau.push_back(delay_tau_n);
    }

    // Нормализация и сортировка задержек
    double minDelay_tau_n = *std::min_element(delays_tau.begin(), delays_tau.end());
    for (int i = 0; i < numberClusters; ++i) 
    {
        delays_tau[i] -= minDelay_tau_n;
    }
    std::sort(delays_tau.begin(), delays_tau.end());

    link.clusterDelays = delays_tau;

    // Масштабирование для LOS
    if (link.isLos)
    {
        double scalingFactor = (0.000017 * riceanK * riceanK * riceanK) + (0.0002 * riceanK * riceanK) - (0.0433 * riceanK) + 0.7705; //(TR38.901v17.0.0 Rel_17, formula 7.5 - 3)
        for (int n = 0; n < numberClusters; ++n) 
        {
            delays_tau[n] /= scalingFactor; //(TR38.901v17.0.0 Rel_17, formula 7.5-4)
        }
    }
    link.clusterScaledDelays = delays_tau;
}



std::pair<std::vector<double>, std::vector<double>> SSP::generateClusterPowers(LinkData& link) 
{
    double delaySpread = link.DS_sec;

    std::vector<double> clusterDelays = link.clusterDelays;
    std::vector<double> clusterPowers;

    double power;
    double sumClusterPowers = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> gaussDist(0.0, 1.0);

    for (int n = 0; n < numberClusters; ++n) 
    {
        double shadowing = (-1) * ksi * gaussDist(gen);

        power = exp(clusterDelays[n] * (1.0 - r_tau) / (r_tau * delaySpread)) * pow(10.0, (shadowing / 10.0));
        clusterPowers.push_back(power);
        sumClusterPowers += power;
    }

    for (int n = 0; n < numberClusters; ++n)
    {
        clusterPowers[n] = clusterPowers[n] / sumClusterPowers;
    }

    std::vector<double> clusterPowersWithScalingFactors = clusterPowers;

    maxPower = 0.0;

    if (link.isLos) 
    {
        double K_lin = pow(10.0, link.K_db / 10.0);

        for (int n = 0; n < numberClusters; ++n) 
        {
            clusterPowersWithScalingFactors[n] = clusterPowersWithScalingFactors[n] / (K_lin + 1.0);
        }
        clusterPowersWithScalingFactors[0] = clusterPowersWithScalingFactors[0] + (K_lin / (K_lin + 1.0));

        maxPower = *std::max_element(clusterPowersWithScalingFactors.begin(), clusterPowersWithScalingFactors.end());

        for (int n = 0; n < numberClusters; ++n) 
        {
            if (log10(clusterPowersWithScalingFactors[n] / maxPower) < (-2.5)) 
            {
                clusterPowers[n] = 0.0;
                clusterPowersWithScalingFactors[n] = 0.0;
            }
        }
    }

    if (!link.isLos)
    {
        maxPower = *std::max_element(clusterPowers.begin(), clusterPowers.end());

        for (int n = 0; n < numberClusters; ++n) 
        {
            if (log10(clusterPowers[n] / maxPower) < (-2.5))
            {
                clusterPowers[n] = 0.0;
                clusterPowersWithScalingFactors[n] = 0.0;
            }
        }
    }



    return { clusterPowers, clusterPowersWithScalingFactors };
}



//Исправленный метод 
void SSP::calculateLosAngles(LinkData& link, Eigen::Vector3d y_local, Eigen::Vector3d z_local)
{
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
//void SSP::generateArrivalAndDepartureAngles(LinkData& link, std::vector<double> clusterPowersWithScalingFactors) {
//    std::random_device rd;
//    std::mt19937 gen(rd());
//
//    std::uniform_int_distribution<size_t> discreteDist(0, 1);
//    std::vector<int> discreteValues = {-1, 1};
//    std::normal_distribution<> gaussDist(0.0, 1.0);
//
//    double kFactor = link.K_db;
//    double ASD = link.ASD_deg;
//    double ASA = link.ASA_deg;
//    double ZSD = link.ZSD_deg;
//    double ZSA = link.ZSA_deg;
//
//    numberClusters = clusterPowersWithScalingFactors.size();
//
//    if (link.isLos) {
//        C_phi *= ((0.0001 * kFactor * kFactor * kFactor) - (0.002 * kFactor * kFactor) - (0.028 * kFactor) + 1.1035);
//        C_theta *= ((0.0002 * kFactor * kFactor * kFactor) - (0.0077 * kFactor * kFactor) + (0.0339 * kFactor) + 1.3086);
//    }
//
//    int numRows = link.isLos ? numberClusters + 1 : numberClusters;
//    Eigen::MatrixXd AOD_n_m(numRows, numberRaysPerCluster);
//    Eigen::MatrixXd AOA_n_m(numRows, numberRaysPerCluster);
//    Eigen::MatrixXd ZOD_n_m(numRows, numberRaysPerCluster);
//    Eigen::MatrixXd ZOA_n_m(numRows, numberRaysPerCluster);
//
//    Eigen::VectorXd AOD_n(numberClusters);
//    Eigen::VectorXd AOA_n(numberClusters);
//    Eigen::VectorXd ZOD_n(numberClusters);
//    Eigen::VectorXd ZOA_n(numberClusters);
//
//    double invC_phi = 1.0 / C_phi;
//    double invC_theta = 1.0 / C_theta;
//    double ASD_1_4 = ASD / 1.4;
//    double ASA_1_4 = ASA / 1.4;
//    double ZSD_7 = ZSD / 7.0;
//    double ZSA_7 = ZSA / 7.0;
//
//    for (int n = 0; n < numberClusters; ++n) {
//        double power_n = clusterPowersWithScalingFactors[n];
//        double logPower = std::log(power_n / maxPower);
//
//        AOD_n(n) = 2.0 * ASD_1_4 * std::sqrt(-logPower) * invC_phi;
//        AOA_n(n) = 2.0 * ASA_1_4 * std::sqrt(-logPower) * invC_phi;
//        ZOD_n(n) = -ZSD * logPower * invC_theta;
//        ZOA_n(n) = -ZSA * logPower * invC_theta;
//
//        AOD_n(n) = AOD_n(n) * discreteValues[discreteDist(gen)] + gaussDist(gen) * ZSD_7 + (link.isLos ? link.AoD_Los : 0.0);
//        AOA_n(n) = AOA_n(n) * discreteValues[discreteDist(gen)] + gaussDist(gen) * ZSA_7 + (link.isLos ? link.AoA_Los : 0.0);
//        ZOD_n(n) = ZOD_n(n) * discreteValues[discreteDist(gen)] + gaussDist(gen) * ZSD_7 + (link.isLos ? link.ZoD_Los : 0.0);
//        ZOA_n(n) = ZOA_n(n) * discreteValues[discreteDist(gen)] + gaussDist(gen) * ZSA_7 + (link.isLos ? link.ZoA_Los : 0.0);
//    }
//
//    if (link.isLos) {
//        AOD_n_m(0, 0) = link.AoD_Los;
//        AOA_n_m(0, 0) = link.AoA_Los;
//        ZOD_n_m(0, 0) = link.ZoD_Los;
//        ZOA_n_m(0, 0) = link.ZoA_Los;
//    }
//
//    for (int n = 0; n < numberClusters; ++n) {
//        for (int m = 0; m < numberRaysPerCluster; ++m) {
//            double rayOffsetAngle = rayOfsetAngles[m];
//            AOD_n_m(n + link.isLos, m) = AOD_n(n) + rayOffsetAngle * clustASD;
//            AOA_n_m(n + link.isLos, m) = AOA_n(n) + rayOffsetAngle * clustASA;
//            ZOD_n_m(n + link.isLos, m) = ZOD_n(n) + rayOffsetAngle * clustZSD;
//            ZOA_n_m(n + link.isLos, m) = ZOA_n(n) + rayOffsetAngle * clustZSA;
//
//            if (ZOD_n_m(n + link.isLos, m) >= 180.0) ZOD_n_m(n + link.isLos, m) = 360.0 - ZOD_n_m(n + link.isLos, m);
//            if (ZOA_n_m(n + link.isLos, m) >= 180.0) ZOA_n_m(n + link.isLos, m) = 360.0 - ZOA_n_m(n + link.isLos, m);
//        }
//    }
//
//    link.AOD_n_m = AOD_n_m;
//    link.AOA_n_m = AOA_n_m;
//    link.ZOD_n_m = ZOD_n_m;
//    link.ZOA_n_m = ZOA_n_m;
//}



void SSP::generateArrivalAndDepartureAngles(LinkData& link, std::vector<double> clusterPowersWithScalingFactors) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> discreteDist(0, 1);
    std::vector<int> discreteValues = { -1, 1 };
    std::normal_distribution<> gaussDist(0.0, 1.0);

    double kFactor = link.K_db;
    double ASD = link.ASD_deg;
    double ASA = link.ASA_deg;
    double ZSD = link.ZSD_deg;
    double ZSA = link.ZSA_deg;

    int numRows = numberClusters + (link.isLos ? 1 : 0);
    Eigen::MatrixXd AOD_n_m(numRows, 20);
    Eigen::MatrixXd AOA_n_m(numRows, 20);
    Eigen::MatrixXd ZOD_n_m(numRows, 20);
    Eigen::MatrixXd ZOA_n_m(numRows, 20);

    Eigen::VectorXd AOD_n(numberClusters);
    Eigen::VectorXd AOA_n(numberClusters);
    Eigen::VectorXd ZOD_n(numberClusters);
    Eigen::VectorXd ZOA_n(numberClusters);

    double ASD_1_4 = ASD / 1.4;
    double ASA_1_4 = ASA / 1.4; 
    double ASD_7 = ASD / 7.0;
    double ASA_7 = ASA / 7.0; 
    double ZSD_7 = ZSD / 7.0; 
    double ZSA_7 = ZSA / 7.0;

    if (link.isLos) 
    {
        C_phi = C_phi * ((0.0001 * kFactor * kFactor * kFactor) - (0.002 * kFactor * kFactor) - (0.028 * kFactor) + 1.1035);
        C_theta = C_theta * ((0.0002 * kFactor * kFactor * kFactor) - (0.0077 * kFactor * kFactor) + (0.0339 * kFactor) + 1.3086);

        double invC_phi = 1.0 / C_phi;
        double invC_theta = 1.0 / C_theta;

        for (int n = 0; n < numberClusters; ++n) 
        {
            double logPower = std::log(clusterPowersWithScalingFactors[n] / maxPower);

            // (7.5-9)
            AOD_n(n) = 2.0 * ASD_1_4 * std::sqrt(-logPower) * invC_phi;
            AOA_n(n) = 2.0 * ASA_1_4 * std::sqrt(-logPower) * invC_phi;
            // (7.5-14)
            ZOD_n(n) = -ZSD * logPower * invC_theta;
            ZOA_n(n) = -ZSA * logPower * invC_theta;

            // (7.5-11) и (7.5-16) without LOS term
            AOD_n(n) = (AOD_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ASD_7);
            AOA_n(n) = (AOA_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ASA_7);
            ZOD_n(n) = (ZOD_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ZSD_7);
            ZOA_n(n) = (ZOA_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ZSA_7);
        }

        // Использование первого кластера для направления LOS
        double AOD_1 = AOD_n(0);
        double AOA_1 = AOA_n(0);
        double ZOD_1 = ZOD_n(0);
        double ZOA_1 = ZOA_n(0);

        for (int n = 0; n < numberClusters; ++n) 
        {
            AOD_n(n) = AOD_n(n) - AOD_1;
            AOA_n(n) = AOA_n(n) - AOA_1;
            ZOD_n(n) = ZOD_n(n) - ZOD_1;
            ZOA_n(n) = ZOA_n(n) - ZOA_1;

            for (int m = 0; m < numberRaysPerCluster; ++m) 
            {
                auto rayOffsetAngle = rayOfsetAngles[m];

                // Нормализация углов в диапазон [-180, 180]
                auto normalizeAngle = [](double angle)
                    {
                        angle = fmod(angle, 360.0);  // Приводим к [0, 360)
                        if (angle < 0)
                        {
                            angle += 360.0;
                        }
                        return angle;
                    };

                AOD_n_m(n + 1, m) = fmod(AOD_n(n) + rayOffsetAngle * clustASD, 180);
                AOA_n_m(n + 1, m) = fmod(AOA_n(n) + rayOffsetAngle * clustASA, 180);
                ZOD_n_m(n + 1, m) = normalizeAngle(ZOD_n(n) + rayOffsetAngle * clustZSD + 90.0);
                ZOA_n_m(n + 1, m) = normalizeAngle(ZOA_n(n) + rayOffsetAngle * clustZSA + (90.0 - link.elevationAngle));
            }
        }
    }
    else
    {
        double invC_phi = 1.0 / C_phi;
        double invC_theta = 1.0 / C_theta;

        for (int n = 0; n < numberClusters; ++n) 
        {
            double logPower = std::log(clusterPowersWithScalingFactors[n] / maxPower);
            AOD_n(n) = 2.0 * ASD_1_4 * std::sqrt(-logPower) * invC_phi;
            AOA_n(n) = 2.0 * ASA_1_4 * std::sqrt(-logPower) * invC_phi;
            ZOD_n(n) = -ZSD * logPower * invC_theta;
            ZOA_n(n) = -ZSA * logPower * invC_theta;

            AOD_n(n) = (AOD_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ASD_7);
            AOA_n(n) = (AOA_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ASA_7);
            ZOD_n(n) = (ZOD_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ZSD_7);
            ZOA_n(n) = (ZOA_n(n) * discreteValues[discreteDist(gen)]) + (gaussDist(gen) * ZSA_7);

            for (int m = 0; m < numberRaysPerCluster; ++m) 
            {
                auto rayOffsetAngle = rayOfsetAngles[m];

                auto normalizeAngle = [](double angle) 
                    {
                    angle = fmod(angle, 360.0);  // Приводим к [0, 360)
                    if (angle < 0) { angle += 360.0; }
                    if (angle > 180 && angle < 360) 
                    { 
                        angle = 360 - angle; 
                    }
                    return angle;
                    };

                AOD_n_m(n, m) = fmod(AOD_n(n) + rayOffsetAngle * clustASD, 180);
                AOA_n_m(n, m) = fmod(AOA_n(n) + rayOffsetAngle * clustASA, 180);
                ZOD_n_m(n, m) = normalizeAngle(ZOD_n(n) + rayOffsetAngle * clustZSD + 90.0);
                ZOA_n_m(n, m) = normalizeAngle(ZOA_n(n) + rayOffsetAngle * clustZSA + (90.0 - link.elevationAngle));
            }
        }
    }

    // Сохранение результатов
    link.AOD_n_m = AOD_n_m;
    link.AOA_n_m = AOA_n_m;
    link.ZOD_n_m = ZOD_n_m;
    link.ZOA_n_m = ZOA_n_m;
}

void SSP::generateXPR(LinkData& link) 
{
    Eigen::MatrixXd XPR(numberClusters, numberRaysPerCluster);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(meanXPR, standardDeviationXPR);

    // Генерация гауссовских значений для всех элементов сразу
    Eigen::VectorXd gaussValues = Eigen::VectorXd::NullaryExpr(numberClusters * numberRaysPerCluster,
        [&]() { return d(gen); });

    // Преобразование из дБ в разы
    XPR = (gaussValues.array() / 10.0).pow(10);

    link.XPR_n_m = XPR;
}














//// Function to calculate angular spread and mean angles
//std::vector<double> SSP::calculateAngularSpreadandMeanAngles(bool los, const std::vector<double>& clusterPowers,
//    Eigen::MatrixXd& AOD, Eigen::MatrixXd& AOA,
//    Eigen::MatrixXd& ZOD, Eigen::MatrixXd& ZOA) {
//    std::vector<double> ASandMeanAnglesforAOD_AOA_ZOD_ZOA;
//    AOD *= M_PI / 180;
//    AOA *= M_PI / 180;
//    ZOD *= M_PI / 180;
//    ZOA *= M_PI / 180;
//
//    std::complex<double> weighted_sumAOA(0.0, 0.0);
//    std::complex<double> weighted_sumAOD(0.0, 0.0);
//    std::complex<double> weighted_sumZOA(0.0, 0.0);
//    std::complex<double> weighted_sumZOD(0.0, 0.0);
//    double weighted_sumPowers = 0.0;
//
//    // Calculate weighted sums
//    for (size_t n = 0; n < clusterPowers.size(); ++n) {
//        for (int m = 0; m < 20; ++m) {
//            double weight = clusterPowers[n] / 20;
//            weighted_sumAOD += weight * std::complex<double>(std::cos(AOD(n + (los ? 1 : 0), m)), std::sin(AOD(n + (los ? 1 : 0), m)));
//            weighted_sumAOA += weight * std::complex<double>(std::cos(AOA(n + (los ? 1 : 0), m)), std::sin(AOA(n + (los ? 1 : 0), m)));
//            weighted_sumZOD += weight * std::complex<double>(std::cos(ZOD(n + (los ? 1 : 0), m)), std::sin(ZOD(n + (los ? 1 : 0), m)));
//            weighted_sumZOA += weight * std::complex<double>(std::cos(ZOA(n + (los ? 1 : 0), m)), std::sin(ZOA(n + (los ? 1 : 0), m)));
//            weighted_sumPowers += weight;
//        }
//    }
//
//    // Calculate angular spreads
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::sqrt(-2.0 * std::log(std::abs(weighted_sumAOD / weighted_sumPowers))));
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::sqrt(-2.0 * std::log(std::abs(weighted_sumAOA / weighted_sumPowers))));
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::sqrt(-2.0 * std::log(std::abs(weighted_sumZOD / weighted_sumPowers))));
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::sqrt(-2.0 * std::log(std::abs(weighted_sumZOA / weighted_sumPowers))));
//
//    // Calculate mean angles
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::atan2(weighted_sumAOD.imag(), weighted_sumAOD.real()));
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::atan2(weighted_sumAOA.imag(), weighted_sumAOA.real()));
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::atan2(weighted_sumZOD.imag(), weighted_sumZOD.real()));
//    ASandMeanAnglesforAOD_AOA_ZOD_ZOA.push_back(std::atan2(weighted_sumZOA.imag(), weighted_sumZOA.real()));
//
//    return ASandMeanAnglesforAOD_AOA_ZOD_ZOA;
//}

// Mixing angle:

// Функция для перемешивания подгруппы в строке 
void SSP::mixingClusterInRow(Eigen::VectorXd& row, const std::vector<size_t>& indices, std::minstd_rand& g) {
    for (size_t i = indices.size() - 1; i > 0; --i) {
        size_t j = std::uniform_int_distribution<size_t>(0, i)(g);
        std::swap(row(indices[i]), row(indices[j]));
    }
}
// Функция для перемешивания подгрупп во всех строках матрицы
void SSP::mixingMatrix(Eigen::MatrixXd& matrix, std::minstd_rand& g) {
#pragma omp parallel for
    for (int i = 0; i < matrix.rows(); ++i) {
        Eigen::VectorXd row = matrix.row(i);
        if (i < 2) {
            mixingClusterInRow(row, subCluster1, g);
            mixingClusterInRow(row, subCluster2, g);
            mixingClusterInRow(row, subCluster3, g);
        }
        else { mixingClusterInRow(row, cluster, g); }

        matrix.row(i) = row;
    }
}
// Основной метод для перемешивания подгрупп в структуре LinkData
void SSP::mixingAngles(LinkData& data, std::minstd_rand& g) {
#pragma omp parallel sections
    {
#pragma omp section
        mixingMatrix(data.AOD_n_m, g); // Обрабатываем AOD_n_m
#pragma omp section
        mixingMatrix(data.AOA_n_m, g); // Обрабатываем AOA_n_m
#pragma omp section
        mixingMatrix(data.ZOD_n_m, g); // Обрабатываем ZOD_n_m
#pragma omp section
        mixingMatrix(data.ZOA_n_m, g); // Обрабатываем ZOA_n_m
    }
}




// power -> max_min
void SSP::sortRelativeToFirstVector(
    std::vector<double>& mainVector,
    std::vector<double>& secondVector,
    Eigen::MatrixXd& matrix1,
    Eigen::MatrixXd& matrix2,
    Eigen::MatrixXd& matrix3,
    Eigen::MatrixXd& matrix4,
    Eigen::MatrixXd& matrix5
) {

    std::vector<size_t> indices(mainVector.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(), [&mainVector](size_t i, size_t j) {
        return mainVector[i] > mainVector[j];
        });


    std::vector<double> sortedMainVector(mainVector.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        sortedMainVector[i] = mainVector[indices[i]];
    }
    mainVector = sortedMainVector;


    std::vector<double> sortedSecondVector(secondVector.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        sortedSecondVector[i] = secondVector[indices[i]];
    }
    secondVector = sortedSecondVector;

    auto sortMatrixRows = [&indices](Eigen::MatrixXd& mat) {
        Eigen::MatrixXd sortedMat(mat.rows(), mat.cols());
        for (size_t i = 0; i < indices.size(); ++i) {
            sortedMat.row(i) = mat.row(indices[i]);
        }
        mat = sortedMat;
        };


    sortMatrixRows(matrix1);
    sortMatrixRows(matrix2);
    sortMatrixRows(matrix3);
    sortMatrixRows(matrix4);
    sortMatrixRows(matrix5);

    // Удаление нулевых элементов в конце главного вектора
    auto it = std::find_if(mainVector.rbegin(), mainVector.rend(), [](double val) {
        return val != 0.0;
        });

    if (it != mainVector.rend()) {
        // Находим индекс последнего ненулевого элемента
        size_t lastNonZeroIndex = std::distance(mainVector.begin(), it.base()) - 1;

        // Укорачиваем главный вектор
        mainVector.resize(lastNonZeroIndex + 1);

        // Укорачиваем второй вектор
        secondVector.resize(lastNonZeroIndex + 1);

        // Укорачиваем матрицы (удаляем строки)
        matrix1.conservativeResize(lastNonZeroIndex + 1, matrix1.cols());
        matrix2.conservativeResize(lastNonZeroIndex + 1, matrix2.cols());
        matrix3.conservativeResize(lastNonZeroIndex + 1, matrix3.cols());
        matrix4.conservativeResize(lastNonZeroIndex + 1, matrix4.cols());
        matrix5.conservativeResize(lastNonZeroIndex + 1, matrix5.cols());
    }
}