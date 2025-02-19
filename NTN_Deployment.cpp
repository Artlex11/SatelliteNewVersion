#include "NTN_Deployment.h"



SatelliteLink::SatelliteLink(int nTiersCore, int nTiresWrArnd, int nUePerCell, double satHeightKm, double elMinDegrees, double elTargetDegrees)
    : nTiersCore(nTiersCore), nTiresWrArnd(nTiresWrArnd), nUePerCell(nUePerCell), satHeightKm(satHeightKm), elMinDegrees(elMinDegrees), elTargetDegrees(elTargetDegrees) {
    nTiers = nTiersCore + nTiresWrArnd;
    nCells = 1 + (nTiers * (nTiers + 1) * 3);
    nUEs = nCells * nUePerCell;

    rEarth = 6371;
    beamWidth_degrees = 4.4127;
    uvStep = std::sin(M_PI / 180 * (beamWidth_degrees));
    uvBeamRadius = uvStep / std::sqrt(3);

    xyzSat = Eigen::Vector3d(rEarth + satHeightKm, 0.0, 0.0);
    xyzSatNew << RandomGenerators::generateGauss(0.0, 1.0) * (rEarth + satHeightKm),
        RandomGenerators::generateGauss(0.0, 1.0)* (rEarth + satHeightKm),
        RandomGenerators::generateGauss(0.0, 1.0)* (rEarth + satHeightKm);

    xyzSatNew.normalize();
    xyzSatNew *= (rEarth + satHeightKm);
}

void SatelliteLink::generateLinks() {
    std::vector<Eigen::Vector3d> xyzUEs_all(nUEs);
    int UEcnt = 0;

    generateUVPlane();

    double phiMax = M_PI - (M_PI / 180.0 * (90.0 + elMinDegrees)) - std::asin(std::sin(M_PI / 180.0 * (90.0 + elMinDegrees)) * rEarth / (rEarth + satHeightKm));
    double phiSat = M_PI - (M_PI / 180.0 * (90.0 + elTargetDegrees)) - std::asin(std::sin(M_PI / 180.0 * (90.0 + elTargetDegrees)) * rEarth / (rEarth + satHeightKm));
    double thetaSat = M_PI - phiSat - (M_PI / 180.0 * (90.0 + elTargetDegrees));
    double treshold = rEarth * std::cos(phiMax);

    Eigen::Matrix3d rotMatrix = Eigen::AngleAxisd(thetaSat, Eigen::Vector3d::UnitZ()).toRotationMatrix();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);

    while (UEcnt < nUEs) {
        std::vector<Eigen::Vector3d> xyzUEs = generateRandomPoints(nUEs * 10);

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(xyzUEs.size()); ++i) { // Исправлено: size_t -> int
            const auto& xyzUE = xyzUEs[i];
            if (xyzUE.x() <= treshold) continue;

            Eigen::Vector3d rSatUe = (xyzSat - xyzUE).normalized();

            if (thetaSat != 0.0) {
                rSatUe = rotMatrix * rSatUe;
            }

            Eigen::Vector2d uvUE(rSatUe(1), rSatUe(2));

            for (const auto& uvCell : uvSet) {
                if ((uvUE - uvCell).squaredNorm() <= uvBeamRadius * uvBeamRadius) {
#pragma omp critical
                    {
                        if (UEcnt < nUEs) {
                            xyzUEs_all[UEcnt++] = xyzUE;
                        }
                    }
                    break;
                }
            }
        }
    }

    for (const auto& user : xyzUEs_all) {
        links.addLink(std::move(LinkData{ user, xyzSat, SatelliteLink::calculateElevetionAngle(user,xyzSat), Antenna{}, Antenna{} }));
    }

    transformCoordinates(0.0);
}

void SatelliteLink::generateUVPlane() {
#pragma omp parallel for
    for (int i = -nTiers; i <= nTiers; ++i) {
        for (int j = -nTiers; j <= nTiers; ++j) {
            double u = i * uvStep + j * (uvStep / 2);
            double v = j * (std::sqrt(3) * uvStep / 2);
            if (u <= (nTiers * uvStep - v / std::sqrt(3)) + 1e-10 && u >= (-nTiers * uvStep - v / std::sqrt(3)) - 1e-10) {
#pragma omp critical
                {
                    uvSet.push_back(Eigen::Vector2d(u, v));
                }
            }
        }
    }
}

double SatelliteLink::calculateElevetionAngle(const Eigen::Vector3d& UE, const Eigen::Vector3d& Sat)
{

    return std::asin((Sat - UE).normalized().dot(UE.normalized())) * 180 / M_PI;
}

std::vector<Eigen::Vector3d> SatelliteLink::generateRandomPoints(int count) {
    std::vector<Eigen::Vector3d> points(count);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);

#pragma omp parallel for
    for (int i = 0; i < count; ++i) { // Исправлено: int вместо size_t
        points[i] = Eigen::Vector3d(std::abs(d(gen)), d(gen), d(gen)).normalized() * rEarth;
    }

    return points;
}

void SatelliteLink::transformCoordinates(double gamma) {
    Eigen::Vector3d d = xyzSatNew;

    double alpha = std::atan2(d.y(), d.x());
    double beta = -M_PI / 2 + std::acos(d.z() / d.norm());

    Eigen::Matrix3d Rz;
    Rz << std::cos(alpha), -std::sin(alpha), 0,
        std::sin(alpha), std::cos(alpha), 0,
        0, 0, 1;

    Eigen::Matrix3d Ry;
    Ry << std::cos(beta), 0, std::sin(beta),
        0, 1, 0,
        -std::sin(beta), 0, std::cos(beta);

    Eigen::Matrix3d Rx;
    Rx << 1, 0, 0,
        0, std::cos(gamma), -std::sin(gamma),
        0, std::sin(gamma), std::cos(gamma);

    Eigen::Matrix3d Rzyx = Rz * Ry * Rx;

    // Получаем ссылку на вектор связей
    std::vector<LinkData>& linkList = links.getLinks();

#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(linkList.size()); ++i) {
        linkList[i].userPosition = Rzyx * linkList[i].userPosition;
        linkList[i].satellitePosition = Rzyx * linkList[i].satellitePosition;
    }



    //#pragma omp parallel for
        //for (auto& link : links.getLinks()) {
        //    link.userPosition = Rzyx * link.userPosition;
        //    link.satellitePosition = Rzyx * link.satellitePosition;
        //}
}