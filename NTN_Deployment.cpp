#include "NTN_Deployment.h"

SatelliteLink::SatelliteLink(int nTiersCore, int nTiresWrArnd, int nUePerCell, double satHeightKm, double elMinDegrees, double elTargetDegrees, double azTargetDegrees)
    : nTiersCore(nTiersCore), nTiresWrArnd(nTiresWrArnd), nUePerCell(nUePerCell), satHeightKm(satHeightKm), elMinDegrees(elMinDegrees), elTargetDegrees(elTargetDegrees), azTargetDegrees(azTargetDegrees) {
    nTiers = nTiersCore + nTiresWrArnd;
    nCells = 1 + (nTiers * (nTiers + 1) * 3);
    nUEs = nCells * nUePerCell;
    rEarth = 6371;
    beamWidth_degrees = 4.4127;
    uvStep = std::sin(M_PI / 180 * (beamWidth_degrees * 0.865));
    uvBeamRadius = uvStep / std::sqrt(3);


    //xyzSat = Eigen::Vector3d(1.0, 0.0, 0.0);
    xyzSat << RandomGenerators::generateGauss(0.0, 1.0),
        RandomGenerators::generateGauss(0.0, 1.0),
        RandomGenerators::generateGauss(0.0, 1.0);

    xyzSat.normalize();
    xyzSatOnEarth = xyzSat * rEarth;
    xyzSat = xyzSat * (rEarth + satHeightKm);
}



void SatelliteLink::generateLinks() {
    std::vector<Eigen::Vector3d> xyzUEs_all(nUEs);
    xyzUEs_all.reserve(nUEs);

    int UEcnt = 0;
    generateUVPlane();
    generateCoreUVPlane();

    const double phiMax = M_PI - (M_PI / 180.0 * (90.0 + elMinDegrees)) - std::asin(std::sin(M_PI / 180.0 * (90.0 + elMinDegrees)) * rEarth / (rEarth + satHeightKm));
    const double phiSat = M_PI - (M_PI / 180.0 * (90.0 + elTargetDegrees)) - std::asin(std::sin(M_PI / 180.0 * (90.0 + elTargetDegrees)) * rEarth / (rEarth + satHeightKm));
    const double thetaSat = M_PI - phiSat - (M_PI / 180.0 * (90.0 + elTargetDegrees));
    const double threshold = rEarth * std::cos(phiMax) * rEarth;
    const double uvBeamRadiusSquared = uvBeamRadius * uvBeamRadius;

    Eigen::Vector3d globalZ(0.0, 0.0, 1.0);


    if (!xyzSat.normalized().isApprox(globalZ) && !xyzSat.normalized().isApprox(-globalZ)) {
        p1 = xyzSat.cross(globalZ).normalized();
    }
    else {
        p1 = Eigen::Vector3d(0.0, -1.0, 0.0);
    }
    p2 = xyzSat.cross(p1).normalized();



    if (thetaSat != 90.0) {
        //const Eigen::Vector3d rotationAxis = (std::cos(azTargetDegrees * M_PI / 180.0) * p1 + std::sin(azTargetDegrees * M_PI / 180.0) * p2).normalized();
        rotMatrix = Eigen::AngleAxisd(thetaSat, (std::cos(azTargetDegrees * M_PI / 180.0) * p1 + std::sin(azTargetDegrees * M_PI / 180.0) * p2).normalized()).toRotationMatrix();
    }


    while (UEcnt < nUEs) {
        std::vector<Eigen::Vector3d> xyzUEs = generateRandomPoints(nUEs * 20);

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(xyzUEs.size()); ++i) {
            if (UEcnt >= nUEs) break; // Останавливаем, если набрали достаточно точек

            const auto& xyzUE = xyzUEs[i];
            const auto& projection = xyzUE.dot(xyzSatOnEarth);
            if (projection > threshold) {

                Eigen::Vector3d rSatUe = rotMatrix * ((xyzUE - xyzSat).normalized());

                const double u = rSatUe.dot(p1);
                const double v = rSatUe.dot(p2);


                for (const auto& uvCell : uvSet) {
                    const double du = u - uvCell.x();
                    const double dv = v - uvCell.y();
                    const double distanceSquared = du * du + dv * dv;

                    if (distanceSquared <= uvBeamRadiusSquared) {
#pragma omp critical
                        {
                            if (UEcnt < nUEs) {
                                xyzUEs_all[UEcnt++] = xyzUE;
                            }
                        }
                        break; // Прерываем цикл по uvSet
                    }
                }
            }
        }
    }

#pragma omp parallel for
    for (int i = 0; i < xyzUEs_all.size(); ++i) {
        const auto& user = xyzUEs_all[i];
#pragma omp critical
        {
            links.addLink({ user, xyzSat, std::asin((xyzSat - user).normalized().dot(user.normalized())) * 180 / M_PI });
        }
    }
}
//for (const auto& user : xyzUEs_all) {

//    links.addLink((LinkData{ user, xyzSat, std::asin((xyzSat - user).normalized().dot(user.normalized())) * 180 / M_PI })); // угол места std::asin((xyzSat - user).normalized().dot(user.normalized())) * 180 / M_PI
//}



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


void SatelliteLink::generateCoreUVPlane() {
#pragma omp parallel for
    for (int i = -nTiersCore; i <= nTiersCore; ++i) {
        for (int j = -nTiersCore; j <= nTiersCore; ++j) {
            double u = i * uvStep + j * (uvStep / 2);
            double v = j * (std::sqrt(3) * uvStep / 2);
            if (u <= (nTiersCore * uvStep - v / std::sqrt(3)) + 1e-10 && u >= (-nTiersCore * uvStep - v / std::sqrt(3)) - 1e-10) {
#pragma omp critical
                {
                    uvSetCore.push_back(Eigen::Vector2d(u, v));
                }
            }
        }
    }
}




std::vector<Eigen::Vector3d> SatelliteLink::generateRandomPoints(int count) {
    std::vector<Eigen::Vector3d> points(count);
    std::atomic<size_t> validPoints = 0;

#pragma omp parallel
    {
        std::random_device rd;
        std::minstd_rand gen(rd());
        std::normal_distribution<> d(0.0, 1.0);

        while (validPoints < static_cast<size_t>(count)) {
            Eigen::Vector3d point(d(gen), d(gen), d(gen));

            if (point.dot(xyzSat) > 0) {
                size_t idx = validPoints++;
                if (idx < static_cast<size_t>(count)) {
                    points[idx] = (point.normalized()) * rEarth;
                }
            }
        }
    }
    return points;
}

//
//std::vector<Eigen::Vector3d> SatelliteLink::getVectorsToCellCenters() {
//    std::vector<Eigen::Vector3d> vectorsToCenters;
//
//    // Получаем центры ячеек в UV-плоскости
//    const std::vector<Eigen::Vector2d>& uvCenters = uvSet;
//
//
//    for (const auto& uv : uvCenters) {
//
// 
//        Eigen::Vector3d centerCellOnEarth = uvTo3D(uv);
//        // Вычисляем вектор от спутника до центра ячейки
//        vectorsToCenters.push_back(centerCellOnEarth); // Нормализуем вектор
//    }
//
//    return vectorsToCenters;
//}
//
//Eigen::Vector3d SatelliteLink::uvTo3D(const Eigen::Vector2d& uv) {
//    // Преобразуем UV-координаты в трёхмерные координаты на поверхности Земли
//    Eigen::Vector3d xyzU = uv(0) * p1;
//    Eigen::Vector3d xyzV = uv(1) * p2;
//
//    // Используем базисные векторы p1 и p2, которые уже определены в generateLinks()
//    Eigen::Vector3d centerOnEarth = -(xyzSatOnEarth/rEarth) * (1 - (uv(0)* uv(0) + uv(1) * uv(1))) + xyzU + xyzV;
//
//    // Нормализуем и умножаем на радиус Земли, чтобы точка лежала на поверхности
//    centerOnEarth.normalize();
//
//    return centerOnEarth;
//}




std::vector<Eigen::Vector3d> SatelliteLink::getVectorsToCellCenters() {
    Eigen::Vector3d xyzSatOnEarthNorm = xyzSatOnEarth.normalized();
    std::vector<Eigen::Vector3d> vectorsToCenters;
    for (const auto& uv : uvSet) {
        Eigen::Vector3d xyzU = uv(0) * p1;
        Eigen::Vector3d xyzV = uv(1) * p2;


        vectorsToCenters.push_back((-(xyzSatOnEarth / rEarth) * (1 - (uv(0) * uv(0) + uv(1) * uv(1))) + xyzU + xyzV).normalized()); // Нормализуем вектор
    }
    return vectorsToCenters;
}