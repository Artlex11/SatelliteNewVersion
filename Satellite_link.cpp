#include "Satellite_link.h"



SatelliteLink::SatelliteLink(int nTiersCore, int nTiresWrArnd, int nUePerCell, double satHeightKm, double elMinDegrees, double elTargetDegrees, double azTargetDegrees)
    : nTiersCore(nTiersCore), nTiresWrArnd(nTiresWrArnd), nUePerCell(nUePerCell), satHeightKm(satHeightKm), elMinDegrees(elMinDegrees), elTargetDegrees(elTargetDegrees), azTargetDegrees(azTargetDegrees) {
    nTiers = nTiersCore + nTiresWrArnd;
    nCells = 1 + (nTiers * (nTiers + 1) * 3);
    nUEs = nCells * nUePerCell;

    rEarth = 6371;
    beamWidth_degrees = 4.4127;
    uvStep = std::sin(M_PI / 180 * (beamWidth_degrees * 0.865));
    uvBeamRadius = uvStep / std::sqrt(3);

    xyzSat = Eigen::Vector3d(1.0, -0.05, 0.0);
    /*xyzSat << RandomGenerators::generateGauss(0.0, 1.0) ,
       RandomGenerators::generateGauss(0.0, 1.0),
       RandomGenerators::generateGauss(0.0, 1.0);*/

    xyzSat.normalize();
    xyzSatOnEarth = xyzSat;
    xyzSat *= (rEarth + satHeightKm);
}


void SatelliteLink::generateLinks() {
    std::vector<Eigen::Vector3d> xyzUEs_all(nUEs);
    int UEcnt = 0;

    generateUVPlane();

    double phiMax = M_PI - (M_PI / 180.0 * (90.0 + elMinDegrees)) - std::asin(std::sin(M_PI / 180.0 * (90.0 + elMinDegrees)) * rEarth / (rEarth + satHeightKm));
    double phiSat = M_PI - (M_PI / 180.0 * (90.0 + elTargetDegrees)) - std::asin(std::sin(M_PI / 180.0 * (90.0 + elTargetDegrees)) * rEarth / (rEarth + satHeightKm));
    double thetaSat = M_PI - phiSat - (M_PI / 180.0 * (90.0 + elTargetDegrees));
    double treshold = rEarth * std::cos(phiMax);

    // ��������� ����������� ������� ������� ������������ ���������� ��� Z
    Eigen::Vector3d globalZ(0, 0, 1);

    // ���� ������� �� ����������� � ���������� ���� Z, ���������� � ��� ���������� p1
    if (!xyzSat.isApprox(globalZ) && !xyzSat.isApprox(-globalZ)) {
        p1 = xyzSat.cross(globalZ).normalized(); // ������ ������������� ������
    }
    else {
        // ���� ������� ����������� � ���� Z, ���������� ���������� ��� X
        p1 = Eigen::Vector3d(1, 0, 0);
    }
    p2 = xyzSat.cross(p1).normalized();

    Eigen::Vector3d rotationAxis = (std::cos(azTargetDegrees * M_PI / 180.0) * p1 + std::sin(azTargetDegrees * M_PI / 180.0) * p2).normalized();


    Eigen::Matrix3d rotMatrix = Eigen::AngleAxisd(thetaSat, rotationAxis).toRotationMatrix();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);



    while (UEcnt < nUEs) {
        std::vector<Eigen::Vector3d> xyzUEs = generateRandomPoints(nUEs * 10);

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(xyzUEs.size()); ++i) {
            const auto& xyzUE = xyzUEs[i];
            if (xyzUE.dot(xyzSatOnEarth * rEarth) / rEarth <= treshold) continue;

            Eigen::Vector3d rSatUe = (xyzUE - xyzSat).normalized();

            if (thetaSat != 0.0) {
                rSatUe = rotMatrix * rSatUe;
            }

            Eigen::Vector2d uvUE(rSatUe.dot(p1), rSatUe.dot(p2)); // u - rSatUe.dot(p1) , v - rSatUe.dot(p2)


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

        links.addLink((LinkData{ user, xyzSat, SatelliteLink::calculateElevetionAngle(user,xyzSat), Antenna{}, Antenna{} }));
    }

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
    for (int i = 0; i < count; ++i) { // ����������: int ������ size_t
        points[i] = Eigen::Vector3d(d(gen), d(gen), d(gen)).normalized() * rEarth;
    }

    return points;
}


//-------------------------�����------------------------------------------------//

std::vector<Eigen::Vector3d> SatelliteLink::getVectorsToCellCenters() {
    std::vector<Eigen::Vector3d> vectorsToCenters;

    // �������� ������ ����� � UV-���������
    const std::vector<Eigen::Vector2d>& uvCenters = uvSet;

    // ����������� UV-���������� � ��������� ����������
    for (const auto& uv : uvCenters) {
        // ����������� UV-���������� � ��������� ���������� �� ����������� �����
        Eigen::Vector3d centerOnEarth = uvTo3D(uv);

        // ��������� ������ �� �������� �� ������ ������
        Eigen::Vector3d vectorToCenter = -centerOnEarth;
        vectorsToCenters.push_back(vectorToCenter.normalized()); // ����������� ������
    }

    return vectorsToCenters;
}

Eigen::Vector3d SatelliteLink::uvTo3D(const Eigen::Vector2d& uv) {
    // ����������� UV-���������� � ��������� ���������� �� ����������� �����
    double u = uv(0);
    double v = uv(1);

    // ���������� �������� ������� p1 � p2, ������� ��� ���������� � generateLinks()
    Eigen::Vector3d centerOnEarth = xyzSatOnEarth + u * p1 + v * p2;

    // ����������� � �������� �� ������ �����, ����� ����� ������ �� �����������
    centerOnEarth.normalize();
    centerOnEarth *= rEarth;

    return centerOnEarth;
}


//std::vector<Eigen::Vector3d> SatelliteLink::generateRandomPoints(int count) {
//    std::vector<Eigen::Vector3d> points(count);
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> disPhi(0, 2 * M_PI);
//    std::uniform_real_distribution<> disTheta(0, M_PI);
//
//#pragma omp parallel for
//    for (int i = 0; i < count; ++i) {
//        double phi = disPhi(gen);
//        double theta = disTheta(gen);
//        double x = sin(theta) * cos(phi);
//        double y = sin(theta) * sin(phi);
//        double z = cos(theta);
//        points[i] = Eigen::Vector3d(x, y, z) * rEarth;
//    }
//
//    return points;
//}




// ������ �� �� ������� ���� � �� �� ����� , � �������� ��������� ��� ������ 
// 
//void SatelliteLink::transformCoordinates(double gamma) {
//    Eigen::Vector3d d = xyzSatOnEarth;
//
//    double alpha = std::atan2(d.y(), d.x());
//    double beta = -M_PI / 2 + std::acos(d.z() / d.norm());
//
//    Eigen::Matrix3d Rz;
//    Rz << std::cos(alpha), -std::sin(alpha), 0,
//        std::sin(alpha), std::cos(alpha), 0,
//        0, 0, 1;
//
//    Eigen::Matrix3d Ry;
//    Ry << std::cos(beta), 0, std::sin(beta),
//        0, 1, 0,
//        -std::sin(beta), 0, std::cos(beta);
//
//    Eigen::Matrix3d Rx;
//    Rx << 1, 0, 0,
//        0, std::cos(gamma), -std::sin(gamma),
//        0, std::sin(gamma), std::cos(gamma);
//
//    Eigen::Matrix3d Rzyx = Rz * Ry * Rx;
//
//    // �������� ������ �� ������ ������
//    std::vector<LinkData>& linkList = links.getLinks();
//
//#pragma omp parallel for
//    for (int i = 0; i < static_cast<int>(linkList.size()); ++i) {
//        linkList[i].userPosition = Rzyx * linkList[i].userPosition;
//        linkList[i].satellitePosition = Rzyx * linkList[i].satellitePosition;
//    }
//
//
//
////#pragma omp parallel for
//    //for (auto& link : links.getLinks()) {
//    //    link.userPosition = Rzyx * link.userPosition;
//    //    link.satellitePosition = Rzyx * link.satellitePosition;
//    //}
//}