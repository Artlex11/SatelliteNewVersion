#ifndef SATELLITE_LINK_H
#define SATELLITE_LINK_H

#include <vector>
#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <cmath>

#include <omp.h> 

#include "Generators.h"
#include "Links.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


class SatelliteLink {


public:
    SatelliteLink(int nTiersCore, int nTiresWrArnd, int nUePerCell, double satHeightKm, double elMinDegrees, double elTargetDegrees, double azTargetDegrees);
    void generateLinks();
    void transformCoordinates(double gamma);

    Links links;
    std::vector<Eigen::Vector2d> uvSet;

    //-------------------------Проба------------------------------------------------//
    Eigen::Vector3d p1, p2;
    std::vector<Eigen::Vector3d> getVectorsToCellCenters();
    Eigen::Vector3d uvTo3D(const Eigen::Vector2d& uv);


private:
    void generateUVPlane();
    double calculateElevetionAngle(const Eigen::Vector3d& UE, const Eigen::Vector3d& Sat);
    std::vector<Eigen::Vector3d> generateRandomPoints(int count);

    int nTiersCore, nTiresWrArnd, nUePerCell, nTiers, nCells, nUEs;
    double satHeightKm, elMinDegrees, elTargetDegrees, azTargetDegrees, rEarth, beamWidth_degrees, uvStep, uvBeamRadius;
    Eigen::Vector3d xyzSat, xyzSatOnEarth;
    //std::vector<Eigen::Vector2d> uvSet;

};

#endif // SATELLITE_LINK_H