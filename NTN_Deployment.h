#ifndef SATELLITE_LINK_H
#define SATELLITE_LINK_H

#include <vector>
#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <cmath>
#include <unordered_set>


#include <omp.h> 

#include "Generators.h"
#include "Links.h"



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Пользовательская хэш-функция для Eigen::Vector2d
struct Vector2dHash {
    std::size_t operator()(const Eigen::Vector2d& v) const {
        std::size_t seed = 0;
        seed ^= std::hash<double>()(v.x()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= std::hash<double>()(v.y()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};



class SatelliteLink {


public:
    SatelliteLink(int nTiersCore, int nTiresWrArnd, int nUePerCell, double satHeightKm, double elMinDegrees, double elTargetDegrees, double azTargetDegrees);
    void generateLinks();
    Links links;

    std::vector<Eigen::Vector2d> uvSet;
    std::vector<Eigen::Vector2d> uvSetCore;

    Eigen::Matrix3d rotMatrix = Eigen::Matrix3d::Identity();

    Eigen::Vector3d p1, p2;
    std::vector<Eigen::Vector3d> getVectorsToCellCenters();
    Eigen::Vector3d uvTo3D(const Eigen::Vector2d& uv);


private:
    void generateUVPlane();
    void generateCoreUVPlane();
    //double calculateElevetionAngle(const Eigen::Vector3d& UE, const Eigen::Vector3d& Sat); // Встроен в generateLinks()
    std::vector<Eigen::Vector3d> generateRandomPoints(int count);

    int nTiersCore, nTiresWrArnd, nUePerCell, nTiers, nCells, nUEs;
    double satHeightKm, elMinDegrees, elTargetDegrees, azTargetDegrees, rEarth, beamWidth_degrees, uvStep, uvBeamRadius;
    Eigen::Vector3d xyzSat, xyzSatOnEarth;


};

#endif // SATELLITE_LINK_H