//#ifndef NTN_DEPLOYMENT_H
//#define NTN_DEPLOYMENT_H
//
//#include <vector>
//#include <Eigen/Dense>
//#include <random>
//#include <iostream>
//#include <cmath>
//
//#include <omp.h> 
//
//#include "Generators.h"
//#include "Links.h"
//
//
//#ifndef M_PI
//#define M_PI 3.14159265358979323846
//#endif
//
//
//class SatelliteLink {
//
//
//public:
//    SatelliteLink(int nTiersCore, int nTiresWrArnd, int nUePerCell, double satHeightKm, double elMinDegrees, double elTargetDegrees);
//    void generateLinks();
//    void transformCoordinates(double gamma);
//    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d >> links;
//    Links links;
//    std::vector<Eigen::Vector3d> getVectorsToCellCenters();
//    Eigen::Vector3d uvTo3D(const Eigen::Vector2d& uv);
//
//private:
//    void generateUVPlane();
//    double calculateElevetionAngle(const Eigen::Vector3d& UE, const Eigen::Vector3d& Sat);
//    std::vector<Eigen::Vector3d> generateRandomPoints(int count);
//
//    int nTiersCore, nTiresWrArnd, nUePerCell, nTiers, nCells, nUEs;
//    double satHeightKm, elMinDegrees, elTargetDegrees, rEarth, beamWidth_degrees, uvStep, uvBeamRadius;
//    Eigen::Vector3d xyzSat, xyzSatNew;
//    std::vector<Eigen::Vector2d> uvSet;
//
//};
//
//#endif // SATELLITE_LINK_H