//#ifndef NTN_DEPLOYMENT_H
//#define NTN_DEPLOYMENT_H
//
//#include <iostream>
//#include <random>
//#include <cmath>
//#include <complex>
//#include <Eigen/Dense>
//#include <vector>
//#include "Generators.h"
//
//using namespace Eigen;
//
//// Глобальные константы
//const double EARTH_RADIUS = 6371000.0; // Радиус Земли в метрах
//
//const double PI = 3.14159265358979323846; // Число Пи
//
//Eigen::MatrixXd rotateAroundAxis(const Eigen::MatrixXd& inVect, const Eigen::Vector3d& ax, double angl); // Поворот Родрига
//
//class Satellite 
//{
//    private:
//        double x, y, z; // Координаты спутника
//        const double beamWidth; // Ширина луча
//        double elMinDegrees = 10.0; // Минимальный угол ( места ) видимости 
//        double elTargetDegrees = 90.0; // Угол цели спутника 
//        Eigen::MatrixXf uvSet; // Матрица для хранения UV-координат
//
//    public:
//        Satellite(double x, double y, double z, double beamWidth)
//            : x(x), y(y), z(z + EARTH_RADIUS), beamWidth(beamWidth) {
//        };
//
//        double getBeamWidth() const;
//
//
//        void printCoordinates() const;
//
//        // метод для расчета угла места
//        void calculateZenithAngle(Eigen::MatrixXd& users);
//
//        Eigen::MatrixXd dropUEs(int nTiers, int nUePerCell);
//};
//#endif