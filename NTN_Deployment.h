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
//// ���������� ���������
//const double EARTH_RADIUS = 6371000.0; // ������ ����� � ������
//
//const double PI = 3.14159265358979323846; // ����� ��
//
//Eigen::MatrixXd rotateAroundAxis(const Eigen::MatrixXd& inVect, const Eigen::Vector3d& ax, double angl); // ������� �������
//
//class Satellite 
//{
//    private:
//        double x, y, z; // ���������� ��������
//        const double beamWidth; // ������ ����
//        double elMinDegrees = 10.0; // ����������� ���� ( ����� ) ��������� 
//        double elTargetDegrees = 90.0; // ���� ���� �������� 
//        Eigen::MatrixXf uvSet; // ������� ��� �������� UV-���������
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
//        // ����� ��� ������� ���� �����
//        void calculateZenithAngle(Eigen::MatrixXd& users);
//
//        Eigen::MatrixXd dropUEs(int nTiers, int nUePerCell);
//};
//#endif