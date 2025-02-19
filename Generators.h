#ifndef GENERATORS_H
#define GENERATORS_H

#include <random>
#include <chrono>
#include <thread>
#include <vector>

namespace RandomGenerators {

    // ��������� ����������� �������������
    double generateGauss(double mean, double stddev);

    // ��������� ������������ ������������� ��� ����� �����
    int generateUniform(int min, int max);

    // ��������� ������������ ������������� ��� ������� �����
    double generateUniform(double min, double max);

    // ��������� ������������ ������������� �� ������� ��� ����� �����
    int generateUniformFromVector(const std::vector<int>& values);

    // ��������� ������������ ������������� �� ������� ��� ������� �����
    double generateUniformFromVector(const std::vector<double>& values);

}

#endif