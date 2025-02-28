#ifndef GENERATORS_H
#define GENERATORS_H

#include <random>
#include <chrono>
#include <thread>
#include <vector>

namespace RandomGenerators {

    // Генератор нормального распределения
    double generateGauss(double mean, double stddev);

    // Генератор равномерного распределения для целых чисел
    int generateUniform(int min, int max);

    // Генератор равномерного распределения для дробных чисел
    double generateUniform(double min, double max);

    // Генератор равномерного распределения из вектора для целых чисел
    int generateUniformFromVector(const std::vector<int>& values);

    // Генератор равномерного распределения из вектора для дробных чисел
    double generateUniformFromVector(const std::vector<double>& values);

}

#endif