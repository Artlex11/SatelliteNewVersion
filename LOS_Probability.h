#ifndef CALCULATE_PROBABILITY_LOS_H
#define CALCULATE_PROBABILITY_LOS_H

#include <iostream>
#include <cmath>
#include "Generators.h"
#include <Eigen/Dense>
#include <string>

using namespace Eigen;

// Вектора с вероятностями LOS для разных сценариев
extern Vector<double, 9> DenseUrban;
extern Vector<double, 9> Urban;
extern Vector<double, 9> Suburban_Rural;

extern Vector<double, 9> angles;

// Массивы double для вероятностей
//double DenseUrban[] = {28.2, 33.1, 39.8, 46.8, 53.7, 61.2, 73.8, 82.0, 98.1};
//double Urban[] = {24.6, 38.6, 49.3, 61.3, 72.6, 80.5, 91.9, 96.8, 99.2};
//double Suburban_Rural[] = {78.2, 86.9, 91.9, 92.9, 93.5, 94.0, 94.9, 95.2, 99.8};
//
//double angles[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

double AngleForLSP(double deg);

bool CalculateLOSProbability(int index, std::string scenario);

#endif #pragma once