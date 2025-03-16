#ifndef LOS_PROBABILITY_H
#define LOS_PROBABILITY_H

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


double AngleForLSP(double deg);

bool CalculateLOSProbability(int index, std::string scenario);

#endif 