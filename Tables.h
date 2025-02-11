#ifndef TABLES_H
#define TABLES_H

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <vector>

using namespace Eigen;

MatrixXd GenerateMatrix(bool los, double f, std::string scenario);

#endif 