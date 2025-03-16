#ifndef TABLES_H
#define TABLES_H

#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <vector>

using namespace Eigen;

MatrixXd GenerateMatrix(bool los, double f, std::string scenario);
double LinearInterpolation(double x, int x_0, int x_1, double y_0, double y_1);

VectorXd InterpolatedParameters(Eigen::MatrixXd Table, double angle);

#endif 