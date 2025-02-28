#ifndef PATHLOSSCALCULATE_H
#define PATHLOSSCALCULATE_H
#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "Generators.h"
#include "Matlab_plot.h"

// для работы с MATLAB
#include "engine.h"
#include "mat.h"
// Глобальные константы
const double EARTH_RADIUS = 6371.0; // Радиус Земли

const double PI = 3.14159265358979323846; // Число Пи

using namespace Eigen;

// Рассчёт дистанции
double CalculateDistance(double Re, double h0, double alpha);

// Выбор СКО для SF
double ChooseSTD(bool los, double f, double alpha, std::string scenario);

// Генерация SF
double GenerateSF(double std);

// Free Space Pathloss
double Calculate_FSPL(double d, double f);

// Cluster Loss
double ChooseCL(bool los, double f, double alpha, std::string scenario);

// Basis Pathloss
double CalculateBasisPathLoss(double FSPL, double SF, double CL);

// Horizontal paths
double CalculateLh(double r, double s, double f, double t);

// Correlation for elevation angle of the path 
double CalculateLe(double elevation);

// Mu_1 and Sigma_1
Vector<double, 2> CalculateParametersForA(double Lh, double Le, double u, double v, double f);

// Mu_2 and Sigma_2
Vector<double, 2> CalculateParametersForB(double w, double x, double y, double z, double f);

// F(P)^(-1)
double CalculateInverseCummulativeNormalDistribution(double P);

// Building entry loss
double CalculateBuildingEntryLoss(double A, double B, double C);

//double CalculatePathLossInGasses(double d, double f);
#endif

//double CalculateDistance(double Re, double h0, double alpha);
//
//double GenerateSF(double std);
//
//double CalculatePathLoss(double d, double f);
