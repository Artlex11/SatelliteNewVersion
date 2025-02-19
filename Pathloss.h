#ifndef PATHLOSS_H
#define PATHLOSS_H
#include <iostream>
#include "Tables.h"
#include "Generators.h"

double CalculateDistance(double Re, double h0, double alpha);

double GenerateSF(double std);

double CalculatePathLoss(double d, double f);

#endif

