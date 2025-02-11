#ifndef PATHLOSS_H
#define PATHLOSS_H
#include <iostream>
//#include "NTN_Deployment.h"
#include "Tables.h"
#include "Generators.h"

double CalculateDistance(double Re, double h0, double alpha);

double GenerateSF(double std);

double CalculatePathLoss(double d, double f);

#endif

