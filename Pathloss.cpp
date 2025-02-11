#include "PathLoss.h"

// PL = PL_b + PL_g + PL_s + PL_e
// PL_b = FSPL + SF + CL(alpha, f)

double CalculateDistance(double Re, double h0, double alpha)
{
	double d = sqrt(pow(Re, 2) * pow(sin(alpha), 2) + pow(h0, 2) + 2 * h0 * Re) - Re * sin(alpha);
	return d;
}

double GenerateSF(double std)
{
	double SF = RandomGenerators::generateGauss(0, std);
	return SF;
}

double CalculatePathLoss(double d, double f, double SF, double CL)
{
	double PL = 0;
	double FSPL = 32.45 + 20 * log10(f) + 20 * log10(d);
	double PL_b = FSPL + SF + CL;
	double PL_e = 0;
	double PL_g = 0;
	double PL_s = 0;
	PL = PL_b + PL_g + PL_s + PL_e;

	return PL;
}

