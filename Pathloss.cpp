#include "PathLoss.h"

// PL = PL_b + PL_g + PL_s + PL_e
// PL = PL_b + PL_g
// PL_b = FSPL + SF + CL(alpha, f)


// Re - ������ �����
// h0 - ������ ��������
// alpha - elevation angle

// SF_std

// LOS

// S-band
Vector<double, 9> SF_std_Dense_Urban_LOS_S = { 3.5, 3.4, 2.9, 3.0, 3.1, 2.7, 2.5, 2.3, 1.2 };
Vector<double, 9> SF_std_Urban_LOS_S = { 4, 4, 4, 4, 4, 4, 4, 4, 4 };
Vector<double, 9> SF_std_Suburban_LOS_S = { 1.79, 1.14, 1.14, 0.92, 1.42, 1.56, 0.85, 0.72, 0.72 };
Vector<double, 9> SF_std_Rural_LOS_S = { 1.79, 1.14, 1.14, 0.92, 1.42, 1.56, 0.85, 0.72, 0.72 };

// Ka-band
Vector<double, 9> SF_std_Dense_Urban_LOS_Ka = { 2.9, 2.4, 2.7, 2.4, 2.4, 2.7, 2.6, 2.8, 0.6 };
Vector<double, 9> SF_std_Urban_LOS_Ka = { 4, 4, 4, 4, 4, 4, 4, 4, 4 };
Vector<double, 9> SF_std_Suburban_LOS_Ka = { 1.9, 1.6, 1.9, 2.3, 2.7, 3.1, 3.0, 3.6, 0.4 };
Vector<double, 9> SF_std_Rural_LOS_Ka = { 1.9, 1.6, 1.9, 2.3, 2.7, 3.1, 3.0, 3.6, 0.4 };

// NLOS

// S-band
Vector<double, 9> SF_std_Dense_Urban_NLOS_S = { 15.5, 13.9, 12.4, 11.7, 10.6, 10.5, 10.1, 9.2, 9.2 };
Vector<double, 9> SF_std_Urban_NLOS_S = { 6, 6, 6, 6, 6, 6, 6, 6, 6 };
Vector<double, 9> SF_std_Suburban_NLOS_S = { 8.93, 9.08, 8.78, 10.25, 10.56, 10.74, 10.17, 11.52, 11.52 };
Vector<double, 9> SF_std_Rural_NLOS_S = { 8.93, 9.08, 8.78, 10.25, 10.56, 10.74, 10.17, 11.52, 11.52 };

// Ka-band
Vector<double, 9> SF_std_Dense_Urban_NLOS_Ka = { 17.1, 17.1, 15.6, 14.6, 14.2, 12.6, 12.1, 12.3, 12.3 };
Vector<double, 9> SF_std_Urban_NLOS_Ka = { 6, 6, 6, 6, 6, 6, 6, 6, 6 };
Vector<double, 9> SF_std_Suburban_NLOS_Ka = { 10.7, 10.0, 11.2, 11.6, 11.8, 10.8, 10.8, 10.8, 10.8 };
Vector<double, 9> SF_std_Rural_NLOS_Ka = { 10.7, 10.0, 11.2, 11.6, 11.8, 10.8, 10.8, 10.8, 10.8 };


// ClusterLoss (for NLOS case)

// S-band
Vector<double, 9> CL_Dense_Urban_S = { 34.3, 30.9, 29.0, 27.7, 26.8, 26.2, 25.8, 25.5, 25.5 };
Vector<double, 9> CL_Urban_S = { 34.3, 30.9, 29.0, 27.7, 26.8, 26.2, 25.8, 25.5, 25.5 };
Vector<double, 9> CL_Suburban_S = { 19.52, 18.17, 18.42, 18.28, 18.63, 17.68, 16.50, 16.30, 16.30 };
Vector<double, 9> CL_Rural_S = { 19.52, 18.17, 18.42, 18.28, 18.63, 17.68, 16.50, 16.30, 16.30 };

// Ka-band
Vector<double, 9> CL_Dense_Urban_Ka = { 44.3, 39.9, 37.5, 35.8, 34.6, 33.8, 33.3, 33.0, 32.9 };
Vector<double, 9> CL_Urban_Ka = { 44.3, 39.9, 37.5, 35.8, 34.6, 33.8, 33.3, 33.0, 32.9 };
Vector<double, 9> CL_Suburban_Ka = { 29.5, 24.6, 21.9, 20.0, 18.7, 17.8, 17.2, 16.9, 16.8 };
Vector<double, 9> CL_Rural_Ka = { 29.5, 24.6, 21.9, 20.0, 18.7, 17.8, 17.2, 16.9, 16.8 };

// ������ ������� ������
double CalculateDistance(const double EARTH_RADIUS, double altitude, double alpha)
{
	double d = sqrt(EARTH_RADIUS * EARTH_RADIUS * sin(alpha) * sin(alpha) + altitude * altitude + 2 * altitude * EARTH_RADIUS) - EARTH_RADIUS * sin(alpha);
	return d;
}


double ChooseSTD(bool los, double f, double alpha, std::string scenario)
{
	int deg = int(alpha / 10);
	double std;
	if (los)
	{
		if (f < 6.0)
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_LOS_S[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_LOS_S[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_LOS_S[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_LOS_S[deg];
			}
		}
		else
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_LOS_Ka[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_LOS_Ka[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_LOS_Ka[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_LOS_Ka[deg];
			}
		}
	}
	else
	{
		if (f < 6.0)
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_NLOS_S[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_NLOS_S[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_NLOS_S[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_NLOS_S[deg];
			}
		}
		else
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_NLOS_Ka[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_NLOS_Ka[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_NLOS_Ka[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_NLOS_Ka[deg];
			}
		}
	}
	return std;
}


double GenerateSF(double std)
{
	double SF = RandomGenerators::generateGauss(0, std);
	return SF;
}


double Calculate_FSPL(double d, double f)
{
	double FSPL = 32.45 + 20 * log10(f) + 20 * log10(d);
	return FSPL;
}


double ChooseCL(bool los, double f, double alpha, std::string scenario)
{
	int deg = int(alpha / 10);
	double CL;
	if (los)
	{
		CL = 0;
	}
	else
	{
		if (f < 6.0)
		{
			if (scenario == "DenseUrban")
			{
				CL = CL_Dense_Urban_S[deg];
			}
			else if (scenario == "Urban")
			{
				CL = CL_Urban_S[deg];
			}
			else if (scenario == "Suburban")
			{
				CL = CL_Suburban_S[deg];
			}
			else if (scenario == "Rural")
			{
				CL = CL_Rural_S[deg];
			}
		}
		else
		{
			if (scenario == "DenseUrban")
			{
				CL = CL_Dense_Urban_Ka[deg];
			}
			else if (scenario == "Urban")
			{
				CL = CL_Urban_Ka[deg];
			}
			else if (scenario == "Suburban")
			{
				CL = CL_Suburban_Ka[deg];
			}
			else if (scenario == "Rural")
			{
				CL = CL_Rural_Ka[deg];
			}
		}
	}
	return CL;
}


double CalculateBasisPathLoss(double FSPL, double SF, double CL)
{
	double PL_b = FSPL; // SF - ����������� �����  dist(0, 0.72) ,  CL  ���� �� ��������� ;
	return PL_b;
}

