#include "LOS_Probability.h"

// Вектора с вероятностями LOS для разных сценариев
Vector<double, 9> DenseUrban = { 28.2, 33.1, 39.8, 46.8, 53.7, 61.2, 73.8, 82.0, 98.1 };
Vector<double, 9> Urban = { 24.6, 38.6, 49.3, 61.3, 72.6, 80.5, 91.9, 96.8, 99.2 };
Vector<double, 9> Suburban_Rural = { 78.2, 86.9, 91.9, 92.9, 93.5, 94.0, 94.9, 95.2, 99.8 };

Vector<double, 9> angles = { 10, 20, 30, 40, 50, 60, 70, 80, 90 };

double AngleForLSP(double deg)
{
    int deg_div;
    double deg_mod;
    deg_div = int(deg / 10);
    deg_mod = int(deg) % 10;
    if (deg_mod < 5)
    {
        deg = deg_div * 10;
    }

    else
    {
        deg = (deg_div + 1) * 10;
    }

    if (deg < 5)
    {
        deg = 10;
    }

    return deg;
}

bool CalculateLOSProbability(int index, std::string scenario)
{
    double p = RandomGenerators::generateUniform(0.0, 1.0);
    double P_LOS;

    if (scenario == "Dense_Urban")
    {
        P_LOS = DenseUrban[index - 1] / 100;
    }
    else if (scenario == "Urban")
    {
        P_LOS = Urban[index - 1] / 100;
    }
    else
    {
        P_LOS = Suburban_Rural[index - 1] / 100;
    }

    std::cout << "P_LOS: " << P_LOS << ", p: " << p;
    bool los = (P_LOS > p);
    std::cout << "\nLink - " << (los ? "LOS" : "NLOS") << std::endl;

    return los;
}