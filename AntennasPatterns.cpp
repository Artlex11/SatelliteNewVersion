#include "AntennasPatterns.h"

namespace NTN
{
	double AntennaPattern_NTN1_V(double el_deg)
	{
		el_deg > 180 ? el_deg = el_deg - 360 : el_deg;
		el_deg < -180 ? el_deg = el_deg + 360 : el_deg;

		//std::cout << "el_deg_NTN_V: " << el_deg << std::endl;

		double AP_V_dB = -std::min(12 * pow((el_deg - 90) / phi_3dB, 2), AP_max);

		return AP_V_dB;
	}

	double AntennaPattern_NTN1_H(double az_deg)
	{
		az_deg > 180 ? az_deg = az_deg - 360 : az_deg;
		az_deg < -180 ? az_deg = az_deg + 360 : az_deg;

		//std::cout << "az_deg_NTN_H: " << az_deg << std::endl;

		double AP_H_dB = -std::min(12 * pow(az_deg / phi_3dB, 2), AP_max);

		return AP_H_dB;
	}

	double AntennaPattern_NTN1_3D(double az_deg, double el_deg)
	{
		double AP_3D_dB = AP_gain - std::min(-AntennaPattern_NTN1_V(el_deg) + AntennaPattern_NTN1_H(az_deg), AP_max);
		return AP_3D_dB;
	}
}

namespace INDOOR
{
	double AntennaPattern_Indoor_V(double el_deg)
	{
		el_deg > 180 ? el_deg = el_deg - 360 : el_deg;
		el_deg < -180 ? el_deg = el_deg + 360 : el_deg;

		//std::cout << "el_deg_INDOOR: " << el_deg << std::endl;

		double AP_V_dB = -std::min(12 * pow((el_deg - 90) / phi_3dB, 2), AP_max);

		return AP_V_dB;
	}

	double AntennaPattern_Indoor_H(double az_deg)
	{
		az_deg > 180 ? az_deg = az_deg - 360 : az_deg;
		az_deg < -180 ? az_deg = az_deg + 360 : az_deg;

		//std::cout << "az_deg_INDOOR: " << az_deg << std::endl;

		double AP_H_dB = -std::min(12 * pow((az_deg - 90) / phi_3dB, 2), AP_max);

		return AP_H_dB;
	}

	double AntennaPattern_Indoor_3D(double az_deg, double el_deg)
	{
		double AP_3D_dB = AP_gain - std::min(-AntennaPattern_Indoor_V(el_deg) + AntennaPattern_Indoor_H(az_deg), AP_max);
		return AP_3D_dB;
	}
}

// Функция для вычисления диаграммы направленности
//void calculateDishPattern(Eigen::VectorXd& teta_rad, Eigen::VectorXd& AP_dB, double rDish_WL)
//{
//	double Gain = 10 * log10(std::pow(M_PI * 2 * rDish_WL, 2)) - 2.4478;
//
//	for (int i = 0; i < teta_rad.size(); ++i)
//	{
//		double teta = teta_rad(i);
//		if (teta == 0)
//		{
//			AP_dB(i) = 1.0;
//		}
//		else
//		{
//			double bessel_arg = 2 * M_PI * rDish_WL * sin(teta);
//			double bessel_val = std::abs(std::cyl_bessel_j(1, bessel_arg));
//			//double bessel_val = std::cyl_bessel_j(1, bessel_arg);
//			double res = 4 * std::pow(bessel_val / bessel_arg, 2);
//			//double res = 4 * std::pow(std::abs(bessel_val / bessel_arg), 2);
//			AP_dB(i) = 10 * log10(res) + Gain;
//			//AP_lin(i) = pow(10, AP_dB(i)); - linear scale
//		}
//	}
//}


