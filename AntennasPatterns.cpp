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



