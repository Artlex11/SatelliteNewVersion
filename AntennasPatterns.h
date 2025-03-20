#pragma once
#ifndef ANTENNAS_PATTERNS_H
#define ANTENNAS_PATTERNS_H

#include <iostream>
#include <cmath>

static double phi_3dB = 65;
static double AP_max = 25;

static double AP_gain = 4;

namespace NTN
{
	double AntennaPattern_NTN1_V(double el_deg);

	double AntennaPattern_NTN1_H(double az_deg);

	double AntennaPattern_NTN1_3D(double az_deg, double el_deg);
}

namespace INDOOR
{
	double AntennaPattern_Indoor_V(double el_deg);

	double AntennaPattern_Indoor_H(double az_deg);

	double AntennaPattern_Indoor_3D(double az_deg, double el_deg);
}

#endif