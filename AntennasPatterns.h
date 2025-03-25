#ifndef ANTENNAS_H
#define ANTENNAS_H
#include <iostream>
#include <Eigen\Dense>
#include "links.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace AntennaPattern {
	//Sat
	double AntennaPattern_NTN1_3D(double az_deg, double el_deg);
	//UE for S-band
	double AntennaPattern_Omni();
	//mb UE
	double AntennaPattern_38901_3D(double az_deg, double el_deg);
	//UE for 30 and 63 GHz
	double AntennaPattern_37885_PedCellUE_3D(double az_deg, double el_deg);
}

namespace AntennaArray {

	//Mg - numbers panels in row, Ng - number panels col, P - numbers polarizations, M - numbers antenna element in row (panel), N - numbers  antenna element in col (panel), Frequency in GHz
	Eigen::MatrixXd createAntennaArray_Mg_Ng_P_M_N(int Mg, int Ng, int M, int N, double frequency, double spacing);

	Eigen::MatrixXd rotateAntennaArray(double alpha, double beta, const Eigen::MatrixXd& antennaArray);

	void addAntennaUE(LinkData& link, int numElemensInCol, int numElemensInRow, std::string pattern, double frequency, std::vector<double> anglesPol, double spacing);
	void addAntennaSat(LinkData& link, int numElemensInCol, int numElemensInRow, std::string pattern, double frequency, std::vector<double> anglesPol, double spacing);

	Eigen::Vector2d FieldPatternWithPolarization(double& antennaPattern, double anglePol);


	//Eigen::Vector2d getFieldPatternUE(LinkData& link);
	//Eigen::Vector2d getFieldPatternSat(LinkData& link);
}

#endif // !ANTENNAS_H


