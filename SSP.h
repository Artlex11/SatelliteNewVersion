#ifndef SMALL_SCALE_PARAMETERS_H
#define SMALL_SCALE_PARAMETERS_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <Eigen/Dense>

#include "Links.h"
#include "Generators.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace SSP {

	extern double r_tau;
	extern double numberClusters;
	extern double ksi;



	void setParameters(Eigen::VectorXd Parameters);

	void generateClusterDelays(LinkData& link);

	std::pair<std::vector<double>, std::vector<double>>  generateClusterPowers(LinkData& link);

	void calculateLosAngles(LinkData& link, Eigen::Vector3d y_local, Eigen::Vector3d z_local);



	//const std::vector<double> los_C = { 8.0, 5.0, 9.0, 0.375 * std::pow(10, -1.43 * std::log(1 + 30) + 2.228) };
	//const std::vector<double> nlos_C = { 11.0, 5.0, 9.0, 0.375 * std::pow(10, 1.08) };

	//const std::vector<double> am = { 0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492,
	//								  0.3715, -0.3715, 0.5129, -0.5129, 0.6797, -0.6797,
	//								  0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195,
	//								  2.1551, -2.1551 };



	//Eigen::MatrixXd generateAOAandAOD_n_m(bool los, const std::vector<double>& clusterPowers,
	//	double ASAorASD, double riceanK, double AOAorAOD, int AOA_0_or_AOD_1);

	//Eigen::MatrixXd generateZOAandZOD_n_m(bool los, const std::vector<double>& clusterPowers,
	//	double ZSAorZSD, double riceanK, double ZOAorZOD, int ZOA_2_or_ZOD_3);

}

#endif 
