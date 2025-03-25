#ifndef SMALL_SCALE_PARAMETERS_H
#define SMALL_SCALE_PARAMETERS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <numeric> 

#include <random>
#include <algorithm> // для std::shuffle
#include <omp.h> // Подключаем OpenMP

#include "links.h"
#include "Generators.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//3GPP TR 38.811 V15.4.0 (2020-09) Table 6.7.2-1aa
const std::vector<double> tableScalingFactorsAOAandAOD = { 0.0, 0.501 , 0.680, 0.779, 0.860, 0.0, 0.0, 1.018, 0.0, 1.090, 1.123, 1.146, 0.0, 1.190, 1.211, 1.226, 0.0, 0.0, 1.273, 1.289 };
//3GPP TR 38.811 V15.4.0 (2020-09) Table 6.7.2-1ab
const std::vector<double> tableScalingFactorsZOAandZOD = { 0.0, 0.430, 0.594, 0.697, 0.0, 0.0, 0.0, 0.889, 0.0, 0.957, 1.031, 1.104, 0.0, 0.0, 1.1088, 0.0, 0.0, 0.0, 1.184, 1.178 };

const std::vector<double> rayOfsetAngles = { 0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492,
											0.3715, -0.3715, 0.5129, -0.5129, 0.6797, -0.6797,
											0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195,
											2.1551, -2.1551 };

const std::vector<int> subCluster1 = { 0, 1, 2, 3, 4, 5, 6, 7, 18, 19 };
const std::vector<int> subCluster2 = { 8, 9, 10, 11, 16, 17 };
const std::vector<int> subCluster3 = { 12, 13, 14, 15 };
const std::vector<int> cluster = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 , 16, 17, 18, 19 };


namespace SSP {

	void setParameters(Eigen::VectorXd Parameters);

	void generateClusterDelays(LinkData& link);

	std::pair<std::vector<double>, std::vector<double>>  generateClusterPowers(LinkData& link);



	void calculateLosAngles(LinkData& link);
	void generateArrivalAndDepartureAngles(LinkData& link, std::vector<double> clusterPowersWithScalingFactors);

	void generateXPR(LinkData& link);

	void mixingClusterInRow(Eigen::VectorXd& row, const std::vector<int>& indices, std::minstd_rand& g);
	void mixingMatrix(Eigen::MatrixXd& matrix, std::minstd_rand& g);
	void mixingAngles(LinkData& link, std::minstd_rand& g);


	void sortRelativeToFirstVector(LinkData& link, std::vector<double>& clusterPowers);

	void transformVectors2MatForLink(LinkData& link, std::vector<double>& vec);


}
#endif 
