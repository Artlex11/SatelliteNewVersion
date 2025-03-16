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


namespace SSP {


	//3GPP TR 38.811 V15.4.0 (2020-09) Table 6.7.2-1aa
	const std::vector<double> tableScalingFactorsAOAandAOD = { 0.0, 0.501 , 0.680, 0.779, 0.860, 0.0, 0.0, 1.090, 1.123, 1.146, 0.0, 1.190, 1.211, 1.226, 0.0, 0.0, 1.273, 1.289 };
	//3GPP TR 38.811 V15.4.0 (2020-09) Table 6.7.2-1ab
	const std::vector<double> tableScalingFactorsZOAandZOD = { 0.0, 0.430, 0.594, 0.697, 0.0, 0.0, 0.0, 0.889, 0.0, 0.957, 1.031, 1.104, 0.0, 0.0, 1.1088, 0.0, 0.0, 0.0, 1.184, 1.178 };

	const std::vector<double> rayOfsetAngles = { 0.0447, -0.0447, 0.1413, -0.1413, 0.2492, -0.2492,
												0.3715, -0.3715, 0.5129, -0.5129, 0.6797, -0.6797,
												0.8844, -0.8844, 1.1481, -1.1481, 1.5195, -1.5195,
												2.1551, -2.1551 };







	void setParameters(Eigen::VectorXd Parameters);
	void generateClusterDelays(LinkData& link);
	std::pair<std::vector<double>, std::vector<double>>  generateClusterPowers(LinkData& link);



	void calculateLosAngles(LinkData& link, Eigen::Vector3d y_local, Eigen::Vector3d z_local);
	void generateArrivalAndDepartureAngles(LinkData& link, std::vector<double> clusterPowersWithScalingFactors);
	void generateXPR(LinkData& link);



	//std::vector<double> calculateAngularSpreadandMeanAngles(bool los, const std::vector<double>& clusterPowers,
		//Eigen::MatrixXd& AOD, Eigen::MatrixXd& AOA,
		//Eigen::MatrixXd& ZOD, Eigen::MatrixXd& ZOA);

	// Функция для перемешивания подгруппы в строке (без копирования данных)




	const std::vector<size_t> subCluster1 = { 0, 1, 2, 3, 4, 5, 6, 7, 18, 19 };
	const std::vector<size_t> subCluster2 = { 8, 9, 10, 11, 16, 17 };
	const std::vector<size_t> subCluster3 = { 12, 13, 14, 15 };
	const std::vector<size_t> cluster = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 , 16, 17, 18, 19 };

	void mixingClusterInRow(Eigen::VectorXd& row, const std::vector<size_t>& indices, std::minstd_rand& g);

	// Функция для перемешивания подгрупп во всех строках матрицы
	void mixingMatrix(Eigen::MatrixXd& matrix, std::minstd_rand& g);

	// Основной метод для перемешивания подгрупп в структуре LinkData
	void mixingAngles(LinkData& data, std::minstd_rand& g);


	void sortRelativeToFirstVector(
		std::vector<double>& mainVector,
		std::vector<double>& secondVector,
		Eigen::MatrixXd& matrix1,
		Eigen::MatrixXd& matrix2,
		Eigen::MatrixXd& matrix3,
		Eigen::MatrixXd& matrix4,
		Eigen::MatrixXd& matrix5
	);





}
#endif 
