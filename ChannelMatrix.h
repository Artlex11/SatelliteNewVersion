#ifndef CHANNEL_MATRIX_H
#define CHANNEL_MATRIX_H

#include "Eigen/Dense"
#include <iostream>
#include "Generators.h"
#include "links.h"
#include <random>
#include "Antennas.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif




namespace ChannelMatrixH {

	Eigen::MatrixXcd generateRayGain(LinkData& link);
}


#endif // !CHANNEL_MATRIX_H

