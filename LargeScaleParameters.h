#ifndef LARGE_SCALE_PARAMETERS_H
#define LARGE_SCALE_PARAMETERS_H

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <omp.h> 

#include "links.h"
#include "Generators.h"

using namespace Eigen;

namespace LSP {

    void initializeParameters(bool isLos, LinkData& link, VectorXd& Parameters);
    void initializeLosParameters(LinkData& link, VectorXd& Parameters);
    void initializeNlosParameters(LinkData& link, VectorXd& Parameters);

};

#endif
