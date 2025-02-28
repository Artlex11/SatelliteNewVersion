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

    extern double StandardDeviationSF, StandardDeviationK, StandardDeviationDS, StandardDeviationASD, StandardDeviationASA, StandardDeviationZSA, StandardDeviationZSD;

    extern double MeanSF, MeanK, MeanDS, MeanASD, MeanASA, MeanZSA, MeanZSD;

    extern double ASDvsDS, ASAvsDS, ASAvsSF, ASDvsSF, DSvsSF, ASDvsASA, ASDvsK, ASAvsK, DSvsK, SFvsK, ZSDvsSF, ZSAvsSF, ZSDvsK, ZSAvsK, ZSDvsDS,
        ZSAvsDS, ZSDvsASD, ZSAvsASD, ZSDvsASA, ZSAvsASA, ZSDvsZSA;

    void initializeParameters(bool isLos, LinkData& link, VectorXd& Parameters);
    void initializeLosParameters(LinkData& link, VectorXd& Parameters);
    void initializeNlosParameters(LinkData& link, VectorXd& Parameters);

};

#endif