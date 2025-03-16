#include "LargeScaleParameters.h"



void LSP::initializeParameters(bool isLos, LinkData& link, VectorXd& Parameters) {
    if (isLos) { LSP::initializeLosParameters(link, Parameters); }
    else { LSP::initializeNlosParameters(link, Parameters); }

}

void LSP::initializeLosParameters(LinkData& link, VectorXd& Parameters) {

    double StandardDeviationSF = Parameters[11];
    double StandardDeviationK = Parameters[13];
    double StandardDeviationDS = Parameters[1];
    double StandardDeviationASD = Parameters[3];
    double StandardDeviationASA = Parameters[5];
    double StandardDeviationZSA = Parameters[7];
    double StandardDeviationZSD = Parameters[9];

    double MeanSF = Parameters[10];
    double MeanK = Parameters[12];
    double MeanDS = Parameters[0];
    double MeanASD = Parameters[2];
    double MeanASA = Parameters[4];
    double MeanZSA = Parameters[6];
    double MeanZSD = Parameters[8];

    double ASDvsDS = Parameters[14]; // +
    double ASAvsDS = Parameters[15]; // +
    double ASAvsSF = Parameters[16]; // +
    double ASDvsSF = Parameters[17]; // +
    double DSvsSF = Parameters[18]; // +
    double ASDvsASA = Parameters[19]; // +
    double ASDvsK = Parameters[20]; // +
    double ASAvsK = Parameters[21]; // +
    double DSvsK = Parameters[22]; // +
    double SFvsK = Parameters[23]; // +
    double ZSDvsSF = Parameters[24]; // +
    double ZSAvsSF = Parameters[25]; // +
    double ZSDvsK = Parameters[26]; // +
    double ZSAvsK = Parameters[27]; // +
    double ZSDvsDS = Parameters[28]; // +
    double ZSAvsDS = Parameters[29]; // +
    double ZSDvsASD = Parameters[30]; //+
    double ZSAvsASD = Parameters[31]; //+
    double ZSDvsASA = Parameters[32]; // +
    double ZSAvsASA = Parameters[33]; //+
    double ZSDvsZSA = Parameters[34]; // +


    Vector<double, 7> value;
#pragma omp parallel for
    for (int i = 0; i < 7; ++i) {
        value(i) = RandomGenerators::generateGauss(0.0, 1.0);
    }
    Vector<double, 7> means;
    Matrix<double, 7, 7> C;


    means << MeanSF, MeanK, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;

    C << StandardDeviationSF * StandardDeviationSF, SFvsK* StandardDeviationSF* StandardDeviationK, DSvsSF* StandardDeviationSF* StandardDeviationDS, ASDvsSF* StandardDeviationSF* StandardDeviationASD, ASAvsSF* StandardDeviationSF* StandardDeviationASA, ZSDvsSF* StandardDeviationSF* StandardDeviationZSD, ZSAvsSF* StandardDeviationSF* StandardDeviationZSA,
        SFvsK* StandardDeviationK* StandardDeviationSF, StandardDeviationK* StandardDeviationK, DSvsK* StandardDeviationK* StandardDeviationDS, ASDvsK* StandardDeviationK* StandardDeviationASD, ASAvsK* StandardDeviationK* StandardDeviationASA, ZSDvsK* StandardDeviationK* StandardDeviationZSD, ZSAvsK* StandardDeviationK* StandardDeviationZSA,
        DSvsSF* StandardDeviationDS* StandardDeviationSF, DSvsK* StandardDeviationDS* StandardDeviationK, StandardDeviationDS* StandardDeviationDS, ASDvsDS* StandardDeviationDS* StandardDeviationASD, ASAvsDS* StandardDeviationDS* StandardDeviationASA, ZSDvsDS* StandardDeviationDS* StandardDeviationZSD, ZSAvsDS* StandardDeviationDS* StandardDeviationZSA,
        ASDvsSF* StandardDeviationASD* StandardDeviationSF, ASDvsK* StandardDeviationASD* StandardDeviationK, ASDvsDS* StandardDeviationASD* StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, ASDvsASA* StandardDeviationASD* StandardDeviationASA, ZSDvsASD* StandardDeviationASD* StandardDeviationZSD, ZSAvsASD* StandardDeviationASD* StandardDeviationZSA,
        ASAvsSF* StandardDeviationASA* StandardDeviationSF, ASAvsK* StandardDeviationASA* StandardDeviationK, ASAvsDS* StandardDeviationASA* StandardDeviationDS, ASDvsASA* StandardDeviationASA* StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, ZSDvsASA* StandardDeviationASA* StandardDeviationZSD, ZSAvsASA* StandardDeviationASA* StandardDeviationZSA,
        ZSDvsSF* StandardDeviationZSD* StandardDeviationSF, ZSDvsK* StandardDeviationZSD* StandardDeviationK, ZSDvsDS* StandardDeviationZSD* StandardDeviationDS, ZSDvsASD* StandardDeviationZSD* StandardDeviationASD, ZSDvsASA* StandardDeviationZSD* StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, ZSDvsZSA* StandardDeviationZSD* StandardDeviationZSA,
        ZSAvsSF* StandardDeviationZSA* StandardDeviationSF, ZSAvsK* StandardDeviationZSA* StandardDeviationK, ZSAvsDS* StandardDeviationZSA* StandardDeviationDS, ZSAvsASD* StandardDeviationZSA* StandardDeviationASD, ZSAvsASA* StandardDeviationZSA* StandardDeviationASA, ZSDvsZSA* StandardDeviationZSA* StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;

    //Matrix<double, 7, 7> L = C.llt().matrixL();

    value = C.llt().matrixL() * value + means;

    link.SF_db = value[0];
    link.K_db = value[1];
    link.DS_sec = pow(10, value[2]);
    link.ASD_deg = pow(10, value[3]);
    link.ASA_deg = pow(10, value[4]);
    link.ZSD_deg = pow(10, value[5]);
    link.ZSA_deg = pow(10, value[6]);



}

void LSP::initializeNlosParameters(LinkData& link, VectorXd& Parameters) {


    double StandardDeviationSF = Parameters[11]; // +
    double StandardDeviationDS = Parameters[1];
    double StandardDeviationASD = Parameters[3];
    double StandardDeviationASA = Parameters[5];
    double StandardDeviationZSA = Parameters[7];
    double StandardDeviationZSD = Parameters[9];

    double MeanSF = Parameters[10];
    double MeanDS = Parameters[0];
    double MeanASD = Parameters[2];
    double MeanASA = Parameters[4];
    double MeanZSA = Parameters[6];
    double MeanZSD = Parameters[8];

    double ASDvsDS = Parameters[12]; // +
    double ASAvsDS = Parameters[13]; // +
    double ASAvsSF = Parameters[14]; // +
    double ASDvsSF = Parameters[15]; // +
    double DSvsSF = Parameters[16]; // +
    double ASDvsASA = Parameters[17]; // +
    double ZSDvsSF = Parameters[18]; // +
    double ZSAvsSF = Parameters[19]; // +
    double ZSDvsDS = Parameters[20]; // +
    double ZSAvsDS = Parameters[21]; // +
    double ZSDvsASD = Parameters[22]; //+
    double ZSAvsASD = Parameters[23]; //+
    double ZSDvsASA = Parameters[24]; // +
    double ZSAvsASA = Parameters[25]; //+
    double ZSDvsZSA = Parameters[26]; // +

    Vector<double, 6> value;
#pragma omp parallel for
    for (int i = 0; i < 6; ++i) {
        value(i) = RandomGenerators::generateGauss(0.0, 1.0);
    }
    Vector<double, 6> means;
    Matrix<double, 6, 6> C;


    means << MeanSF, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;

    C << StandardDeviationSF * StandardDeviationSF, DSvsSF* StandardDeviationSF* StandardDeviationDS, ASDvsSF* StandardDeviationSF* StandardDeviationASD, ASAvsSF* StandardDeviationSF* StandardDeviationASA, ZSDvsSF* StandardDeviationSF* StandardDeviationZSD, ZSAvsSF* StandardDeviationSF* StandardDeviationZSA,
        DSvsSF* StandardDeviationDS* StandardDeviationSF, StandardDeviationDS* StandardDeviationDS, ASDvsDS* StandardDeviationDS* StandardDeviationASD, ASAvsDS* StandardDeviationDS* StandardDeviationASA, ZSDvsDS* StandardDeviationDS* StandardDeviationZSD, ZSAvsDS* StandardDeviationDS* StandardDeviationZSA,
        ASDvsSF* StandardDeviationASD* StandardDeviationSF, ASDvsDS* StandardDeviationASD* StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, ASDvsASA* StandardDeviationASD* StandardDeviationASA, ZSDvsASD* StandardDeviationASD* StandardDeviationZSD, ZSAvsASD* StandardDeviationASD* StandardDeviationZSA,
        ASAvsSF* StandardDeviationASA* StandardDeviationSF, ASAvsDS* StandardDeviationASA* StandardDeviationDS, ASDvsASA* StandardDeviationASA* StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, ZSDvsASA* StandardDeviationASA* StandardDeviationZSD, ZSAvsASA* StandardDeviationASA* StandardDeviationZSA,
        ZSDvsSF* StandardDeviationZSD* StandardDeviationSF, ZSDvsDS* StandardDeviationZSD* StandardDeviationDS, ZSDvsASD* StandardDeviationZSD* StandardDeviationASD, ZSDvsASA* StandardDeviationZSD* StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, ZSDvsZSA* StandardDeviationZSD* StandardDeviationZSA,
        ZSAvsSF* StandardDeviationZSA* StandardDeviationSF, ZSAvsDS* StandardDeviationZSA* StandardDeviationDS, ZSAvsASD* StandardDeviationZSA* StandardDeviationASD, ZSAvsASA* StandardDeviationZSA* StandardDeviationASA, ZSDvsZSA* StandardDeviationZSA* StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;

    //Matrix<double, 6, 6> L = C.llt().matrixL();

    value = C.llt().matrixL() * value + means;

    link.SF_db = value[0];
    link.K_db = -INFINITY;
    link.DS_sec = pow(10, value[1]);
    link.ASD_deg = pow(10, value[2]);
    link.ASA_deg = pow(10, value[3]);
    link.ZSD_deg = pow(10, value[4]);
    link.ZSA_deg = pow(10, value[5]);


}
