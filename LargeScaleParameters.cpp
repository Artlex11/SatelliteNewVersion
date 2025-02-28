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
    link.DS_sec = value[2];
    link.ASD_deg = value[3];
    link.ASA_deg = value[4];
    link.ZSA_deg = value[5];
    link.ZSD_deg = value[6];


}

void LSP::initializeNlosParameters(LinkData& link, VectorXd& Parameters) {


    double StandardDeviationSF = Parameters[11]; // +
    //StandardDeviationK = Parameters[13];
    double StandardDeviationDS = Parameters[1];
    double StandardDeviationASD = Parameters[3];
    double StandardDeviationASA = Parameters[5];
    double StandardDeviationZSA = Parameters[7];
    double StandardDeviationZSD = Parameters[9];

    double MeanSF = Parameters[10];
    //MeanK = Parameters[12];
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
    //ASDvsK = Parameters[20]; // +
    //ASAvsK = Parameters[21]; // +
    //DSvsK = Parameters[22]; // +
    //SFvsK = Parameters[23]; // +
    double ZSDvsSF = Parameters[18]; // +
    double ZSAvsSF = Parameters[19]; // +
    //ZSDvsK = Parameters[26]; // +
    //ZSAvsK = Parameters[27]; // +
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
    link.DS_sec = value[1];
    link.ASD_deg = value[2];
    link.ASA_deg = value[3];
    link.ZSA_deg = value[4];
    link.ZSD_deg = value[5];

    /* double StandardDeviationSF, StandardDeviationDS, StandardDeviationASD, StandardDeviationASA, StandardDeviationZSA, StandardDeviationZSD;
     double MeanSF, MeanDS, MeanASD, MeanASA, MeanZSA, MeanZSD;
     double ASDvsDS, ASAvsDS, ASAvsSF, ASDvsSF, DSvsSF, ASDvsASA, ZSDvsSF, ZSAvsSF, ZSDvsDS,
         ZSAvsDS, ZSDvsASD, ZSAvsASD, ZSDvsASA, ZSAvsASA, ZSDvsZSA;



     Eigen::VectorXd value(6);
     Eigen::VectorXd means(6);
     Eigen::MatrixXd C(6, 6);


     means << MeanSF, MeanDS, MeanASD, MeanASA, MeanZSD, MeanZSA;

     C << StandardDeviationSF * StandardDeviationSF, -0.5 * StandardDeviationSF * StandardDeviationDS, 0.0 * StandardDeviationSF * StandardDeviationASD, -0.4 * StandardDeviationSF * StandardDeviationASA, 0.0 * StandardDeviationSF * StandardDeviationZSD, 0.0 * StandardDeviationSF * StandardDeviationZSA,
         -0.5 * StandardDeviationDS * StandardDeviationSF, StandardDeviationDS* StandardDeviationDS, 0.4 * StandardDeviationDS * StandardDeviationASD, 0.0 * StandardDeviationDS * StandardDeviationASA, -0.27 * StandardDeviationDS * StandardDeviationZSD, -0.06 * StandardDeviationDS * StandardDeviationZSA,
         0.0 * StandardDeviationASD * StandardDeviationSF, 0.4 * StandardDeviationASD * StandardDeviationDS, StandardDeviationASD* StandardDeviationASD, 0.0 * StandardDeviationASD * StandardDeviationASA, 0.35 * StandardDeviationASD * StandardDeviationZSD, 0.23 * StandardDeviationASD * StandardDeviationZSA,
         -0.4 * StandardDeviationASA * StandardDeviationSF, 0.0 * StandardDeviationASA * StandardDeviationDS, 0.0 * StandardDeviationASA * StandardDeviationASD, StandardDeviationASA* StandardDeviationASA, -0.08 * StandardDeviationASA * StandardDeviationZSD, 0.43 * StandardDeviationASA * StandardDeviationZSA,
         0.0 * StandardDeviationZSD * StandardDeviationSF, -0.27 * StandardDeviationZSD * StandardDeviationDS, 0.35 * StandardDeviationZSD * StandardDeviationASD, -0.08 * StandardDeviationZSD * StandardDeviationASA, StandardDeviationZSD* StandardDeviationZSD, 0.42 * StandardDeviationZSD * StandardDeviationZSA,
         0.0 * StandardDeviationZSA * StandardDeviationSF, -0.06 * StandardDeviationZSA * StandardDeviationDS, 0.23 * StandardDeviationZSA * StandardDeviationASD, 0.43 * StandardDeviationZSA * StandardDeviationASA, 0.42 * StandardDeviationZSA * StandardDeviationZSD, StandardDeviationZSA* StandardDeviationZSA;

     Eigen::MatrixXd L = C.llt().matrixL();

     for (int i = 0; i < 6; ++i) {
         value(i) = RandomGenerators::generateGauss(0.0, 1.0);
     }

     Eigen::VectorXd value_new = L * value + means;

     shadowFading = value_new(0);
     delaySpread = value_new(1);
     azimuthSpreadDeparture = std::min(value_new(2), log10(104.0));
     azimuthSpreadArrival = std::min(value_new(3), log10(104.0));
     zenithSpreadDeparture = std::min(value_new(4), log10(52.0));
     zenithSpreadArrival = std::min(value_new(5), log10(52.0));
 }

 void LargeScaleParameters::showParameters() const {
     if (los) {
         std::cout << "SF [dB] : " << shadowFading << ",\nK [dB] : " << riceanK << ",\nDS [log10(DS/1s)] : " << delaySpread << ",\nASD [log10(ASD/ 1* degree] : " << azimuthSpreadDeparture << ",\nASA [log10(ASA/ 1* degree] : " << azimuthSpreadArrival << ",\nZSD [log10(ZSD/ 1* degree] : " << zenithSpreadDeparture << ",\nZSA [log10(ZSA/ 1* degree] : " << zenithSpreadArrival << std::endl << std::endl;
     }
     else {
         std::cout << "SF [dB] : " << shadowFading << ",\nDS [log10(DS/1s)] : " << delaySpread << ",\nASD [log10(ASD/ 1* degree] : " << azimuthSpreadDeparture << ",\nASA [log10(ASA/ 1* degree] : " << azimuthSpreadArrival << ",\nZSD [log10(ZSD/ 1* degree] : " << zenithSpreadDeparture << ",\nZSA [log10(ZSA/ 1* degree] : " << zenithSpreadArrival << std::endl << std::endl;
     }*/
}
