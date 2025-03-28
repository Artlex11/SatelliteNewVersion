#ifndef LINKS_H
#define LINKS_H

#include <Eigen/Dense>
#include <vector>

struct Antenna 
{
    Eigen::MatrixXd antennaArray;
    std::string pattern;
    std::vector<double> anglesPol;

    double  bearingAntArray;
    double  elevationAntArray;
};

// ��������� ��� �������� ���� ����������� ������
struct LinkData 
{
    Eigen::Vector3d userPosition;
    Eigen::Vector3d satellitePosition;
    double elevationAngle;

    bool isLos;
    double PL_db;

    //Large scale parameters:
    double SF_db;
    double K_db;
    double DS_sec;
    double ASA_deg; // ASA � ��������
    double ASD_deg; // ASD � ��������
    double ZSA_deg; // ZSA � ��������
    double ZSD_deg; // ZSD � ��������

    //Small scale parameters:
    double AoD_Los; // AoD � ��������
    double AoA_Los; // AoA � ��������
    double ZoD_Los; // ZoD � ��������
    double ZoA_Los; // ZoA � ��������
    std::vector<double> clusterDelays;
    std::vector<double> clusterScaledDelays; // ������ ��� ������ ��������� , �� �� ������������ � ����������� 

    Eigen::MatrixXd powerInRays;
    Eigen::MatrixXd delayInRays;



    Eigen::MatrixXd AOD_n_m;
    Eigen::MatrixXd AOA_n_m;
    Eigen::MatrixXd ZOD_n_m;
    Eigen::MatrixXd ZOA_n_m;
    Eigen::MatrixXd XPR_n_m;

    Antenna antennaUE;
    Antenna antennaSat;

};

// ����� ��� ������ � ������� ������
class Links 
{

public:

    std::vector<LinkData>& getLinks();

    // ���������������� ������: ���������� � �������������� ����������� (move semantics)
    void addLink(LinkData&& link);

private:
    std::vector<LinkData> links;  // ������ ��� �������� ���� ������
};

#endif // LINKS_H