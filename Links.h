#ifndef LINKS_H
#define LINKS_H

#include <Eigen/Dense>
#include <vector>

// ��������� ��� �������� ������ ������� ( ���� ������ ��� ������� ) 
struct Antenna {
    double gain;               // ����������� �������� �������
};

// ��������� ��� �������� ���� ����������� ������
struct LinkData {
    Eigen::Vector3d userPosition;
    Eigen::Vector3d satellitePosition;
    double elevationAngle;
    Antenna userAntenna;
    Antenna satelliteAntenna;
    double PL_db;
    // distance and LOS_p �������� �������� �������
    double distance;
    bool isLOS;
    double K_db;
    double Ds_sec;
    double ASA_deg;
    double ASD_deg;
    double ZSA_deg;
    double ZSD_deg;

};

// ����� ��� ������ � ������� ������
class Links {
public:
    // ���������� ����� ����� (�������� �� ������ ��� ��������� �����������)
    void addLink(LinkData& link);

    // ��������� ���� ������ (������� �� ������ ��� ��������� �����������)
    std::vector<LinkData>& getLinks();

    // ���������������� ������: ���������� � �������������� ����������� (move semantics)
    void addLink(LinkData&& link);


    std::vector<LinkData> links;  // ������ ��� �������� ���� ������
};

#endif // LINKS_H