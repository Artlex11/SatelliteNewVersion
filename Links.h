#ifndef LINKS_H
#define LINKS_H

#include <Eigen/Dense>
#include <vector>


// Структура для хранения всех необходимых данных
struct LinkData {
    Eigen::Vector3d userPosition;
    Eigen::Vector3d satellitePosition;
    double elevationAngle;

    bool isLos;
    double PL_db;

    //Large scale parameters:
    double SF_db;
    double K_db;
    double DS_sec;
    double ASA_deg; // ASA в градусах
    double ASD_deg; // ASD в градусах
    double ZSA_deg; // ZSA в градусах
    double ZSD_deg; // ZSD в градусах

    //Small scale parameters:
    double AoD_Los; // AoD в градусах
    double AoA_Los; // AoA в градусах
    double ZoD_Los; // ZoD в градусах
    double ZoA_Los; // ZoA в градусах
    std::vector<double> clusterDelays;
    std::vector<double> clusterScaledDelays; // только для прямой видимости , но не используется в вычислениях 



};

// Класс для работы с данными связей
class Links {

public:

    std::vector<LinkData>& getLinks();

    // Оптимизированная версия: добавление с использованием перемещения (move semantics)
    void addLink(LinkData&& link);

private:
    std::vector<LinkData> links;  // Вектор для хранения всех связей
};

#endif // LINKS_H