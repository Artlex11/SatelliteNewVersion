#ifndef LINKS_H
#define LINKS_H

#include <Eigen/Dense>
#include <vector>

// Структура для хранения данных антенны ( пока просто как затычка ) 
struct Antenna {
    double gain;               // Коэффициент усиления антенны
};

// Структура для хранения всех необходимых данных
struct LinkData {
    Eigen::Vector3d userPosition;
    Eigen::Vector3d satellitePosition;
    double elevationAngle;
    Antenna userAntenna;
    Antenna satelliteAntenna;
    double PL_db;
    // distance and LOS_p возможно поменять местами
    double distance;
    bool isLOS;
    double K_db;
    double Ds_sec;
    double ASA_deg;
    double ASD_deg;
    double ZSA_deg;
    double ZSD_deg;

};

// Класс для работы с данными связей
class Links {
public:
    // Добавление новой связи (передача по ссылке для избежания копирования)
    void addLink(LinkData& link);

    // Получение всех связей (возврат по ссылке для избежания копирования)
    std::vector<LinkData>& getLinks();

    // Оптимизированная версия: добавление с использованием перемещения (move semantics)
    void addLink(LinkData&& link);


    std::vector<LinkData> links;  // Вектор для хранения всех связей
};

#endif // LINKS_H