#include "LOS_Probability.h"
#include "Tables.h"

#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <Eigen/Dense>
#include <iomanip> // Для std::setprecision

// Глобальные константы
const double EARTH_RADIUS = 6371000.0; // Радиус Земли в метрах

const double PI = 3.14159265358979323846; // Число Пи

// Функция поворота Родрига
Eigen::MatrixXd rotateAroundAxis(const Eigen::MatrixXd& inVect, const Eigen::Vector3d& ax, double angl)
{
    // Check if the input vector is 3-by-N
    if (inVect.rows() != 3) {
        throw std::invalid_argument("Input vector should be 3-by-N.");
    }

    // Normalize the rotation axis
    Eigen::Vector3d normalizedAx = ax.normalized();

    // Calculate the components of the rotation axis
    double x = normalizedAx(0);
    double y = normalizedAx(1);
    double z = normalizedAx(2);
    double c = cos(angl);
    double s = sin(angl);

    // Construct the rotation matrix
    Eigen::Matrix3d rotMatrix;
    rotMatrix << c + (1 - c) * x * x, (1 - c)* x* y - s * z, (1 - c)* x* z + s * y,
        (1 - c)* y* x + s * z, c + (1 - c) * y * y, (1 - c)* y* z - s * x,
        (1 - c)* z* x - s * y, (1 - c)* z* y + s * x, c + (1 - c) * z * z;

    // Perform the rotation
    Eigen::MatrixXd outVect = rotMatrix * inVect;

    return outVect;
}

class Satellite
{

private:

    double x, y, z; // Координаты спутника
    const double beamWidth; // Ширина луча
    double elMinDegrees = 10.0; // Минимальный угол ( места ) видимости 
    double elTargetDegrees = 90.0; // Угол цели спутника 
    Eigen::MatrixXd uvSet; // Матрица для хранения UV-координат

public:

    // Конструктор класс и инициализация параметров
    Satellite(double x, double y, double z, double beamWidth)
        : x(x), y(y), z(z + EARTH_RADIUS), beamWidth(beamWidth)
    {

    }

    double getBeamWidth() const
    {
        return beamWidth;
    }

    void printCoordinates() const
    {
        std::cout << "Satellite Coordinates: (" << x << ", " << y << ", " << z << ")\n";
    }

    // метод для расчета угла места
    void calculateZenithAngle(Eigen::MatrixXd& users)
    {
        // Проверка, что матрица имеет 3 строки
        if (users.rows() != 3)
        {
            std::cerr << "Error: The input matrix must have 3 rows for x, y, z coordinates." << std::endl;
            return;
        }

        // Проверка, что матрица имеет хотя бы одну колонку
        int nUEs = users.cols();
        if (nUEs == 0)
        {
            std::cerr << "Error: The input matrix must have at least one column." << std::endl;
            return;
        }

        // Добавляем четвертую строку для углов
        users.conservativeResize(4, nUEs);

        for (int col = 0; col < nUEs; ++col)
        {
            // Разностный вектор
            double dx = x - users(0, col);
            double dy = y - users(1, col);
            double dz = z - users(2, col);

            // Нормируем разностный вектор
            double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
            double norm_dx = dx / distance;
            double norm_dy = dy / distance;
            double norm_dz = dz / distance;

            // Нормальный вектор пользователя (от центра Земли)
            double normal_x = users(0, col) / EARTH_RADIUS;
            double normal_y = users(1, col) / EARTH_RADIUS;
            double normal_z = users(2, col) / EARTH_RADIUS;

            // Скалярное произведение
            double dot_product = normal_x * norm_dx + normal_y * norm_dy + normal_z * norm_dz;
            std::cout << dot_product << std::endl;

            // Находим угол места (аркосинус)
            double elevation_angle = std::acos(dot_product) * 180 / PI; // почти все углы в интервале от 90 до 180 - Земля прозрачная, надо отсечь

            // Устанавливаем угол в четвертую строку
            users(3, col) = elevation_angle;
        }
    }

    // разброс пользователей
    Eigen::MatrixXd dropUEs(int nTiers, int nUePerCell) {
        double beamWidthRadians = beamWidth * (PI / 180.0);
        double uvStep = sin(beamWidthRadians * 0.865); // Шаг UV
        double r = nTiers * uvStep;

        // Генерация координат
        int size = (2 * nTiers + 1);
        Eigen::MatrixXd uvSet(size * size, 2); // 2 столбца для U и V
        int index = 0;

        for (int i = -nTiers; i <= nTiers; ++i) {
            for (int j = -nTiers; j <= nTiers; ++j) {
                // Расчет U и V координат
                double u_val = i * uvStep - (j * (uvStep / 2));
                double v_val = (j * (sqrt(3) * uvStep / 2));

                // Проверка, находится ли точка внутри шестиугольника
                if (fabs(u_val + v_val / sqrt(3)) <= r + 1e-10 &&
                    fabs(u_val - v_val / sqrt(3)) <= r + 1e-10) {
                    uvSet(index, 0) = u_val;
                    uvSet(index, 1) = v_val;
                    index++;
                }
            }
        } // для чего в проверке поправка 1e-10

        // Обрезаем матрицу до фактического размера
        uvSet.conservativeResize(index, Eigen::NoChange);

        std::cout << "uvStep: " << uvStep << std::endl;
        std::cout << "uvSet rows: " << uvSet.rows() << std::endl;
        //std::cout << "uvSet:\n" << uvSet << std::endl;

        // Три угла из Matlaba
        double phiMax = PI - (PI / 180.0 * (90.0 + elMinDegrees)) - std::asin(std::sin(PI / 180.0 * (90.0 + elMinDegrees)) * EARTH_RADIUS / (EARTH_RADIUS + z));
        double phiSat = PI - (PI / 180.0 * (90.0 + elTargetDegrees)) - std::asin(std::sin(PI / 180.0 * (90.0 + elTargetDegrees)) * EARTH_RADIUS / (EARTH_RADIUS + z));
        double thetaSat = PI - phiSat - (PI / 180.0 * (90 + elTargetDegrees));

        // Ещё один угол
        double threshold = EARTH_RADIUS * std::cos(phiMax);

        int nUEs = nUePerCell * index; // Количество пользователей на основе сгенерированных UV координат
        Eigen::MatrixXd xyzUEs_all(3, nUEs * 10); // Вектор для хранения координат пользователей
        Eigen::MatrixXd finalxyzUEs_all(3, nUEs);

        int UEcnt = 0;

        while (UEcnt != nUEs)
        {

            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<double> dist(0.0, 1.0);

            for (int i = 0; i < nUEs * 10; ++i)
            {
                double xyz[3] = { dist(gen), dist(gen), dist(gen) }; // Трёхмерный Гаусс
                double VectorLenght = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]); // радиус вектор до сгенерированного пользователя
                xyzUEs_all(0, i) = EARTH_RADIUS * (xyz[0] / VectorLenght);
                xyzUEs_all(1, i) = EARTH_RADIUS * (xyz[1] / VectorLenght);
                xyzUEs_all(2, i) = EARTH_RADIUS * (xyz[2] / VectorLenght);
            }

            // Фильтрация столбцов
            std::vector<int> valid_indices;
            for (int i = 0; i < nUEs; ++i)
            {
                if (abs(xyzUEs_all(2, i)) >= threshold)
                {
                    valid_indices.push_back(i);
                }
            }
            // я тут
            // 
            // 
            // 
            // 
            // Создание новой матрицы с отфильтрованными столбцами
            Eigen::MatrixXd filtered_xyzUEs_all(3, valid_indices.size());
            for (size_t i = 0; i < valid_indices.size(); ++i)
            {
                filtered_xyzUEs_all.col(i) = xyzUEs_all.col(valid_indices[i]);
            }

            Eigen::Vector3d xyzSat(x, y, z);

            Eigen::MatrixXd xyzUEstoSat(3, filtered_xyzUEs_all.cols());
            Eigen::MatrixXd uvUEstoSat(2, filtered_xyzUEs_all.cols());
            Eigen::MatrixXd uvDistance(uvSet.rows(), uvUEstoSat.cols());
            // Вычисляем разности
            for (int i = 0; i < xyzUEstoSat.cols(); ++i)
            {
                xyzUEstoSat.col(i) = xyzSat - xyzUEs_all.col(i);
                xyzUEstoSat.col(i).normalize();
                xyzUEstoSat.col(i) = rotateAroundAxis(xyzUEstoSat.col(i), Eigen::Vector3d(0, 0, 1), thetaSat);
                uvUEstoSat(0, i) = xyzUEstoSat(1, i); // y-координата
                uvUEstoSat(1, i) = xyzUEstoSat(2, i); // z-координата
            }

            // Вычисляем расстояния
            for (int i = 0; i < uvSet.rows(); ++i) {
                for (int j = 0; j < uvUEstoSat.cols(); ++j)
                {
                    // Вычисляем евклидово расстояние между uvSet и uvUEstoSat
                    double deltaU = uvUEstoSat(0, j) - uvSet(i, 0);
                    double deltaV = uvUEstoSat(1, j) - uvSet(i, 1);
                    uvDistance(i, j) = std::sqrt(deltaU * deltaU + deltaV * deltaV); // Расстояние
                }
            }

            // Применяем условие abs(uvDistance) <= uvBeamRadius
            Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> indCellUE(1, uvUEstoSat.cols());
            indCellUE.setZero(); // Инициализация
            double uvBeamRadius = uvStep * sqrt(3);
            for (size_t j = 0; j < uvUEstoSat.cols(); ++j)
            {
                indCellUE(0, j) = false; // Изначально считаем, что ничего не найдено
                for (size_t i = 0; i < uvSet.rows(); ++i) {
                    if (uvDistance(i, j) <= uvBeamRadius) {
                        indCellUE(0, j) = true; // Нашли хотя бы одно значение удовлетворяющее условию
                        break; // Можно прекратить поиск
                    }
                }
            }

            // Извлечение хороших UE
            Eigen::MatrixXd xyzGoodUEs; // Для хранения допустимых UE
            std::vector<int> goodUEIndices;

            for (int j = 0; j < uvUEstoSat.cols(); ++j) {
                if (indCellUE(0, j)) {
                    goodUEIndices.push_back(j);
                }
            }

            // Получаем количество хороших UE
            int nGoodUEs = goodUEIndices.size();
            int nUEsToInclude = std::min(nUEs - UEcnt, nGoodUEs);

            // Копируем координаты хороших UE в xyzUEs_all
            for (int k = 0; k < nUEsToInclude; ++k) {
                finalxyzUEs_all.col(k + UEcnt) = filtered_xyzUEs_all.col(goodUEIndices[k]);
            }

            UEcnt += nUEsToInclude;
        }
        // После завершения заполнения матрицы xyzUEs_all
        /*std::cout << "UE Coordinates (x, y, z):" << std::endl;
        for (int col = 0; col < finalxyzUEs_all.cols(); ++col) {
            std::cout << "Column " << col + 1 << ": ("
                << finalxyzUEs_all(0, col) << ", "
                << finalxyzUEs_all(1, col) << ", "
                << finalxyzUEs_all(2, col) << ")" << std::endl;
        }*/
        return finalxyzUEs_all;
    }


};

int main() {
    
    // сферические координаты для задания случайного положения спутника

    double altitude = 1000000.0;
    double r = EARTH_RADIUS + altitude;
    double theta = RandomGenerators::generateUniform(0.0, PI);
    double phi = RandomGenerators::generateUniform(0.0, 2 * PI);

    double x = r * cos(phi) * sin(theta);
    double y = r * sin(phi) * sin(theta);
    double z = r * cos(theta);

    std::cout << sqrt(x * x + y * y + z * z) << std::endl;

    // 2 вариант: 
    /*Можно сгенерировать случайные координаты на единичной сфере(с помощью нормализации) и затем масштабировать их на нужный радиус.*/

    Satellite satellite(x, y, z, 4.4127); // Спутник над Землёй
    satellite.printCoordinates();

    int nTiers; // Переменная для хранения количества уровней
    std::cout << "Input number of levels (nTiers): ";
    std::cin >> nTiers; // Ввод количества уровней

    Eigen::MatrixXd users = satellite.dropUEs(nTiers, 10); // Генерируем UV-плоскость
    satellite.calculateZenithAngle(users);
    std::cout << "UE Coordinates (x, y, z):" << std::endl;
    for (int col = 0; col < users.cols(); ++col) {
        std::cout << "Column " << col + 1 << ": ("
            << users(0, col) << ", "
            << users(1, col) << ", "
            << users(2, col) << ", "
            << users(3, col) << " grad )" << std::endl;
    }

    std::string scenario;
    double f;
    bool los;

    std::cout << "\nInput scenario (Dense_Urban/Urban/Suburban/Rural): ";
    std::cin >> scenario;
    std::cout << "\nInput frequency(f) in GHz: ";
    std::cin >> f;

    for (int col = 0; col < users.cols(); ++col)
    {
        double deg = users(3, col);

        int index = int(AngleForLSP(deg)) / 10;
        std::cout << "User " << col << ": index: " << index << ", angle: " << AngleForLSP(deg) << "\n";
        los = CalculateLOSProbability(index, scenario);

        MatrixXd Table = GenerateMatrix(los, f, scenario);

        if (los)
        {
            VectorXd Parameters{ 45 };
            Parameters = Table.col(index - 1);
        }
        else
        {
            VectorXd Parameters{ 38 };
            Parameters = Table.col(index - 1);
        }

    }

    return 0;
}
