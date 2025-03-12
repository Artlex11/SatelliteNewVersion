#include "PathLoss.h"

// S-band
Vector<double, 9> SF_std_Dense_Urban_LOS_S = { 3.5, 3.4, 2.9, 3.0, 3.1, 2.7, 2.5, 2.3, 1.2 };
Vector<double, 9> SF_std_Urban_LOS_S = { 4, 4, 4, 4, 4, 4, 4, 4, 4 };
Vector<double, 9> SF_std_Suburban_LOS_S = { 1.79, 1.14, 1.14, 0.92, 1.42, 1.56, 0.85, 0.72, 0.72 };
Vector<double, 9> SF_std_Rural_LOS_S = { 1.79, 1.14, 1.14, 0.92, 1.42, 1.56, 0.85, 0.72, 0.72 };

// Ka-band
Vector<double, 9> SF_std_Dense_Urban_LOS_Ka = { 2.9, 2.4, 2.7, 2.4, 2.4, 2.7, 2.6, 2.8, 0.6 };
Vector<double, 9> SF_std_Urban_LOS_Ka = { 4, 4, 4, 4, 4, 4, 4, 4, 4 };
Vector<double, 9> SF_std_Suburban_LOS_Ka = { 1.9, 1.6, 1.9, 2.3, 2.7, 3.1, 3.0, 3.6, 0.4 };
Vector<double, 9> SF_std_Rural_LOS_Ka = { 1.9, 1.6, 1.9, 2.3, 2.7, 3.1, 3.0, 3.6, 0.4 };

// NLOS

// S-band
Vector<double, 9> SF_std_Dense_Urban_NLOS_S = { 15.5, 13.9, 12.4, 11.7, 10.6, 10.5, 10.1, 9.2, 9.2 };
Vector<double, 9> SF_std_Urban_NLOS_S = { 6, 6, 6, 6, 6, 6, 6, 6, 6 };
Vector<double, 9> SF_std_Suburban_NLOS_S = { 8.93, 9.08, 8.78, 10.25, 10.56, 10.74, 10.17, 11.52, 11.52 };
Vector<double, 9> SF_std_Rural_NLOS_S = { 8.93, 9.08, 8.78, 10.25, 10.56, 10.74, 10.17, 11.52, 11.52 };

// Ka-band
Vector<double, 9> SF_std_Dense_Urban_NLOS_Ka = { 17.1, 17.1, 15.6, 14.6, 14.2, 12.6, 12.1, 12.3, 12.3 };
Vector<double, 9> SF_std_Urban_NLOS_Ka = { 6, 6, 6, 6, 6, 6, 6, 6, 6 };
Vector<double, 9> SF_std_Suburban_NLOS_Ka = { 10.7, 10.0, 11.2, 11.6, 11.8, 10.8, 10.8, 10.8, 10.8 };
Vector<double, 9> SF_std_Rural_NLOS_Ka = { 10.7, 10.0, 11.2, 11.6, 11.8, 10.8, 10.8, 10.8, 10.8 };


// ClusterLoss (for NLOS case)

// S-band
Vector<double, 9> CL_Dense_Urban_S = { 34.3, 30.9, 29.0, 27.7, 26.8, 26.2, 25.8, 25.5, 25.5 };
Vector<double, 9> CL_Urban_S = { 34.3, 30.9, 29.0, 27.7, 26.8, 26.2, 25.8, 25.5, 25.5 };
Vector<double, 9> CL_Suburban_S = { 19.52, 18.17, 18.42, 18.28, 18.63, 17.68, 16.50, 16.30, 16.30 };
Vector<double, 9> CL_Rural_S = { 19.52, 18.17, 18.42, 18.28, 18.63, 17.68, 16.50, 16.30, 16.30 };

// Ka-band
Vector<double, 9> CL_Dense_Urban_Ka = { 44.3, 39.9, 37.5, 35.8, 34.6, 33.8, 33.3, 33.0, 32.9 };
Vector<double, 9> CL_Urban_Ka = { 44.3, 39.9, 37.5, 35.8, 34.6, 33.8, 33.3, 33.0, 32.9 };
Vector<double, 9> CL_Suburban_Ka = { 29.5, 24.6, 21.9, 20.0, 18.7, 17.8, 17.2, 16.9, 16.8 };
Vector<double, 9> CL_Rural_Ka = { 29.5, 24.6, 21.9, 20.0, 18.7, 17.8, 17.2, 16.9, 16.8 };

// Расчёт базовых потерь
double CalculateDistance(const double EARTH_RADIUS, double altitude, double alpha)
{
	double d = sqrt(EARTH_RADIUS * EARTH_RADIUS * sin(alpha) * sin(alpha) + altitude * altitude + 2 * altitude * EARTH_RADIUS) - EARTH_RADIUS * sin(alpha);
	return d;
}

double ChooseSTD(bool los, double f, double alpha, std::string scenario)
{
	int deg = int(alpha / 10);
	double std;
	if (los)
	{
		if (f < 6.0)
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_LOS_S[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_LOS_S[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_LOS_S[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_LOS_S[deg];
			}
		}
		else
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_LOS_Ka[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_LOS_Ka[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_LOS_Ka[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_LOS_Ka[deg];
			}
		}
	}
	else
	{
		if (f < 6.0)
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_NLOS_S[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_NLOS_S[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_NLOS_S[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_NLOS_S[deg];
			}
		}
		else
		{
			if (scenario == "DenseUrban")
			{
				std = SF_std_Dense_Urban_NLOS_Ka[deg];
			}
			else if (scenario == "Urban")
			{
				std = SF_std_Urban_NLOS_Ka[deg];
			}
			else if (scenario == "Suburban")
			{
				std = SF_std_Suburban_NLOS_Ka[deg];
			}
			else if (scenario == "Rural")
			{
				std = SF_std_Rural_NLOS_Ka[deg];
			}
		}
	}
	return std;
}

double GenerateSF(double std)
{
	double SF = RandomGenerators::generateGauss(0, std);
	return SF;
}

double Calculate_FSPL(double d, double f)
{
	double FSPL = 32.45 + 20 * log10(f) + 20 * log10(d);
	return FSPL;
}

double ChooseCL(bool los, double f, double alpha, std::string scenario)
{
	int deg = int(alpha / 10) - 1;
	double CL;
	if (los)
	{
		CL = 0;
	}
	else
	{
		if (f < 6.0)
		{
			if (scenario == "DenseUrban")
			{
				CL = CL_Dense_Urban_S[deg];
			}
			else if (scenario == "Urban")
			{
				CL = CL_Urban_S[deg];
			}
			else if (scenario == "Suburban")
			{
				CL = CL_Suburban_S[deg];
			}
			else if (scenario == "Rural")
			{
				CL = CL_Rural_S[deg];
			}
		}
		else
		{
			if (scenario == "DenseUrban")
			{
				CL = CL_Dense_Urban_Ka[deg];
			}
			else if (scenario == "Urban")
			{
				CL = CL_Urban_Ka[deg];
			}
			else if (scenario == "Suburban")
			{
				CL = CL_Suburban_Ka[deg];
			}
			else if (scenario == "Rural")
			{
				CL = CL_Rural_Ka[deg];
			}
		}
	}
	return CL;
}

double CalculateBasisPathLoss(double FSPL, double SF, double CL)
{
	double PL_b = FSPL;// +SF + CL;
	return PL_b;
}

// Расчёт газовых потерь(неправильно), пока не учитываем.

// Расчёт газовых потерь 
//double CalculatePathLossInGasses(double d, double f) 
//{
//	Engine* ep;
//	MATFile* pmat;
//	mxArray* pa_xData;
//	mxArray* pa_yData;
//	double* xData;
//	double* yData;
//	int numElements;
//
//	// Запускаем MATLAB Engine
//	ep = engOpen(NULL);
//	if (ep == NULL) {
//		std::cerr << "Не удалось запустить MATLAB Engine!" << std::endl;
//		return 1;
//	}
//
//	// Открываем .mat файл
//	pmat = matOpen("D:\\Всё моё\\Радиофак\\Практика\\РАДИО МОДУЛЬ\\ЗАДАНИЕ 2 (СПУТНИК)\\ChannelModel_NTN\\LSP\\gasAttTable.mat", "r");
//	if (pmat == NULL) {
//		std::cerr << "Не удалось открыть .mat файл!" << std::endl;
//		engClose(ep);
//		return 1;
//	}
//
//	// Читаем xData из .mat файла
//	pa_xData = matGetVariable(pmat, "xData");
//	if (pa_xData == NULL) {
//		std::cerr << "Не удалось прочитать xData из .mat файла!" << std::endl;
//		matClose(pmat);
//		engClose(ep);
//		return 1;
//	}
//
//	// Читаем yData из .mat файла
//	pa_yData = matGetVariable(pmat, "yData");
//	if (pa_yData == NULL) {
//		std::cerr << "Не удалось прочитать yData из .mat файла!" << std::endl;
//		mxDestroyArray(pa_xData); // Освобождаем память, выделенную для xData
//		matClose(pmat);
//		engClose(ep);
//		return 1;
//	}
//
//	// Получаем данные из mxArray для xData
//	numElements = mxGetNumberOfElements(pa_xData); // Должно быть 1000
//	xData = mxGetPr(pa_xData);
//
//	// Получаем данные из mxArray для yData
//	int numElementsY = mxGetNumberOfElements(pa_yData);  // Должно быть 1000
//	yData = mxGetPr(pa_yData);
//
//	// Проверка, что оба массива имеют одинаковую длину (опционально)
//	if (numElements != numElementsY) {
//		std::cerr << "Внимание: массивы xData и yData имеют разную длину!" << std::endl;
//	}
//
//	// Выводим несколько элементов для проверки
//	/*std::cout << "Первые 5 элементов xData:" << std::endl;
//	for (int i = 0; i < 5; ++i) {
//		std::cout << xData[i] << " ";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Первые 5 элементов yData:" << std::endl;
//	for (int i = 0; i < 5; ++i) {
//		std::cout << yData[i] << " ";
//	}
//	std::cout << std::endl;*/
//
//	double logFC = log10(f);
//	//std::cout << "logFC: " << logFC << std::endl;
//	int ind; // Переменная, которая ищет элемент из mat таблицы
//
//	for (int i = 0; i < 1000; ++i)
//	{
//		if (logFC > xData[i])
//		{
//			continue;
//		}
//		else
//		{
//			ind = i;
//			break;
//		}
//
//	}
//
//	//std::cout << "Index of element: " << ind << std::endl;
//	//std::cout << "xData[index]: " << xData[ind] << std::endl;
//	//std::cout << "yData[index]: " << yData[ind] << std::endl;
//
//	double PL_g;
//
//	PL_g = pow(10, yData[ind]) * d;
//
//	// Закрываем .mat файл и MATLAB Engine
//	mxDestroyArray(pa_xData);
//	mxDestroyArray(pa_yData);
//	matClose(pmat);
//	engClose(ep);
//
//	return PL_g;
//}


// Расчёт потерь при прохождении через здания


//Vector<double, 2> r = { 12.64, 28.19 };
//Vector<double, 2> s = { 3.72, -3.00 };
//Vector<double, 2> t = { 0.96, 8.48 };
//Vector<double, 2> u = { 9.6, 13.5 };
//Vector<double, 2> v = { 2.0, 3.8 };
//Vector<double, 2> w = { 9.1, 27.8 };
//Vector<double, 2> x = { -3.0, -2.9 };
//Vector<double, 2> y = { 4.5, 9.4 };
//Vector<double, 2> z = { -2.0, -2.1 };
//
//double C = -3.0;
//
//double CalculateLh(double r, double s, double f, double t)
//{
//	double Lh = r + s * log10(f) + t * pow(log10(f), 2);
//	return Lh;
//}
//
//double CalculateLe(double elevation)
//{
//	double Le = 0.212 * abs(elevation);
//	return Le;
//}
//
//Vector<double, 2> CalculateParametersForA(double Lh, double Le, double u, double v, double f)
//{
//	double mu_one = Lh + Le;
//	double sigma_one = u + v * log10(f);
//	Vector<double, 2> Apar = { mu_one, sigma_one };
//	return Apar;
//}
//
//Vector<double, 2> CalculateParametersForB(double w, double x, double y, double z, double f)
//{
//	double mu_two = w + x * log10(f);
//	double sigma_two = y + z * log10(f);
//	Vector<double, 2> Bpar = { mu_two, sigma_two };
//	return Bpar;
//}

//double CalculateInverseCummulativeNormalDistribution(double P)
//{
//
//}

//double CalculateBuildingEntryLoss(double A, double B, double C)
//{
//	double PL_e = 10 * log10(pow(10, 0.1 * A) + pow(10, 0.1 * B) + pow(10, 0.1 * C));
//	return PL_e;
//}





// PL = PL_b + PL_g + PL_s + PL_e
// PL_b = FSPL + SF + CL(alpha, f)

//double CalculateDistance(double Re, double h0, double alpha)
//{
//	double d = sqrt(pow(Re, 2) * pow(sin(alpha), 2) + pow(h0, 2) + 2 * h0 * Re) - Re * sin(alpha);
//	return d;
//}
//
//double GenerateSF(double std)
//{
//	double SF = RandomGenerators::generateGauss(0, std);
//	return SF;
//}
//
//double CalculatePathLoss(double d, double f, double SF, double CL)
//{
//	double PL = 0;
//	double FSPL = 32.45 + 20 * log10(f) + 20 * log10(d);
//	double PL_b = FSPL + SF + CL;
//	double PL_e = 0;
//	double PL_g = 0;
//	double PL_s = 0;
//	PL = PL_b + PL_g + PL_s + PL_e;
//
//	return PL;
//}

