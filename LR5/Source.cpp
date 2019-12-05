#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <iomanip>
#include <fstream>

std::ofstream noise("noisedFunction.txt");
std::ofstream filter("filteredFunction.txt");

double originalFunction(double x) // Оригинальная функция
{
	return sin(x) + 0.5;
}

int experimentCounter(double p, double epsilon) // Подсчет необходимого числа экспериментов
{
	int n = (log(1 - p) / log(1 - (epsilon / M_PI)) + 1);
	return n;
}

std::vector<double> noisedFunction(std::vector<double> x) // Зашумленная функция
{
	std::vector<double> noised;
	std::mt19937 engine(std::random_device{}());
	for (int i = 0; i < x.size(); i++)
	{
		auto generator = std::uniform_real_distribution<double>(-0.25, 0.25);
		double alpha = generator(engine);
		noised.push_back(originalFunction(x[i]) + alpha);
	}
	for (int i = 0; i < noised.size(); i++)
	{
		noise << '(' << x[i] << ";" << noised[i] << ')' << std::endl;
	}
	return noised;
}

std::vector<double> vectorX(double xmin, double xmax, int K) // Исходный вектор аргументов
{
	std::vector<double> x;
	x.resize(K);
	for (int k = 0; k < K; k++)
	{
		x[k] = xmin + k * (xmax - xmin) / K;
	}
	return x;
}

std::vector<double> alphaGenerator(int r) // Генерация вектора альфа
{
	std::mt19937 engine(std::random_device{}());
	if (r == 3)
	{
		std::vector<double> alpha(3);
		auto generator = std::uniform_real_distribution<double>(0, 1);
		alpha[1] = generator(engine);
		alpha[0] = alpha[2] = 0.5*(1 - alpha[1]);
		return alpha;
	}
	else
	{
		std::vector<double> alpha(5);
		auto generator = std::uniform_real_distribution<double>(0, 1);
		alpha[2] = generator(engine);
		generator = std::uniform_real_distribution<double>(0, 1 - alpha[2]);
		alpha[1] = alpha[3] = 0.5 * generator(engine);
		double a = 0;
		for (int i = 1; i < 4; i++)
		{
			a += alpha[i];
		}
		alpha[0] = alpha[4] = 0.5 * (1 - a);
		return alpha;
	}
	
}

std::vector<double> filtration(int K, int M, std::vector<double> function, std::vector<double> alpha, std::vector<double> x) // Фильтрация
{
	std::vector<double> filteredFunction;
	double f = 1;
	for (int i = 0; i < K; i++)
	{
		for (int j = i - M + 1; j <= i + M - 1; j++)
		{
			f *= pow(function[i], alpha[j + M - i]);
		}
		filteredFunction.push_back(f);
		f = 1;
	}
	for (int i = 0; i < filteredFunction.size(); i++)
	{
		filter << "(" << x[i] << ";" << filteredFunction[i] << ")" << std::endl;
	}
	return filteredFunction;
}

double omega(std::vector<double> function, int K) // Омега
{
	double maximum = 0;
	for (int i = 1; i < K; i++)
	{
		maximum = std::max(maximum, abs(function[i - 1] - function[i]));
	}
	return maximum;
}

double delta(std::vector<double> originalFunction, std::vector<double> noisedFunction, int K) // Дельта
{
	double maximum = 0;
	for (int i = 0; i < K; i++)
	{
		maximum = std::max(maximum, abs(noisedFunction[i] - originalFunction[i]));
	}
	return maximum;
}

double distance(double omega, double delta) // Расстояние
{
	return std::max(omega, delta);
}

std::pair<std::vector<double>, std::vector<double>> optimalForLambda(int n, int lambda, int K, int r, std::vector<double> x,std::vector<double> noised) // Поиск оптимального решения для определенного лямбда
{
	double minimum = 100;
	std::pair<std::vector<double>, std::vector<double>> optimal;
	for (int i = 0; i < n; i++)
	{
		std::vector<double> alpha = alphaGenerator(r);
		double w = omega(filtration(K, (r - 1) / 2, noised, alpha, x), K);
		double d = delta(noised, filtration(K, (r - 1) / 2, noised, alpha, x), K);
		double J1 = lambda * w;
		double J2 = (1 - lambda) * d;
		if (J1 + J2 < minimum)
		{
			optimal.first = alpha;
			optimal.second.push_back(J1 + J2);
			optimal.second.push_back(w);
			optimal.second.push_back(d);
		}
	}
	return optimal;
}

double lambdaOptimal(int L, int K, int r, std::vector<double> x) // Поиск оптимального значения лямбда
{
	double minimum = 100000000;
	double lOpt = -1;
	std::vector<double> noised = noisedFunction(x);
	for (int i = 0; i <= L; i++)
	{
		double lambda = i * 1.0 / L;
		std::pair<std::vector<double>, std::vector<double>> optimal;
		optimal = optimalForLambda(experimentCounter(0.95, 0.01), lambda, K, r, x, noised);
		if (optimal.second[0] < minimum)
		{
			minimum = optimal.second[0];
			lOpt = lambda;
		}
		std::cout << std::fixed << std::setprecision(4) << lambda << "\t|" << optimal.second[0] << "\t| [ ";
		for (int j = 0; j < optimal.first.size(); j++)
		{
			std::cout << optimal.first[j] << ", ";
		}
		std::cout << "]\t|" << optimal.second[1] << "\t|" << optimal.second[2] << "\n";
	}
	return lOpt;
}


int main()
{
	std::cout << "lambda\t|dis\t|alpha\t\t\t\t|w\t|d\n";
	std::cout << "Optimal lambda: " << lambdaOptimal(10, 100, 3, vectorX(0, M_PI, 100)) << "\n\n\n";
	std::cout << "lambda\t|dis\t|alpha\t\t\t\t\t\t|w\t|d\n";
	std::cout << "Optimal lambda: " << lambdaOptimal(10, 100, 5, vectorX(0, M_PI, 100)) << std::endl;
	system("pause");
	return 0;
}
