#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

#include <vector>
#include <array>
#include <functional>
#include <cstdlib>

class IntegralConfiguration {
public:
	IntegralConfiguration(double start, double end, size_t n);
	double start, end;
	size_t n;
};

class GaussQuadrature {
public:
	GaussQuadrature();
	double integrate(std::function<double(double)> f, double a, double b, size_t n);
	double integrate(std::function<double(double)> f, IntegralConfiguration conf);

	std::vector<IntegralConfiguration> configurations = { { 0, 0.1, 8 }, { 0.1, 3., 30 }, { 3., 6., 12 } };

private:
	std::array<std::array<double, 21>, 21> roots;
	std::array<std::array<double, 21>, 21> weights;
};

//double gauss_legendre_integrate(std::function<double(double)> f, size_t n=6);
double gauss_legendre_integrate(std::function<double(double)> f);

#endif // GAUSSQUADRATURE_H
