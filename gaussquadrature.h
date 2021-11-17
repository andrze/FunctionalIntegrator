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
	double partial_integrate(std::function<double(double)> f, size_t configuration_index) const;
	double integrate(std::function<double(double)> f) const;

	std::vector<IntegralConfiguration> configurations = { { 0, 0.1, 8 }, { 0.1, 3., 30 }, { 3., 6., 12 } };
	std::vector<std::vector<double> > evaluation_points;

private:
	std::array<std::array<double, 41>, 41> roots;
	std::array<std::array<double, 41>, 41> weights;
};

//double gauss_legendre_integrate(std::function<double(double)> f, size_t n=6);
//double gauss_legendre_integrate(std::function<double(double)> f);

#endif // GAUSSQUADRATURE_H
