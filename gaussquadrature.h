#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

#include <vector>
#include <array>
#include <functional>
#include <cstdlib>
#include "realvector.h"

class IntegralConfiguration {
public:
	IntegralConfiguration(PhysicalDouble start, PhysicalDouble end, size_t n);
	PhysicalDouble start, end;
	size_t n;
};

class GaussQuadrature {
public:
	GaussQuadrature();
	PhysicalDouble partial_integrate(std::function<PhysicalDouble(PhysicalDouble)> f, size_t configuration_index) const;
	PhysicalDouble partial_integrate(std::function<PhysicalDouble(PhysicalDouble,PhysicalDouble)> f, size_t configuration_index) const;
	PhysicalDouble integrate(std::function<PhysicalDouble(PhysicalDouble)> f) const;
	PhysicalDouble integrate(std::function<PhysicalDouble(PhysicalDouble,PhysicalDouble)> f) const;

	std::vector<IntegralConfiguration> configurations = { { 0, 0.1, 8 }, { 0.1, 3., 30 }, { 3., 6., 12 } };
	std::vector<std::vector<PhysicalDouble> > evaluation_points;
	std::vector<std::vector<PhysicalDouble> > expm1_points;

private:
	std::array<std::array<PhysicalDouble, 41>, 41> roots;
	std::array<std::array<PhysicalDouble, 41>, 41> weights;
};


#endif // GAUSSQUADRATURE_H
