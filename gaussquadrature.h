#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

#include <vector>
#include <array>
#include <functional>
#include <cstdlib>
#include "realvector.h"

typedef std::array<PhysicalDouble, 3> IntegrandValue;
typedef std::array<PhysicalDouble, 6> IntegrandArgument;
typedef std::function<IntegrandValue(IntegrandArgument)> IntegrandFunction;

class IntegralConfiguration {
public:
	IntegralConfiguration(PhysicalDouble start, PhysicalDouble end, size_t n);
	PhysicalDouble start, end;
	size_t n;
};

class GaussQuadrature {
public:
	GaussQuadrature();
	GaussQuadrature(PhysicalDouble d, PhysicalDouble a);
	IntegrandValue partial_integrate(IntegrandFunction &f, size_t configuration_index) const;
	IntegrandValue integrate(IntegrandFunction &f) const;

private:
	std::array<std::array<PhysicalDouble, 41>, 41> roots;
	std::array<std::array<PhysicalDouble, 41>, 41> weights;
	std::vector<IntegralConfiguration> configurations;
	std::vector<std::vector<IntegrandArgument> > evaluation_points;

	PhysicalDouble d = 0, a = 0;
};

#endif // GAUSSQUADRATURE_H
