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
	GaussQuadrature(PhysicalDouble d);
	//PhysicalDouble partial_integrate(std::function<PhysicalDouble(PhysicalDouble)> f, size_t configuration_index) const;
	PhysicalDouble partial_integrate(std::function<PhysicalDouble(std::array<PhysicalDouble,6>)>& f,
			size_t configuration_index) const;
	//PhysicalDouble integrate(std::function<PhysicalDouble(PhysicalDouble)> f) const;
	PhysicalDouble integrate(std::function<PhysicalDouble(std::array<PhysicalDouble,6>)>& f) const;

private:
	std::array<std::array<PhysicalDouble, 41>, 41> roots;
	std::array<std::array<PhysicalDouble, 41>, 41> weights;
	std::vector<IntegralConfiguration> configurations = { { 0, 0.1, 25 }, { 0.1, 3., 40 }, { 3., 6., 20 } };
	std::vector<std::vector<std::array<PhysicalDouble,6> > > evaluation_points;

	PhysicalDouble d;
};

#endif // GAUSSQUADRATURE_H
