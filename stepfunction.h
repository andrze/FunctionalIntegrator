#ifndef STEPFUNCTION_H
#define STEPFUNCTION_H
#include <vector>
#include <map>
#include "realvector.h"

class StepFunction {
public:
	StepFunction();
	StepFunction(PhysicalDouble step_size, std::vector<PhysicalDouble> y, PhysicalDouble domain_begin = 0);
	StepFunction(PhysicalDouble step_size, size_t num_points, std::function<PhysicalDouble(PhysicalDouble)> generator, PhysicalDouble domain_begin = 0);

	PhysicalDouble step_size;
	PhysicalDouble domain_begin = 0;
	PhysicalDouble domain_end();
	size_t num_points;
	std::vector<PhysicalDouble> vals;
	std::vector<PhysicalDouble> xs();
	StepFunction x_func();

	PhysicalDouble& operator[](size_t i);
	PhysicalDouble operator()(PhysicalDouble x);
	PhysicalDouble derivative(size_t n, size_t pos);
	StepFunction derivative(size_t n = 1);
	RealVector interval(size_t begin, size_t end);
	StepFunction zoom_in(PhysicalDouble begin, PhysicalDouble end);

	PhysicalDouble integral();
	PhysicalDouble norm();
	const std::pair<PhysicalDouble, PhysicalDouble> kappa_u();
	std::pair<PhysicalDouble, PhysicalDouble> minmax();
};

StepFunction operator+(StepFunction lhs, StepFunction rhs);
StepFunction operator*(StepFunction lhs, StepFunction rhs);
StepFunction operator-(StepFunction lhs, StepFunction rhs);
StepFunction operator/(StepFunction lhs, StepFunction rhs);

StepFunction operator+(PhysicalDouble lhs, StepFunction rhs);
StepFunction operator+(StepFunction lhs, PhysicalDouble rhs);
StepFunction operator-(PhysicalDouble lhs, StepFunction rhs);
StepFunction operator-(StepFunction lhs, PhysicalDouble rhs);
StepFunction operator*(PhysicalDouble lhs, StepFunction rhs);
StepFunction operator*(StepFunction lhs, PhysicalDouble rhs);
StepFunction operator/(PhysicalDouble lhs, StepFunction rhs);
StepFunction operator/(StepFunction lhs, PhysicalDouble rhs);

std::ostream& operator<<(std::ostream &out, StepFunction f);

#endif // STEPFUNCTION_H
