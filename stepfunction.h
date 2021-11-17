#ifndef STEPFUNCTION_H
#define STEPFUNCTION_H
#include <vector>
#include <map>
#include "realvector.h"

class StepFunction {
public:
	StepFunction();
	StepFunction(double step_size, std::vector<double> y, double domain_begin = 0);
	StepFunction(double step_size, size_t num_points, std::function<double(double)> generator, double domain_begin = 0);

	double step_size;
	double domain_begin = 0;
	double domain_end();
	size_t num_points;
	std::vector<double> vals;
	std::vector<double> xs();
	StepFunction x_func();

	double& operator[](size_t i);
	double operator()(double x);
	double derivative(size_t n, size_t pos);
	StepFunction derivative(size_t n = 1);
	RealVector interval(size_t begin, size_t end);
	StepFunction zoom_in(double begin, double end);

	double integral();
	double norm();
	const std::pair<double, double> kappa_u();
	std::pair<double, double> minmax();
};

StepFunction operator+(StepFunction lhs, StepFunction rhs);
StepFunction operator*(StepFunction lhs, StepFunction rhs);
StepFunction operator-(StepFunction lhs, StepFunction rhs);
StepFunction operator/(StepFunction lhs, StepFunction rhs);

StepFunction operator+(double lhs, StepFunction rhs);
StepFunction operator+(StepFunction lhs, double rhs);
StepFunction operator-(double lhs, StepFunction rhs);
StepFunction operator-(StepFunction lhs, double rhs);
StepFunction operator*(double lhs, StepFunction rhs);
StepFunction operator*(StepFunction lhs, double rhs);
StepFunction operator/(double lhs, StepFunction rhs);
StepFunction operator/(StepFunction lhs, double rhs);

std::ostream& operator<<(std::ostream &out, StepFunction f);

#endif // STEPFUNCTION_H
