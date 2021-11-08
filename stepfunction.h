#ifndef STEPFUNCTION_H
#define STEPFUNCTION_H
#include "realvector.h"
#include <vector>
#include <map>

class StepFunction {
public:
	StepFunction(double step_size, std::vector<double> y);
	StepFunction(double step_size, size_t num_points, std::function<double(double)> generator);

	double step_size;
	double domain;
	size_t num_points;
	std::vector<double> vals;
	std::vector<double> xs();

	double& operator[](size_t i);
	double operator()(double x);
	double derivative(size_t n, size_t pos);
	StepFunction derivative(size_t n = 1);
	RealVector interval(size_t begin, size_t end);

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
