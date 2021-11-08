#include "stepfunction.h"
#include <utility>
#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iterator>
#include <array>
/*
static const std::array<double, 5> der1_coefficients { { 1. / 12, -2. / 3, 0, 2. / 3, -1. / 12 } };
static const std::array<double, 5> der1_border_coefficients { { -25. / 12, 4., -3., 4. / 3, -1. / 4 } };
static const std::array<double, 5> der1_semiborder_coefficients { { -1. / 4, -5. / 6, 3. / 2, -1. / 2, 1. / 12 } };
static const std::array<double, 5> der2_coefficients { { -1. / 12, 4. / 3, -5. / 2, 4. / 3, -1. / 12 } };

static const std::array<double, 6> der2_border_coefficients { { 15. / 4, -77. / 6, 107. / 6, -13., 61. / 12, -5. / 6 } };
static const std::array<double, 6> der2_semiborder_coefficients { { 5. / 6, -5. / 4, -1. / 3, 7. / 6, -1. / 2, 1. / 12 } };
*/
static const std::array<double, 3> der1_coefficients { { -1./2, 0., 1./2 } };
static const std::array<double, 3> der1_border_coefficients { { -3./2, 2., -1/2. } };
static const std::array<double, 3> der2_coefficients { { 1., -2., 1. } };

static const std::array<double, 4> der2_border_coefficients { { 2., -5, 4, -1. } };

StepFunction::StepFunction(double step_size, std::vector<double> y) {
	this->step_size = step_size;
	vals = y;
	domain = step_size * double(vals.size());
	num_points = vals.size();
}

StepFunction::StepFunction(double step_size, size_t num_points, std::function<double(double)> generator) {
	this->step_size = step_size;
	domain = step_size * double(num_points);
	vals.reserve(num_points);
	for (size_t i = 0; i < num_points; i++) {
		vals.push_back(generator(double(i) * step_size));
	}
	this->num_points = num_points;
}

std::vector<double> StepFunction::xs() {
	std::vector<double> xs;
	double x = 0;
	for (size_t i = 0; i < vals.size(); i++) {
		xs.push_back(x);
		x += step_size;
	}
	return xs;
}

double& StepFunction::operator[](size_t i) {

	if (i > vals.size()) {
		throw std::invalid_argument("Argument out of function domain");
	}
	return vals[i];
}

double StepFunction::operator()(double x) {
	size_t pos = size_t(std::floor(x / step_size));
	if (pos > vals.size()) {
		throw std::invalid_argument("Argument out of function domain");
	}

	double remainder = x / step_size - double(pos);

	return vals[pos] * (1 - remainder) + vals[pos + 1] * remainder;
}
/*
double StepFunction::derivative(size_t n, size_t pos) {
	if(pos>=num_points){
		throw std::invalid_argument("The argument of the derivative lies out of the function domain.");
	}
	std::vector<double> coefficients, border_coefficients, semiborder_coefficients;
	if (n == 0) {
		return vals[pos];
	} else if (n == 1) {
		coefficients = std::vector<double>(der1_coefficients.begin(), der1_coefficients.end());
		border_coefficients = std::vector<double>(der1_border_coefficients.begin(), der1_border_coefficients.end());
		semiborder_coefficients = std::vector<double>(der1_semiborder_coefficients.begin(),
				der1_semiborder_coefficients.end());
	} else if (n == 2) {
		coefficients = std::vector<double>(der2_coefficients.begin(), der2_coefficients.end());
		border_coefficients = std::vector<double>(der2_border_coefficients.begin(), der2_border_coefficients.end());
		semiborder_coefficients = std::vector<double>(der2_semiborder_coefficients.begin(),
				der2_semiborder_coefficients.end());
	} else {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	double normalization = std::pow(step_size, -double(n));

	RealVector coefs(coefficients), border_coefs(border_coefficients), semiborder_coefs(semiborder_coefficients);

	RealVector used_coefs;
	RealVector slice;

	if (pos > 1 && pos < vals.size() - 2) {
		used_coefs = coefs;
		size_t size = coefs.size();
		slice = interval(pos - size / 2, pos + 1 + size / 2);
	} else if (pos == 0) {
		slice = interval(pos, pos + border_coefs.size());
		used_coefs = border_coefs;
	} else if (pos == 1) {
		slice = interval(pos - 1, pos - 1 + semiborder_coefs.size());
		used_coefs = semiborder_coefs;
	} else if (pos == vals.size() - 2) {
		slice = interval(pos + 2, pos + 2 - semiborder_coefs.size());
		used_coefs = semiborder_coefs;
		if (n % 2 == 1)
			used_coefs *= -1;
	} else if (pos == vals.size() - 1) {
		slice = interval(pos + 1, pos + 1 - border_coefs.size());
		used_coefs = border_coefs;
		if (n % 2 == 1)
			used_coefs *= -1;
	} else {
		throw std::runtime_error("");
	}

	return slice * used_coefs * normalization;
}

StepFunction StepFunction::derivative(size_t n) {
	if (n == 0) {
		return *this;
	}
	std::vector<double> new_vals;

	std::vector<double> coefficients, border_coefficients, semiborder_coefficients;
	if (n == 1) {
		coefficients = std::vector<double>(der1_coefficients.begin(), der1_coefficients.end());
		border_coefficients = std::vector<double>(der1_border_coefficients.begin(), der1_border_coefficients.end());
		semiborder_coefficients = std::vector<double>(der1_semiborder_coefficients.begin(),
				der1_semiborder_coefficients.end());
	} else if (n == 2) {
		coefficients = std::vector<double>(der2_coefficients.begin(), der2_coefficients.end());
		border_coefficients = std::vector<double>(der2_border_coefficients.begin(), der2_border_coefficients.end());
		semiborder_coefficients = std::vector<double>(der2_semiborder_coefficients.begin(),
				der2_semiborder_coefficients.end());
	} else {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	double normalization = std::pow(step_size, -double(n));
	for (auto &c : coefficients) {
		c *= normalization;
	}
	for (auto &c : border_coefficients) {
		c *= normalization;
	}
	for (auto &c : semiborder_coefficients) {
		c *= normalization;
	}

	RealVector coefs(coefficients), border_coefs(border_coefficients), semiborder_coefs(semiborder_coefficients);

	for (size_t i = 0; i < vals.size(); i++) {
		RealVector used_coefs;
		RealVector slice;

		if (i > 1 && i < vals.size() - 2) {
			used_coefs = coefs;
			size_t size = coefs.size();
			slice = interval(i - size / 2, i + 1 + size / 2);
		} else if (i == 0) {
			slice = interval(i, i + border_coefs.size());
			used_coefs = border_coefs;
		} else if (i == 1) {
			slice = interval(i - 1, i - 1 + semiborder_coefs.size());
			used_coefs = semiborder_coefs;
		} else if (i == vals.size() - 2) {
			slice = interval(i + 2, i + 2 - semiborder_coefs.size());
			used_coefs = semiborder_coefs;
			if (n % 2 == 1)
				used_coefs *= -1;
		} else if (i == vals.size() - 1) {
			slice = interval(i + 1, i + 1 - border_coefs.size());
			used_coefs = border_coefs;
			if (n % 2 == 1)
				used_coefs *= -1;
		} else {
			throw std::runtime_error("");
		}

		new_vals.push_back(slice * used_coefs);
	}

	return StepFunction(step_size, new_vals);
}*/

double StepFunction::derivative(size_t n, size_t pos) {
	std::vector<double> coefficients, border_coefficients;
	if (n == 0) {
		return vals[pos];
	} else if (n == 1) {
		coefficients = std::vector<double>(der1_coefficients.begin(), der1_coefficients.end());
		border_coefficients = std::vector<double>(der1_border_coefficients.begin(), der1_border_coefficients.end());
	} else if (n == 2) {
		coefficients = std::vector<double>(der2_coefficients.begin(), der2_coefficients.end());
		border_coefficients = std::vector<double>(der2_border_coefficients.begin(), der2_border_coefficients.end());
	} else {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	double normalization = std::pow(step_size, -double(n));

	RealVector coefs(coefficients), border_coefs(border_coefficients);

	RealVector used_coefs;
	RealVector slice;

	if (pos > 1 && pos < vals.size() - 2) {
		used_coefs = coefs;
		size_t size = coefs.size();
		slice = interval(pos - size / 2, pos + 1 + size / 2);
	} else if (pos == 0) {
		slice = interval(pos, pos + border_coefs.size());
		used_coefs = border_coefs;
	} else if (pos == vals.size() - 1) {
		slice = interval(pos + 1, pos + 1 - border_coefs.size());
		used_coefs = border_coefs;
		if (n % 2 == 1)
			used_coefs *= -1;
	} else {
		throw std::runtime_error("");
	}

	return slice * used_coefs * normalization;
}

StepFunction StepFunction::derivative(size_t n) {
	if (n == 0) {
		return *this;
	}
	std::vector<double> new_vals;

	std::vector<double> coefficients, border_coefficients;
	if (n == 1) {
		coefficients = std::vector<double>(der1_coefficients.begin(), der1_coefficients.end());
		border_coefficients = std::vector<double>(der1_border_coefficients.begin(), der1_border_coefficients.end());
	} else if (n == 2) {
		coefficients = std::vector<double>(der2_coefficients.begin(), der2_coefficients.end());
		border_coefficients = std::vector<double>(der2_border_coefficients.begin(), der2_border_coefficients.end());
	} else {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	double normalization = std::pow(step_size, -double(n));
	for (auto &c : coefficients) {
		c *= normalization;
	}
	for (auto &c : border_coefficients) {
		c *= normalization;
	}

	RealVector coefs(coefficients), border_coefs(border_coefficients);

	for (size_t i = 0; i < vals.size(); i++) {
		RealVector used_coefs;
		RealVector slice;

		if (i > 0 && i < vals.size() - 1) {
			used_coefs = coefs;
			size_t size = coefs.size();
			slice = interval(i - size / 2, i + 1 + size / 2);
		} else if (i == 0) {
			slice = interval(i, i + border_coefs.size());
			used_coefs = border_coefs;
		} else if (i == vals.size() - 1) {
			slice = interval(i + 1, i + 1 - border_coefs.size());
			used_coefs = border_coefs;
			if (n % 2 == 1)
				used_coefs *= -1;
		} else {
			throw std::runtime_error("");
		}

		new_vals.push_back(slice * used_coefs);
	}

	return StepFunction(step_size, new_vals);
}

RealVector StepFunction::interval(size_t begin, size_t end) {
	bool reverse = false;
	if (end < begin) {
		reverse = true;
		size_t temp = begin;
		begin = end;
		end = temp;
	}

	auto it = vals.begin();
	auto interval_end = it + int(end);
	if (interval_end > vals.end()) {
		throw std::invalid_argument("Interval cannot reach outside arguments range");
	}

	auto interval_begin = it + int(begin);

	std::vector<double> interval_vals(interval_begin, interval_end);
	if (reverse) {
		std::reverse(interval_vals.begin(), interval_vals.end());
	}

	return RealVector(interval_vals);
}

double StepFunction::integral() {
	double sum = 0;

	for (auto &&v : vals) {
		sum += v;
	}

	return sum * step_size;
}

double StepFunction::norm() {

	return std::sqrt(((*this) * (*this)).integral());
}

const std::pair<double, double> StepFunction::kappa_u() {
	int l_bound = 0, u_bound = int(vals.size() - 1);

	while (vals[size_t(l_bound)] > 0.) {
		l_bound++;
	}

	int guess;
	size_t i = 0;
	while (vals[size_t(l_bound + 1)] <= 0) {
		double f_diff = (vals[size_t(u_bound)] - vals[size_t(l_bound)]);

		double a = f_diff / (u_bound - l_bound);
		guess = int(l_bound - vals[size_t(l_bound)] / a);
		if (vals[size_t(guess)] < 0.) {
			if (guess <= l_bound) {
				l_bound += 1;
			} else {
				l_bound = guess;
			}
		} else {
			if (guess >= u_bound) {
				if (vals[size_t(u_bound) - 1] >= 0) {
					u_bound -= 1;
				} else {
					l_bound = u_bound - 1;
				}
			} else {
				u_bound = guess;
			}
		}
		i++;
		if (i > 100) {
			break;
		}
	}

	auto interpolate = [](double y1, double y2, double x) {
		return y1 * (1 - x) + y2 * x;
	};

	double x = vals[size_t(l_bound)] / (vals[size_t(l_bound)] - vals[size_t(l_bound + 1)]);
	double kappa = step_size * (l_bound + x);

	double u1 = derivative(1, size_t(l_bound));
	double u2 = derivative(1, size_t(l_bound + 1));
	double u = interpolate(u1, u2, x);

	return std::make_pair(kappa, u);
}

std::pair<double, double> StepFunction::minmax() {
	auto mini_maxi = std::minmax_element(vals.begin(), vals.end());
	return std::make_pair(*(mini_maxi.first), *(mini_maxi.second));
}

StepFunction operator+(StepFunction lhs, StepFunction rhs) {
	if (std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<double>::epsilon()) {
		throw std::invalid_argument("Step sizes of added functions should be equal");
	}
	if (lhs.vals.size() != rhs.vals.size()) {
		throw std::invalid_argument("Values sizes of added functions should be equal");
	}
	std::vector<double> new_vals;

	for (size_t i = 0; i < lhs.vals.size(); i++) {
		new_vals.push_back(lhs[i] + rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals);
}

StepFunction operator*(StepFunction lhs, StepFunction rhs) {
	if (std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<double>::epsilon()) {
		throw std::invalid_argument("Step sizes of multiplied functions should be equal");
	}
	if (lhs.vals.size() != rhs.vals.size()) {
		throw std::invalid_argument("Values sizes of multiplied functions should be equal");
	}
	std::vector<double> new_vals;

	for (size_t i = 0; i < lhs.vals.size(); i++) {
		new_vals.push_back(lhs[i] * rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals);
}

StepFunction operator-(StepFunction lhs, StepFunction rhs) {

	return lhs + (-1.) * rhs;
}

StepFunction operator/(StepFunction lhs, StepFunction rhs) {
	if (std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<double>::epsilon()) {
		throw std::invalid_argument("Step sizes of multiplied functions should be equal");
	}
	if (lhs.vals.size() != rhs.vals.size()) {
		throw std::invalid_argument("Values sizes of multiplied functions should be equal");
	}
	std::vector<double> new_vals;

	for (size_t i = 0; i < lhs.vals.size(); i++) {
		new_vals.push_back(lhs[i] / rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals);
}

StepFunction operator+(double lhs, StepFunction rhs) {
	std::vector<double> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs + rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals);
}

StepFunction operator+(StepFunction lhs, double rhs) {
	return rhs + lhs;
}

StepFunction operator-(double lhs, StepFunction rhs) {
	std::vector<double> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs - rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals);
}

StepFunction operator-(StepFunction lhs, double rhs) {
	return lhs + (-rhs);
}

StepFunction operator*(double lhs, StepFunction rhs) {
	std::vector<double> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs * rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals);
}

StepFunction operator*(StepFunction lhs, double rhs) {
	return rhs * lhs;
}

StepFunction operator/(double lhs, StepFunction rhs) {
	std::vector<double> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs / rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals);
}

StepFunction operator/(StepFunction lhs, double rhs) {
	rhs = 1 / rhs;

	return lhs * rhs;
}

std::ostream& operator<<(std::ostream &out, StepFunction f) {
	for (auto &&v : f.vals) {
		out << v << ' ';
	}
	out << std::endl;
	return out;
}

