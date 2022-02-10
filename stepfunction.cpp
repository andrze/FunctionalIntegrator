#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "stepfunction.h"
#include "realvector.h"

const std::vector<PhysicalDouble> DER1_COEFFICIENTS { 1.l / 12, -2.l / 3, 0, 2.l / 3, -1.l / 12 };
const std::vector<PhysicalDouble> DER1_BORDER_COEFFICIENTS { -25.l / 12, 4.l, -3.l, 4.l / 3, -1.l / 4 };
const std::vector<PhysicalDouble> DER1_SEMIBORDER_COEFFICIENTS { -1.l / 4, -5.l / 6, 3.l / 2, -1.l / 2, 1.l / 12 };

const std::vector<PhysicalDouble> DER2_COEFFICIENTS = { -1.l / 12, 4.l / 3, -5.l / 2, 4.l / 3, -1.l / 12 };
const std::vector<PhysicalDouble> DER2_BORDER_COEFFICIENTS = { 35.l / 12, -26.l / 3, 19.l / 2, -14.l / 3, 11.l / 12 };
const std::vector<PhysicalDouble> DER2_SEMIBORDER_COEFFICIENTS = { 11.l / 12, -5.l / 3, 1.l / 2, 1.l / 3, -1.l / 12 };

PhysicalDouble lagrange_polynomial(PhysicalDouble x, std::vector<PhysicalDouble> x_values, std::vector<PhysicalDouble> y_values) {
	PhysicalDouble res = 0;

	for (size_t i = 0; i < x_values.size(); i++) {
		PhysicalDouble li = 1;
		for (size_t j = 0; j < x_values.size(); j++) {
			if (j == i) {
				continue;
			}
			li *= x - x_values[j];
			li /= x_values[i] - x_values[j];
		}
		res += li * y_values[i];
	}
	return res;
}

StepFunction::StepFunction() {
	step_size = 0;
	num_points = 0;
}

StepFunction::StepFunction(PhysicalDouble step_size, std::vector<PhysicalDouble> y, PhysicalDouble domain_begin) {
	this->domain_begin = domain_begin;
	this->step_size = step_size;
	vals = y;
	num_points = vals.size();
}

StepFunction::StepFunction(PhysicalDouble step_size, size_t num_points, std::function<PhysicalDouble(PhysicalDouble)> generator,
		PhysicalDouble domain_begin) {
	this->domain_begin = domain_begin;
	this->step_size = step_size;
	vals.reserve(num_points);
	for (size_t i = 0; i < num_points; i++) {
		vals.push_back(generator(domain_begin + PhysicalDouble(i) * step_size));
	}
	this->num_points = num_points;
}

PhysicalDouble StepFunction::domain_end() {
	return domain_begin + step_size * PhysicalDouble(vals.size() - 1);
}

std::vector<PhysicalDouble> StepFunction::xs() {
	std::vector<PhysicalDouble> xs;
	PhysicalDouble x = domain_begin;
	for (size_t i = 0; i < vals.size(); i++) {
		xs.push_back(x);
		x += step_size;
	}
	return xs;
}

StepFunction StepFunction::x_func() {
	std::vector<PhysicalDouble> xs = this->xs();
	return StepFunction(step_size, xs, domain_begin);
}

PhysicalDouble& StepFunction::operator[](size_t i) {

	if (i >= vals.size()) {
		throw std::invalid_argument("Argument out of function domain");
	}
	return vals[i];
}

PhysicalDouble StepFunction::operator()(PhysicalDouble x) {
	const int INTERPOLATION_ORDER = 5;

	size_t pos = size_t(std::floor((x - domain_begin) / step_size));

	if (x < domain_begin || x > domain_end() + step_size) {
		throw std::invalid_argument("Argument out of function domain");
	}

	int slice_begin = std::max(0, int(pos) - INTERPOLATION_ORDER / 2);
	if (slice_begin + INTERPOLATION_ORDER > int(num_points)) {
		slice_begin = int(num_points) - INTERPOLATION_ORDER;
	}

	std::vector<PhysicalDouble> x_values;
	for (int k = 0; k < INTERPOLATION_ORDER; k++) {
		x_values.push_back(PhysicalDouble(k + slice_begin) * step_size + domain_begin);
	}

	RealVector slice = interval(size_t(slice_begin), size_t(slice_begin + INTERPOLATION_ORDER));

	return lagrange_polynomial(x, x_values, slice.coords);
}

PhysicalDouble StepFunction::derivative(size_t n, size_t pos) {
	if (n == 0) {
		return vals[pos];
	}

	std::vector<PhysicalDouble> coefficients, border_coefficients, semiborder_coefficients;
	if (n == 1) {
		coefficients = DER1_COEFFICIENTS;
		border_coefficients = DER1_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER1_SEMIBORDER_COEFFICIENTS;
	} else if (n == 2) {
		coefficients = DER2_COEFFICIENTS;
		border_coefficients = DER2_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER2_SEMIBORDER_COEFFICIENTS;
	} else if (n > 2) {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	PhysicalDouble normalization = std::pow(step_size, -PhysicalDouble(n));

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
	std::vector<PhysicalDouble> new_vals;

	std::vector<PhysicalDouble> coefficients, border_coefficients, semiborder_coefficients;
	if (n == 1) {
		coefficients = DER1_COEFFICIENTS;
		border_coefficients = DER1_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER1_SEMIBORDER_COEFFICIENTS;
	} else if (n == 2) {
		coefficients = DER2_COEFFICIENTS;
		border_coefficients = DER2_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER2_SEMIBORDER_COEFFICIENTS;
	} else if (n > 2) {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	PhysicalDouble normalization = std::pow(step_size, -PhysicalDouble(n));
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

	return StepFunction(step_size, new_vals, domain_begin);
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
	auto interval_end = std::next(it, int(end));
	if (interval_end > vals.end()) {
		throw std::invalid_argument("Interval cannot reach outside arguments range");
	}

	auto interval_begin = std::next(it, int(begin));

	std::vector<PhysicalDouble> interval_vals(interval_begin, interval_end);
	if (reverse) {
		std::reverse(interval_vals.begin(), interval_vals.end());
	}

	return RealVector(interval_vals);
}

StepFunction StepFunction::zoom_in(PhysicalDouble begin, PhysicalDouble end) {
	std::vector<PhysicalDouble> new_vals;

	PhysicalDouble new_step = (end - begin) / PhysicalDouble(num_points - 1);

	for (size_t i = 0; i < num_points; i++) {
		PhysicalDouble new_x = begin + new_step * PhysicalDouble(i);

		new_vals.push_back((*this)(new_x));
	}

	return StepFunction(new_step, new_vals, begin);
}

PhysicalDouble StepFunction::integral() {
	PhysicalDouble sum = 0;

	for (auto &&v : vals) {
		sum += v;
	}

	return sum * step_size;
}

PhysicalDouble StepFunction::norm() {

	return std::sqrt(((*this) * (*this)).integral());
}

const std::pair<PhysicalDouble, PhysicalDouble> StepFunction::kappa_u() {
	PhysicalDouble x_lower = domain_begin;
	PhysicalDouble x_upper = domain_end();

	const size_t MAX_ITERATIONS = 100;
	const PhysicalDouble THRESHOLD = 1e-8;

	if ((*this)(x_lower) > 0 || (*this)(x_upper) < 0) {
		return {-1, -1};
	}

	for (size_t i = 0; i < MAX_ITERATIONS; i++) {

		PhysicalDouble x_new = (x_lower + x_upper) / 2;
		PhysicalDouble f_new = (*this)(x_new);

		if (f_new > 0) {
			x_upper = x_new;
		} else {
			x_lower = x_new;
		}

		if (std::abs(x_upper - x_lower) < THRESHOLD) {
			break;
		}
	}

	PhysicalDouble root = (x_upper + x_lower) / 2;

	return {root, derivative(1)(root)};
}

std::pair<PhysicalDouble, PhysicalDouble> StepFunction::minmax() {
	auto mini_maxi = std::minmax_element(vals.begin(), vals.end());
	return std::make_pair(*(mini_maxi.first), *(mini_maxi.second));
}

StepFunction operator+(StepFunction lhs, StepFunction rhs) {
	if (std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("Step sizes of added functions should be equal");
	}
	if (std::abs(lhs.domain_begin - rhs.domain_begin) > std::numeric_limits<PhysicalDouble>::epsilon()
			|| std::abs(lhs.domain_end() - rhs.domain_end()) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("Function domains are incompatible");
	}
	if (lhs.vals.size() != rhs.vals.size()) {
		throw std::invalid_argument("Values sizes of added functions should be equal");
	}
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < lhs.vals.size(); i++) {
		new_vals.push_back(lhs[i] + rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals, lhs.domain_begin);
}

StepFunction operator*(StepFunction lhs, StepFunction rhs) {
	if (std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("Step sizes of multiplied functions should be equal");
	}
	if (std::abs(lhs.domain_begin - rhs.domain_begin) > std::numeric_limits<PhysicalDouble>::epsilon()
			|| std::abs(lhs.domain_end() - rhs.domain_end()) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("Function domains are incompatible");
	}
	if (lhs.vals.size() != rhs.vals.size()) {
		throw std::invalid_argument("Values sizes of multiplied functions should be equal");
	}
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < lhs.vals.size(); i++) {
		new_vals.push_back(lhs[i] * rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals, lhs.domain_begin);
}

StepFunction operator-(StepFunction lhs, StepFunction rhs) {

	return lhs + (-1.) * rhs;
}

StepFunction operator/(StepFunction lhs, StepFunction rhs) {
	if (std::abs(lhs.step_size - rhs.step_size) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("Step sizes of multiplied functions should be equal");
	}
	if (std::abs(lhs.domain_begin - rhs.domain_begin) > std::numeric_limits<PhysicalDouble>::epsilon()
			|| std::abs(lhs.domain_end() - rhs.domain_end()) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("Function domains are incompatible");
	}
	if (lhs.vals.size() != rhs.vals.size()) {
		throw std::invalid_argument("Values sizes of multiplied functions should be equal");
	}
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < lhs.vals.size(); i++) {
		new_vals.push_back(lhs[i] / rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals, lhs.domain_begin);
}

StepFunction operator+(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs + rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator+(StepFunction lhs, PhysicalDouble rhs) {
	return rhs + lhs;
}

StepFunction operator-(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs - rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator-(StepFunction lhs, PhysicalDouble rhs) {
	return lhs + (-rhs);
}

StepFunction operator*(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs * rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator*(StepFunction lhs, PhysicalDouble rhs) {
	return rhs * lhs;
}

StepFunction operator/(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.vals.size(); i++) {
		new_vals.push_back(lhs / rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator/(StepFunction lhs, PhysicalDouble rhs) {
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
