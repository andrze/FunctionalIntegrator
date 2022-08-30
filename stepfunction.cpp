#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "stepfunction.h"
#include "realvector.h"

//const std::vector<PhysicalDouble> DER1_COEFFICIENTS { 1.l / 12, -2.l / 3, 0, 2.l / 3, -1.l / 12 }; // 5 point first derivative
//const std::vector<PhysicalDouble> DER1_BORDER_COEFFICIENTS { -25.l / 12, 4.l, -3.l, 4.l / 3, -1.l / 4 }; // 5 point first derivative
//const std::vector<PhysicalDouble> DER1_SEMIBORDER_COEFFICIENTS { -1.l / 4, -5.l / 6, 3.l / 2, -1.l / 2, 1.l / 12 }; // 5 point first derivative

//const std::vector<PhysicalDouble> DER1_BORDER_COEFFICIENTS { -(137.l / 60), 5.l, -5.l, 10.l / 3, -(5.l / 4), 1.l / 5 }; // 6 point first derivative
//const std::vector<PhysicalDouble> DER1_SEMIBORDER_COEFFICIENTS { -(1.l / 5), -(13.l / 12), 2.l, -1.l, 1.l / 3, -(1.l
//	/ 20) }; // 6 point first derivative
//const std::vector<PhysicalDouble> DER1_SEMISEMIBORDER_COEFFICIENTS = { 1.l / 20, -(1.l / 2), -(1.l / 3), 1.l,
//	-(1.l / 4), 1.l / 30 }; // 6 point first derivative

const std::vector<PhysicalDouble> DER1_COEFFICIENTS { -(1.l / 60), 3.l / 20, -(3.l / 4), 0, 3.l / 4, -(3.l / 20), 1.l
	/ 60 }; // 7 point first derivative
const std::vector<PhysicalDouble> DER1_BORDER_COEFFICIENTS { -(49.l / 20), 6.l, -(15.l / 2), 20.l / 3, -(15.l / 4), 6.l
	/ 5, -(1.l / 6) }; // 7 point first derivative
const std::vector<PhysicalDouble> DER1_SEMIBORDER_COEFFICIENTS { -(1.l / 6), -(77.l / 60), 5.l / 2, -(5.l / 3), 5.l / 6,
	-(1.l / 4), 1.l / 30 }; // 7 point first derivative
const std::vector<PhysicalDouble> DER1_SEMISEMIBORDER_COEFFICIENTS = { 1.l / 30, -(2.l / 5), -(7.l / 12), 4.l / 3, -(1.l
	/ 2), 2.l / 15, -(1.l / 60) }; // 7 point first derivative

//const std::vector<PhysicalDouble> DER2_COEFFICIENTS = { -1.l / 12, 4.l / 3, -5.l / 2, 4.l / 3, -1.l / 12 }; // 5/6 point second derivative
//const std::vector<PhysicalDouble> DER2_BORDER_COEFFICIENTS = { 15.l / 4, -77.l / 6, 107.l / 6, -13.l, 61.l / 12, -5.l/ 6 }; // 6 point second derivative
//const std::vector<PhysicalDouble> DER2_SEMIBORDER_COEFFICIENTS = { 5.l / 6, -5.l / 4, -1.l / 3, 7.l / 6, -1.l / 2, 1.l/ 12 }; // 6 point second derivative
//const std::vector<PhysicalDouble> DER2_SEMISEMIBORDER_COEFFICIENTS = {5.l/6, -(5.l/4), -(1.l/3), 7.l/6, -(1.l/2), 1.l/12}; // 6 point second derivative

//const std::vector<PhysicalDouble> DER2_COEFFICIENTS = { -1.l / 12, 4.l / 3, -5.l / 2, 4.l / 3, -1.l / 12 }; // 5 point second derivative
//const std::vector<PhysicalDouble> DER2_BORDER_COEFFICIENTS = { 35.l / 12, -26.l / 3, 19.l / 2, -14.l / 3, 11.l / 12 }; // 5 point second derivative
//const std::vector<PhysicalDouble> DER2_SEMIBORDER_COEFFICIENTS = { 11.l / 12, -5.l / 3, 1.l / 2, 1.l / 3, -1.l / 12 }; // 5 point second derivative
//const std::vector<PhysicalDouble> DER2_SEMISEMIBORDER_COEFFICIENTS =	{ -1.l / 12, 4.l / 3, -5.l / 2, 4.l / 3, -1.l / 12 }; // 5 point second derivative

const std::vector<PhysicalDouble> DER2_COEFFICIENTS = { 1.l / 90, -(3.l / 20), 3.l / 2, -(49.l / 18), 3.l / 2, -(3.l
	/ 20), 1.l / 90 }; // 7 point second derivative
const std::vector<PhysicalDouble> DER2_BORDER_COEFFICIENTS = { 203.l / 45, -(87.l / 5), 117.l / 4, -(254.l / 9), 33.l
	/ 2, -(27.l / 5), 137.l / 180 }; // 7 point second derivative
const std::vector<PhysicalDouble> DER2_SEMIBORDER_COEFFICIENTS = { 137.l / 180, -(49.l / 60), -(17.l / 12), 47.l / 18,
	-(19.l / 12), 31.l / 60, -(13.l / 180) }; // 7 point second derivative
const std::vector<PhysicalDouble> DER2_SEMISEMIBORDER_COEFFICIENTS = { -(13.l / 180), 19.l / 15, -(7.l / 3), 10.l / 9,
	1.l / 12, -(1.l / 15), 1.l / 90 }; // 7 point second derivative

const std::vector<std::vector<PhysicalDouble> > ALL_COEFS { DER1_COEFFICIENTS, DER1_BORDER_COEFFICIENTS,
	DER1_SEMIBORDER_COEFFICIENTS, DER1_SEMISEMIBORDER_COEFFICIENTS, DER2_COEFFICIENTS, DER2_BORDER_COEFFICIENTS,
	DER2_SEMIBORDER_COEFFICIENTS, DER2_SEMISEMIBORDER_COEFFICIENTS };

PhysicalDouble lagrange_polynomial(PhysicalDouble x, std::vector<PhysicalDouble> x_values,
	std::vector<PhysicalDouble> y_values) {
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
	num_points_derivatives = 0;
	for (auto &&c : ALL_COEFS) {
		num_points_derivatives = std::max(num_points_derivatives, c.size());
	}
}

StepFunction::StepFunction(PhysicalDouble step_size, std::vector<PhysicalDouble> y, PhysicalDouble domain_begin) {
	this->domain_begin = domain_begin;
	this->step_size = step_size;
	vals = y;
	num_points = vals.size();

	num_points_derivatives = 0;
	for (auto &&c : ALL_COEFS) {
		num_points_derivatives = std::max(num_points_derivatives, c.size());
	}
}

StepFunction::StepFunction(PhysicalDouble step_size, size_t num_points,
	std::function<PhysicalDouble(PhysicalDouble)> generator, PhysicalDouble domain_begin) {
	this->domain_begin = domain_begin;
	this->step_size = step_size;
	vals.reserve(num_points);
	for (size_t i = 0; i < num_points; i++) {
		vals.push_back(generator(domain_begin + PhysicalDouble(i) * step_size));
	}
	this->num_points = num_points;

	num_points_derivatives = 0;
	for (auto &&c : ALL_COEFS) {
		num_points_derivatives = std::max(num_points_derivatives, c.size());
	}
}

void StepFunction::check_compatibility(StepFunction f) {
	if (std::abs(step_size - f.step_size) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("StepFunction: Operation on StepFunctions of non-equal step sizes");
	}
	if (std::abs(domain_begin - f.domain_begin) > std::numeric_limits<PhysicalDouble>::epsilon()
		|| std::abs(domain_end() - f.domain_end()) > std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::invalid_argument("StepFunction: Operation on StepFunctions with incompatible domains");
	}

	if (num_points != f.num_points) {
		throw std::invalid_argument("StepFunction: Operation on StepFunctions with different numbers of points");
	}
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
	const int INTERPOLATION_ORDER = 7;

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

	PhysicalDouble inverse_step = 1 / step_size;
	PhysicalDouble normalization = 1;
	std::vector<PhysicalDouble> coefficients, border_coefficients, semiborder_coefficients, semisemiborder_coefficients;
	if (n == 1) {
		coefficients = DER1_COEFFICIENTS;
		border_coefficients = DER1_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER1_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER1_SEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step;
	} else if (n == 2) {
		coefficients = DER2_COEFFICIENTS;
		border_coefficients = DER2_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER2_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER2_SEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step * inverse_step;
	} else if (n > 2) {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	RealVector coefs(coefficients), border_coefs(border_coefficients), semiborder_coefs(semiborder_coefficients),
		semisemiborder_coefs(semisemiborder_coefficients);

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
	} else if (pos == 2) {
		slice = interval(pos - 2, pos - 2 + semisemiborder_coefs.size());
		used_coefs = semisemiborder_coefs;
	} else if (pos == vals.size() - 3) {
		slice = interval(pos + 2, pos + 2 - semisemiborder_coefs.size());
		used_coefs = semisemiborder_coefs;
		if (n % 2 == 1)
			used_coefs *= -1;
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

	PhysicalDouble inverse_step = 1 / step_size;
	PhysicalDouble normalization = 1;
	std::vector<PhysicalDouble> coefficients, border_coefficients, semiborder_coefficients, semisemiborder_coefficients;
	if (n == 1) {
		coefficients = DER1_COEFFICIENTS;
		border_coefficients = DER1_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER1_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER1_SEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step;
	} else if (n == 2) {
		coefficients = DER2_COEFFICIENTS;
		border_coefficients = DER2_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER2_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER2_SEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step * inverse_step;
	} else if (n > 2) {
		throw std::invalid_argument("Derivative of order larger than 2 are not implemented\n");
	}

	for (auto &c : coefficients) {
		c *= normalization;
	}
	for (auto &c : border_coefficients) {
		c *= normalization;
	}
	for (auto &c : semiborder_coefficients) {
		c *= normalization;
	}
	for (auto &c : semisemiborder_coefficients) {
		c *= normalization;
	}

	RealVector coefs(coefficients), border_coefs(border_coefficients), semiborder_coefs(semiborder_coefficients),
		semisemiborder_coefs(semisemiborder_coefficients);

	for (size_t i = 0; i < vals.size(); i++) {
		RealVector used_coefs;
		RealVector slice;

		if (i > 2 && i < vals.size() - 3) {
			used_coefs = coefs;
			size_t size = coefs.size();
			slice = interval(i - size / 2, i + 1 + size / 2);
		} else if (i == 0) {
			slice = interval(0, border_coefs.size());
			used_coefs = border_coefs;
		} else if (i == 1) {
			slice = interval(0, semiborder_coefs.size());
			used_coefs = semiborder_coefs;
		} else if (i == 2) {
			slice = interval(0, semisemiborder_coefs.size());
			used_coefs = semisemiborder_coefs;
		} else if (i == vals.size() - 3) {
			slice = interval(vals.size(), vals.size() - semisemiborder_coefs.size());
			used_coefs = semisemiborder_coefs;
			if (n % 2 == 1)
				used_coefs *= -1;
		} else if (i == vals.size() - 2) {
			slice = interval(vals.size(), vals.size() - semiborder_coefs.size());
			used_coefs = semiborder_coefs;
			if (n % 2 == 1)
				used_coefs *= -1;
		} else if (i == vals.size() - 1) {
			slice = interval(vals.size(), vals.size() - border_coefs.size());
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

StepFunction StepFunction::cut_domain(size_t begin, size_t end) {
	std::vector<PhysicalDouble> new_vals(vals.begin() + int(begin), vals.begin() + int(end));

	PhysicalDouble new_domain_begin = domain_begin + step_size * int(begin);
	return StepFunction(step_size, new_vals, new_domain_begin);
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
	const PhysicalDouble THRESHOLD = 1e-16;

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
			std::cout << i << '\n';
			break;
		}
	}

	PhysicalDouble root = (x_upper + x_lower) / 2;

	return {root, derivative(1)(root)};
}

const std::pair<PhysicalDouble, PhysicalDouble> StepFunction::kappa_u(StepFunction der) {
	PhysicalDouble x_lower = domain_begin;
	PhysicalDouble x_upper = domain_end();
	PhysicalDouble guess = (2 * domain_end() - domain_begin) / 3;

	const size_t MAX_ITERATIONS = 100;
	const PhysicalDouble THRESHOLD = 1e-18;

	if ((*this)(x_lower) > 0 || (*this)(x_upper) < 0) {
		return {-1, -1};
	}
	PhysicalDouble damping = .001;

	for (size_t i = 0; i < MAX_ITERATIONS; i++) {
		PhysicalDouble val = (*this)(guess);
		PhysicalDouble der_val = der(guess);

		PhysicalDouble new_guess = guess - (1 - damping) * val / der_val;

		if (std::abs(new_guess - guess) < THRESHOLD) {
			guess = new_guess;
			break;
		}
		if(new_guess > x_lower && new_guess < x_upper){
			guess = new_guess;
		} else {
			damping *= 2;
		}
	}

	return {guess, der(guess)};
}

std::pair<PhysicalDouble, PhysicalDouble> StepFunction::minmax() {
	auto mini_maxi = std::minmax_element(vals.begin(), vals.end());
	return std::make_pair(*(mini_maxi.first), *(mini_maxi.second));
}

std::vector<PhysicalDouble>::iterator StepFunction::begin() {
	return vals.begin();
}

std::vector<PhysicalDouble>::iterator StepFunction::end() {
	return vals.end();
}

StepFunction operator+(StepFunction lhs, StepFunction rhs) {
	lhs.check_compatibility(rhs);
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < lhs.num_points; i++) {
		new_vals.push_back(lhs[i] + rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals, lhs.domain_begin);
}

StepFunction operator*(StepFunction lhs, StepFunction rhs) {
	lhs.check_compatibility(rhs);
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < lhs.num_points; i++) {
		new_vals.push_back(lhs[i] * rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals, lhs.domain_begin);
}

StepFunction operator-(StepFunction lhs, StepFunction rhs) {

	return lhs + (-1.) * rhs;
}

StepFunction operator/(StepFunction lhs, StepFunction rhs) {
	lhs.check_compatibility(rhs);
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < lhs.num_points; i++) {
		new_vals.push_back(lhs[i] / rhs[i]);
	}

	return StepFunction(lhs.step_size, new_vals, lhs.domain_begin);
}

StepFunction operator+(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.num_points; i++) {
		new_vals.push_back(lhs + rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator+(StepFunction lhs, PhysicalDouble rhs) {
	return rhs + lhs;
}

StepFunction operator-(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.num_points; i++) {
		new_vals.push_back(lhs - rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator-(StepFunction lhs, PhysicalDouble rhs) {
	return lhs + (-rhs);
}

StepFunction operator*(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.num_points; i++) {
		new_vals.push_back(lhs * rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator*(StepFunction lhs, PhysicalDouble rhs) {
	return rhs * lhs;
}

StepFunction operator/(PhysicalDouble lhs, StepFunction rhs) {
	std::vector<PhysicalDouble> new_vals;

	for (size_t i = 0; i < rhs.num_points; i++) {
		new_vals.push_back(lhs / rhs[i]);
	}

	return StepFunction(rhs.step_size, new_vals, rhs.domain_begin);
}

StepFunction operator/(StepFunction lhs, PhysicalDouble rhs) {
	rhs = 1 / rhs;

	return lhs * rhs;
}

std::ostream& operator<<(std::ostream &out, StepFunction f) {
	for (auto v = f.begin(); v != f.end(); std::advance(v, 1)) {
		out << *v << ' ';
	}
	out << std::endl;
	return out;
}
