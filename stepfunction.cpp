#include <limits>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include "stepfunction.h"
#include "realvector.h"
#include "numericalmethods.h"

const std::vector<PhysicalDouble> DER1_BORDER_COEFFICIENTS = { -761.l / 280, 8.l / 1, -14.l / 1, 56.l / 3, -35.l / 2,
	56.l / 5, -14.l / 3, 8.l / 7, -1.l / 8 };
const std::vector<PhysicalDouble> DER1_SEMIBORDER_COEFFICIENTS = { -1.l / 8, -223.l / 140, 7.l / 2, -7.l / 2, 35.l / 12,
	-7.l / 4, 7.l / 10, -1.l / 6, 1.l / 56 };
const std::vector<PhysicalDouble> DER1_SEMISEMIBORDER_COEFFICIENTS = { 1.l / 56, -2.l / 7, -19.l / 20, 2.l / 1, -5.l
	/ 4, 2.l / 3, -1.l / 4, 2.l / 35, -1.l / 168 };
const std::vector<PhysicalDouble> DER1_SEMISEMISEMIBORDER_COEFFICIENTS = { -1.l / 168, 1.l / 14, -1.l / 2, -9.l / 20,
	5.l / 4, -1.l / 2, 1.l / 6, -1.l / 28, 1.l / 280 };
const std::vector<PhysicalDouble> DER1_CENTRAL_COEFFICIENTS = { 1.l / 280, -4.l / 105, 1.l / 5, -4.l / 5, 0.l / 1, 4.l
	/ 5, -1.l / 5, 4.l / 105, -1.l / 280 };

const std::vector<PhysicalDouble> DER2_BORDER_COEFFICIENTS = { 29531.l / 5040, -962.l / 35, 621.l / 10, -4006.l / 45,
	691.l / 8, -282.l / 5, 2143.l / 90, -206.l / 35, 363.l / 560 };
const std::vector<PhysicalDouble> DER2_SEMIBORDER_COEFFICIENTS = { 363.l / 560, 8.l / 315, -83.l / 20, 153.l / 20,
	-529.l / 72, 47.l / 10, -39.l / 20, 599.l / 1260, -29.l / 560 };
const std::vector<PhysicalDouble> DER2_SEMISEMIBORDER_COEFFICIENTS = { -29.l / 560, 39.l / 35, -331.l / 180, 1.l / 5,
	9.l / 8, -37.l / 45, 7.l / 20, -3.l / 35, 47.l / 5040 };
const std::vector<PhysicalDouble> DER2_SEMISEMISEMIBORDER_COEFFICIENTS = { 47.l / 5040, -19.l / 140, 29.l / 20, -118.l
	/ 45, 11.l / 8, -1.l / 20, -7.l / 180, 1.l / 70, -1.l / 560 };
const std::vector<PhysicalDouble> DER2_CENTRAL_COEFFICIENTS = { -1.l / 560, 8.l / 315, -1.l / 5, 8.l / 5, -205.l / 72,
	8.l / 5, -1.l / 5, 8.l / 315, -1.l / 560 };

const std::vector<PhysicalDouble> DER3_BORDER_COEFFICIENTS = { -801.l / 80, 349.l / 6, -18353.l / 120, 2391.l / 10,
	-1457.l / 6, 4891.l / 30, -561.l / 8, 527.l / 30, -469.l / 240 };
const std::vector<PhysicalDouble> DER3_SEMIBORDER_COEFFICIENTS = { -469.l / 240, 303.l / 40, -731.l / 60, 269.l / 24,
	-57.l / 8, 407.l / 120, -67.l / 60, 9.l / 40, -1.l / 48 };
const std::vector<PhysicalDouble> DER3_SEMISEMIBORDER_COEFFICIENTS = { -1.l / 48, -53.l / 30, 273.l / 40, -313.l / 30,
	103.l / 12, -9.l / 2, 197.l / 120, -11.l / 30, 3.l / 80 };
const std::vector<PhysicalDouble> DER3_SEMISEMISEMIBORDER_COEFFICIENTS = { 3.l / 80, -43.l / 120, -5.l / 12, 147.l / 40,
	-137.l / 24, 463.l / 120, -27.l / 20, 7.l / 24, -7.l / 240 };
const std::vector<PhysicalDouble> DER3_CENTRAL_COEFFICIENTS = { -7.l / 240, 3.l / 10, -169.l / 120, 61.l / 30, 0.l / 1,
	-61.l / 30, 169.l / 120, -3.l / 10, 7.l / 240 };

const std::vector<PhysicalDouble> DER4_BORDER_COEFFICIENTS = { 1069.l / 80, -1316.l / 15, 15289.l / 60, -2144.l / 5,
	10993.l / 24, -4772.l / 15, 2803.l / 20, -536.l / 15, 967.l / 240 };
const std::vector<PhysicalDouble> DER4_SEMIBORDER_COEFFICIENTS = { 967.l / 240, -229.l / 10, 3439.l / 60, -2509.l / 30,
	631.l / 8, -1489.l / 30, 1219.l / 60, -49.l / 10, 127.l / 240 };
const std::vector<PhysicalDouble> DER4_SEMISEMIBORDER_COEFFICIENTS = { 127.l / 240, -11.l / 15, -77.l / 20, 193.l / 15,
	-407.l / 24, 61.l / 5, -311.l / 60, 19.l / 15, -11.l / 80 };
const std::vector<PhysicalDouble> DER4_SEMISEMISEMIBORDER_COEFFICIENTS = { -11.l / 80, 53.l / 30, -341.l / 60, 77.l
	/ 10, -107.l / 24, 11.l / 30, 13.l / 20, -7.l / 30, 7.l / 240 };
const std::vector<PhysicalDouble> DER4_CENTRAL_COEFFICIENTS = { 7.l / 240, -2.l / 5, 169.l / 60, -122.l / 15, 91.l / 8,
	-122.l / 15, 169.l / 60, -2.l / 5, 7.l / 240 };

const std::vector<std::vector<PhysicalDouble> > ALL_COEFS { DER1_CENTRAL_COEFFICIENTS, DER1_BORDER_COEFFICIENTS,
	DER1_SEMIBORDER_COEFFICIENTS, DER1_SEMISEMIBORDER_COEFFICIENTS, DER1_SEMISEMISEMIBORDER_COEFFICIENTS,
	DER2_CENTRAL_COEFFICIENTS, DER2_BORDER_COEFFICIENTS, DER2_SEMIBORDER_COEFFICIENTS, DER2_SEMISEMIBORDER_COEFFICIENTS,
	DER2_SEMISEMISEMIBORDER_COEFFICIENTS, DER3_CENTRAL_COEFFICIENTS, DER3_BORDER_COEFFICIENTS,
	DER3_SEMIBORDER_COEFFICIENTS, DER3_SEMISEMIBORDER_COEFFICIENTS, DER3_SEMISEMISEMIBORDER_COEFFICIENTS,
	DER4_CENTRAL_COEFFICIENTS, DER4_BORDER_COEFFICIENTS, DER4_SEMIBORDER_COEFFICIENTS, DER4_SEMISEMIBORDER_COEFFICIENTS,
	DER4_SEMISEMISEMIBORDER_COEFFICIENTS };

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
		throw std::invalid_argument("StepFunction: Argument out of function domain");
	}
	return vals[i];
}

PhysicalDouble StepFunction::operator()(PhysicalDouble x) {
	const int INTERPOLATION_ORDER = 7;

	size_t pos = size_t(std::floor((x - domain_begin) / step_size));

	if (x < domain_begin || x > domain_end() + step_size) {
		throw std::invalid_argument("StepFunction: Argument out of function domain");
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

	return lagrange_polynomial(x, x_values, slice);
}

PhysicalDouble StepFunction::derivative(size_t n, size_t pos) {
	if (n == 0) {
		return vals[pos];
	}

	PhysicalDouble inverse_step = 1 / step_size;
	PhysicalDouble normalization = 1;
	std::vector<PhysicalDouble> coefficients, border_coefficients, semiborder_coefficients, semisemiborder_coefficients,
		semisemisemiborder_coefficients;
	if (n == 1) {
		coefficients = DER1_CENTRAL_COEFFICIENTS;
		border_coefficients = DER1_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER1_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER1_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER1_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step;
	} else if (n == 2) {
		coefficients = DER2_CENTRAL_COEFFICIENTS;
		border_coefficients = DER2_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER2_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER2_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER2_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step * inverse_step;
	} else if (n == 3) {
		coefficients = DER3_CENTRAL_COEFFICIENTS;
		border_coefficients = DER3_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER3_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER3_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER3_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = std::pow(inverse_step, 3);
	} else if (n == 4) {
		coefficients = DER4_CENTRAL_COEFFICIENTS;
		border_coefficients = DER4_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER4_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER4_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER4_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = std::pow(inverse_step, 4);
	} else {
		throw std::invalid_argument("StepFunction: Derivatives of order larger than 4 are not implemented\n");
	}

	RealVector coefs(coefficients), border_coefs(border_coefficients), semiborder_coefs(semiborder_coefficients),
		semisemiborder_coefs(semisemiborder_coefficients), semisemisemiborder_coefs(semisemisemiborder_coefficients);

	RealVector used_coefs;
	RealVector slice;

	if (pos > 2 && pos < vals.size() - 3) {
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
	} else if (pos == 3) {
		slice = interval(pos - 3, pos - 3 + semisemisemiborder_coefs.size());
		used_coefs = semisemisemiborder_coefs;
	} else if (pos == vals.size() - 4) {
		slice = interval(pos + 3, pos + 3 - semisemisemiborder_coefs.size());
		used_coefs = semisemisemiborder_coefs;
		if (n % 2 == 1)
			used_coefs *= -1;
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
	std::vector<PhysicalDouble> coefficients, border_coefficients, semiborder_coefficients, semisemiborder_coefficients,
		semisemisemiborder_coefficients;
	if (n == 1) {
		coefficients = DER1_CENTRAL_COEFFICIENTS;
		border_coefficients = DER1_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER1_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER1_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER1_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step;
	} else if (n == 2) {
		coefficients = DER2_CENTRAL_COEFFICIENTS;
		border_coefficients = DER2_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER2_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER2_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER2_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = inverse_step * inverse_step;
	} else if (n == 3) {
		coefficients = DER3_CENTRAL_COEFFICIENTS;
		border_coefficients = DER3_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER3_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER3_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER3_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = std::pow(inverse_step, 3);
	} else if (n == 4) {
		coefficients = DER4_CENTRAL_COEFFICIENTS;
		border_coefficients = DER4_BORDER_COEFFICIENTS;
		semiborder_coefficients = DER4_SEMIBORDER_COEFFICIENTS;
		semisemiborder_coefficients = DER4_SEMISEMIBORDER_COEFFICIENTS;
		semisemisemiborder_coefficients = DER4_SEMISEMISEMIBORDER_COEFFICIENTS;
		normalization = std::pow(inverse_step, 4);
	} else {
		throw std::invalid_argument("StepFunction: Derivatives of order larger than 4 are not implemented\n");
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
	for (auto &c : semisemisemiborder_coefficients) {
		c *= normalization;
	}

	RealVector coefs(coefficients), border_coefs(border_coefficients), semiborder_coefs(semiborder_coefficients),
		semisemiborder_coefs(semisemiborder_coefficients), semisemisemiborder_coefs(semisemisemiborder_coefficients);

	for (size_t i = 0; i < vals.size(); i++) {
		RealVector used_coefs;
		RealVector slice;

		if (i > 3 && i < vals.size() - 4) {
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
		} else if (i == 3) {
			slice = interval(0, semisemisemiborder_coefs.size());
			used_coefs = semisemisemiborder_coefs;
		} else if (i == vals.size() - 4) {
			slice = interval(vals.size(), vals.size() - semisemisemiborder_coefs.size());
			used_coefs = semisemisemiborder_coefs;
			if (n % 2 == 1)
				used_coefs *= -1;
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
		throw std::invalid_argument("StepFunction: Requested interval extends outside arguments range");
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
		if (new_guess > x_lower && new_guess < x_upper) {
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
