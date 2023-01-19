/*
 * rungekutta.cpp
 *
 *  Created on: Sep 23, 2021
 *      Author: andrzej
 *
 *      For reference: https://personal.math.ubc.ca/~cbm/mscthesis/cbm-mscthesis.pdf
 *      			   https://en.wikipedia.org/wiki/List_of_Runge–Kutta_methods
 */

#include <stdexcept>
#include "rungekutta.h"
#include "realvector.h"
#include <cmath>
#include <iomanip>
#include <iostream>

ButcherTable::ButcherTable() {

}

ButcherTable::ButcherTable(size_t steps, std::vector<PhysicalDouble> a_data, std::vector<PhysicalDouble> b_data,
	std::vector<PhysicalDouble> b_err_data) :
	b(b_data), b_err(b_err_data) {
	if (steps != b_data.size()) {
		throw std::invalid_argument("ButcherTable: Wrong size of Runge-Kutta initializing parameters");
	}
	if (b_err_data.size() == b_data.size()) {
		embedded = true;
		for (size_t i = 0; i < steps; i++) {
			b_trunc.push_back(b[i] - b_err[i]);
		}
	} else if (b_err_data.size() > 0) {
		throw std::invalid_argument("ButcherTable: Wrong size of Runge-Kutta initializing parameters");
	}

	size_t counter = 0;
	for (size_t i = 0; i < steps; i++) {
		a.push_back(std::vector<PhysicalDouble>());
		for (size_t j = 0; j < steps; j++) {
			if (j >= i) {
				a[i].push_back(0);
			} else {
				if (counter > a_data.size()) {
					throw std::invalid_argument("ButcherTable: Wrong size of Runge-Kutta initializing parameters");
				}
				a[i].push_back(a_data[counter]);
				counter++;
			}
		}
	}
	usable = true;
}

ButcherTable::~ButcherTable() {
}

void ButcherTable::runge_kutta_step(System &initial) const {
	if (!usable) {
		throw std::runtime_error("ButcherTable: Runge-Kutta method not properly initialized");
	}
	PhysicalDouble h = initial.delta_t;
	bool adaptive_step = embedded && h > low_delta_cutoff;
	std::vector<SystemDelta> derivatives;
	std::vector<PhysicalDouble> etas;
	derivatives.reserve(a.size());
	for (size_t i = 0; i < a.size(); i++) {
		System point = initial;
		for (size_t j = 0; j < i; j++) {
			point += h * a[i][j] * derivatives[j];
			//point.fix_v();
		}
		if (i == 0) { // The point of this if statement is to execute time_derivative on System& initial
			derivatives.push_back(initial.time_derivative());
			etas.push_back(initial.eta);
		} else {
			derivatives.push_back(point.time_derivative());
			etas.push_back(point.eta);
		}
	}
	SystemDelta total_der = b[0] * derivatives[0];
	SystemDelta truncation_err;

	if (adaptive_step) {
		truncation_err = b_trunc[0] * derivatives[0];
	}
	PhysicalDouble z_der = b[0] * etas[0];
	for (size_t i = 1; i < b.size(); i++) {
		total_der += b[i] * derivatives[i];
		z_der += b[i] * etas[i];
		if (adaptive_step) {
			truncation_err += b_trunc[i] * derivatives[i];
		}
	}

	PhysicalDouble h_new = 0;
	if (adaptive_step) {
		RealVector truncation_error_vector = truncation_err.vector_representation(),
			initial_vector = initial.full_vector_representation();
		PhysicalDouble trunc_error = 0;
		for(size_t i=0; i < truncation_error_vector.size(); i++){
			if(std::abs(initial_vector[i]) > 1e+2 ){
				trunc_error += std::pow(truncation_error_vector[i]/initial_vector[i],2);
			} else {
				trunc_error += std::pow(truncation_error_vector[i],2);
			}
		}
		trunc_error = h * std::sqrt(trunc_error/truncation_error_vector.size());
		h_new = 0.9l * h * std::pow(error_tolerance / trunc_error, 0.2l);
		if (h_new < low_delta_cutoff) {
			h_new = low_delta_cutoff;
			std::cout << "Could not reach expected precision with adaptive step, defaulting to minimum allowed step "
				<< std::scientific << low_delta_cutoff << std::fixed << std::endl;
		}

		if (trunc_error > error_tolerance) {
			initial.delta_t = h_new;
			runge_kutta_step(initial);
			return;
		}
	}
	initial += h * total_der;
	initial.z_dim *= (1 + h * z_der);

	//initial.fix_v();
	initial.time += h;
	initial.step++;
	if (adaptive_step) {
		initial.delta_t = h_new;
	}
}

std::ostream& operator<<(std::ostream &out, ButcherTable t) {

	for (size_t i = 0; i < t.a.size(); i++) {
		for (size_t j = 0; j < i; j++) {
			out << t.a[i][j] << ' ';
		}
		out << '\n';
	}
	for (size_t j = 0; j < t.b.size(); j++) {
		out << t.b[j] << ' ';
	}
	out << '\n';

	if (t.embedded) {
		for (size_t j = 0; j < t.b_err.size(); j++) {
			out << t.b_err[j] << ' ';
		}
		out << '\n';
	}

	return out;
}

//Euler
const ButcherTable euler { 1, { 0 }, { 1.l } };

//RK4
const ButcherTable rk4 { 4, { .5l, 0, .5l, 0, 0, 1 }, { 1.l / 6, 1.l / 3, 1.l / 3, 1.l / 6 } };

//SSP RK3
std::array<PhysicalDouble, 3> a_ssprk3_a { { 1, 1.l / 4, 1.l / 4 } };
std::array<PhysicalDouble, 3> b_ssprk3_a { { 1.l / 6, 1.l / 6, 2.l / 3 } };
const ButcherTable ssp_rk3 { 3, { 1, 1.l / 4, 1.l / 4 }, { 1.l / 6, 1.l / 6, 2.l / 3 } };

//SSP RK(5,4)
const ButcherTable ssp_rk4 { 5, { 0.39175l, 0.21767l, 0.36841l, 0.082692l, 0.13996l, 0.25189l, 0.067966l, 0.11503l,
	0.20703l, 0.54497l }, { 0.14681l, 0.24848l, 0.10426l, 0.27444l, 0.22601l } };

//Runge–Kutta–Fehlberg
const ButcherTable rkf { 6, { 0.25l, 3.l / 32, 9.l / 32, 1932.l / 2197, -7200.l / 2197, 7296.l / 2197, 439.l / 216, -8,
	3680.l / 513, -845.l / 4104, -8.l / 27, 2, -3544.l / 2565, 1859.l / 4104, -11.l / 40 }, { 16.l / 135, 0, 6656.l
	/ 12825, 28561.l / 56430, -9.l / 50, 2.l / 55 }, { 25.l / 216, 0, 1408.l / 2565, 2197.l / 4104, -1.l / 5, 0 } };

