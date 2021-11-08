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

ButcherTable::ButcherTable() {

}

ButcherTable::ButcherTable(size_t steps, std::vector<double> a_data, std::vector<double> b_data) :
		b(b_data) {
	if (steps != b_data.size()) {
		throw std::invalid_argument("Wrong size of Runge-Kutta initializing parameters");
	}
	size_t counter = 0;
	for (size_t i = 0; i < steps; i++) {
		a.push_back(std::vector<double>());
		for (size_t j = 0; j < steps; j++) {
			if (j >= i) {
				a[i].push_back(0);
			} else {
				a[i].push_back(a_data[counter]);
				counter++;
			}
		}
	}
	usable = true;
}

ButcherTable::~ButcherTable() {
}

System ButcherTable::runge_kutta_step(System &initial) const {
	if (!usable) {
		throw std::runtime_error("Runge-Kutta method not properly initialized");
	}
	double h = initial.delta_t;
	std::vector<System> derivatives;
	derivatives.reserve(a.size());
	for (size_t i = 0; i < a.size(); i++) {
		System point = initial;
		for (size_t j = 0; j < i; j++) {
			point += h * a[i][j] * derivatives[j];
			//point.fix_v();
		}
		if (i == 0) { // The point of this if statement is to execute time_derivative on System& initial
			derivatives.push_back(initial.time_derivative());
		} else {
			derivatives.push_back(point.time_derivative());
		}
	}

	for (size_t i = 0; i < b.size(); i++) {
		initial += h * b[i] * derivatives[i];
		initial.z_dim *= (1 + h * b[i] * derivatives[i].eta);
	}

	//initial.fix_v();
	initial.time += h;
	initial.step++;

	return initial;
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

	return out;
}

//Euler
const ButcherTable euler { 1, { 0 }, { 1. } };

//RK4
const ButcherTable rk4 { 4, { .5, 0, .5, 0, 0, 1 }, { 1. / 6, 1. / 3, 1. / 3, 1. / 6 } };

//SSP RK3
std::array<double, 3> a_ssprk3_a { { 1, 1. / 4, 1. / 4 } };
std::array<double, 3> b_ssprk3_a { { 1. / 6, 1. / 6, 2. / 3 } };
const ButcherTable ssp_rk3 { 3, { 1, 1. / 4, 1. / 4 }, { 1. / 6, 1. / 6, 2. / 3 } };

//SSP RK(5,4)
const ButcherTable ssp_rk4 { 5, { 0.39175, 0.21767, 0.36841, 0.082692, 0.13996, 0.25189, 0.067966, 0.11503, 0.20703,
		0.54497 }, { 0.14681, 0.24848, 0.10426, 0.27444, 0.22601 } };

