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

ButcherTable::ButcherTable() {

}

ButcherTable::ButcherTable(size_t steps, std::vector<PhysicalDouble> a_data, std::vector<PhysicalDouble> b_data) :
		b(b_data) {
	if (steps != b_data.size()) {
		throw std::invalid_argument("ButcherTable: Wrong size of Runge-Kutta initializing parameters");
	}
	size_t counter = 0;
	for (size_t i = 0; i < steps; i++) {
		a.push_back(std::vector<PhysicalDouble>());
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

void ButcherTable::runge_kutta_step(System &initial) const {
	if (!usable) {
		throw std::runtime_error("ButcherTable: Runge-Kutta method not properly initialized");
	}
	PhysicalDouble h = initial.delta_t;
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

	double z_der = b[0] * derivatives[0].eta;
	derivatives[0] *= b[0];
	for (size_t i = 0; i < b.size(); i++) {
		derivatives[0] += b[i] * derivatives[i];
		z_der += b[i] * derivatives[i].eta;
	}
	initial += h * derivatives[0];
	initial.z_dim *= (1 + h * z_der);

	//initial.fix_v();
	initial.time += h;
	initial.step++;
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
const ButcherTable euler { 1, { 0 }, { 1.l } };

//RK4
const ButcherTable rk4 { 4, { .5l, 0, .5l, 0, 0, 1 }, { 1.l / 6, 1.l / 3, 1.l / 3, 1.l / 6 } };

//SSP RK3
std::array<PhysicalDouble, 3> a_ssprk3_a { { 1, 1.l / 4, 1.l / 4 } };
std::array<PhysicalDouble, 3> b_ssprk3_a { { 1.l / 6, 1.l / 6, 2.l / 3 } };
const ButcherTable ssp_rk3 { 3, { 1, 1.l / 4, 1.l / 4 }, { 1.l / 6, 1.l / 6, 2.l / 3 } };

//SSP RK(5,4)
const ButcherTable ssp_rk4 { 5, { 0.39175l, 0.21767l, 0.36841l, 0.082692l, 0.13996l, 0.25189l, 0.067966l, 0.11503l, 0.20703l,
		0.54497l }, { 0.14681l, 0.24848l, 0.10426l, 0.27444l, 0.22601l } };

