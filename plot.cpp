/*
 * plot.cpp
 *
 *  Created on: 15 gru 2018
 *      Author: andrzej
 */

#include "plot.h"
#include <stdexcept>
#include <cmath>

Plot::Plot() {

}

Plot::~Plot() {
}

size_t Plot::size() {
	size_t size = values.size();
	if (size != derivatives.size() || size != times.size()) {
		throw std::invalid_argument("Sizes of components in plot differ");
	}
	return size;
}

double Plot::exp_time_log_der(size_t k) {
	if (k >= size()) {
		throw std::invalid_argument("Argument k out of plot range");
	}
	return -times[k] * derivatives[k] / values[k];
}

std::array<double, 3> Plot::operator[](size_t k) {
	if (k >= size()) {
		throw std::invalid_argument("Argument k out of plot range");
	}
	std::array<double, 3> slice;
	slice[0] = times[k];
	slice[1] = values[k];
	slice[2] = derivatives[k];

	return slice;
}

void Plot::push(double t, double v, double d) {
	times.push_back(t);
	values.push_back(v);
	derivatives.push_back(d);
}

std::ostream& operator<<(std::ostream &out, Plot p) {
	for (size_t i = 0; i < p.size(); i++) {
		out << p.times[i] << ", " << p.values[i] << ", " << p.derivatives[i] << "\n";
	}
	return out;
}

