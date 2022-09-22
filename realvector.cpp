#include "realvector.h"
#include <stdexcept>
#include <cmath>
#include <iostream>

RealVector::RealVector() {
}

RealVector::RealVector(std::vector<PhysicalDouble> coords) {
	this->coords = coords;
}

RealVector& RealVector::operator+=(RealVector rhs) {
	if (coords.size() != rhs.coords.size()) {
		throw std::invalid_argument("RealVector::Addition of RealVectors of different dimensions");
	}

	for (size_t i = 0; i < coords.size(); i++) {
		coords[i] += rhs[i];
	}

	return *this;
}

size_t RealVector::size() {
	return coords.size();
}

RealVector& RealVector::operator*=(PhysicalDouble rhs) {
	for (size_t i = 0; i < coords.size(); i++) {
		coords[i] *= rhs;
	}
	return *this;
}

RealVector& RealVector::operator-=(RealVector rhs) {
	return (*this) += ((-1) * rhs);
}

RealVector& RealVector::operator/=(PhysicalDouble rhs) {
	return (*this) *= (1 / rhs);
}

PhysicalDouble& RealVector::operator[](size_t i) {
	return coords[i];
}

std::vector<PhysicalDouble>::iterator RealVector::begin() {
	return coords.begin();
}

std::vector<PhysicalDouble>::iterator RealVector::end() {
	return coords.end();
}

PhysicalDouble RealVector::norm() {
	PhysicalDouble norm_sq = 0;
	for (auto &&c : coords) {
		norm_sq += c * c;
	}
	return std::sqrt(norm_sq/size());
}

RealVector operator +(RealVector lhs, RealVector rhs) {
	return lhs += rhs;
}

RealVector operator -(RealVector lhs, RealVector rhs) {
	return lhs -= rhs;
}

RealVector operator *(PhysicalDouble lhs, RealVector rhs) {
	return rhs *= lhs;
}

RealVector operator *(RealVector lhs, PhysicalDouble rhs) {
	return lhs *= rhs;
}

RealVector operator /(RealVector lhs, PhysicalDouble rhs) {
	return lhs /= rhs;
}

PhysicalDouble operator *(RealVector lhs, RealVector rhs) {
	if (lhs.size() != rhs.size()) {
		throw std::invalid_argument("RealVector::Product of RealVectors of different dimensions");
	}

	PhysicalDouble dot_prod = 0.;
	size_t half = lhs.size() / 2;
	if (lhs.size() % 2 == 0) {
		for (size_t i = 0; i < half; i++) {
			dot_prod += (lhs[half - i - 1] * rhs[half - i - 1] + lhs[half + i] * rhs[half + i]);
		}
	} else {
		dot_prod = lhs[half] * rhs[half];
		for (size_t i = 1; i <= half; i++) {
			dot_prod += (lhs[half + i] * rhs[half + i] + lhs[half - i] * rhs[half - i]);
		}
	}

	return dot_prod;
}

std::ostream& operator <<(std::ostream &out, RealVector v) {
	for (size_t i = 0; i < v.size(); i++) {
		out << v[i] << " ";
	}
	return out;
}

RealVector filter(RealVector v, std::vector<bool> filter_vec) {
	if (v.size() != filter_vec.size()) {
		throw std::invalid_argument("RealVector::Vector and filter have different dimensions");
	}

	for (size_t i = 0; i < v.size(); i++) {
		if (!filter_vec[i]) {
			v[i] = 0.;
		}
	}
	return v;
}

