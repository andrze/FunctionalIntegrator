/*
 * numericalmethods.cpp
 *
 *  Created on: Aug 31, 2022
 *      Author: andrzej
 */

#include "numericalmethods.h"

PhysicalDouble pow2(PhysicalDouble x) {
	return x * x;
}

PhysicalDouble lagrange_polynomial(PhysicalDouble x, RealVector x_values, RealVector y_values) {
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
