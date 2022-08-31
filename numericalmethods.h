/*
 * numericalmethods.h
 *
 *  Created on: Aug 31, 2022
 *      Author: andrzej
 */

#ifndef NUMERICALMETHODS_H_
#define NUMERICALMETHODS_H_

#include "realvector.h"

PhysicalDouble pow2(PhysicalDouble x);

PhysicalDouble lagrange_polynomial(PhysicalDouble x, RealVector x_values, RealVector y_values);

#endif /* NUMERICALMETHODS_H_ */
