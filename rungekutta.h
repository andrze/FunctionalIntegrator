/*
 * rungekutta.h
 *
 *  Created on: Sep 23, 2021
 *      Author: andrzej
 */

#ifndef RUNGEKUTTA_H_
#define RUNGEKUTTA_H_

#include <vector>
#include <iostream>
#include "system.h"

class ButcherTable {
public:
	ButcherTable();
	ButcherTable(size_t steps, std::vector<double> a_data, std::vector<double> b_data);

	virtual ~ButcherTable();

	std::vector<std::vector<double> > a;
	std::vector<double> b;
	System runge_kutta_step(System& initial) const;
	bool usable=false;

};

std::ostream& operator<<(std::ostream &out, ButcherTable t);

extern const ButcherTable euler, rk4, ssp_rk3, ssp_rk4;

#endif /* RUNGEKUTTA_H_ */
