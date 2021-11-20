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
#include "realvector.h"

class ButcherTable {
public:
	ButcherTable();
	ButcherTable(size_t steps, std::vector<PhysicalDouble> a_data, std::vector<PhysicalDouble> b_data);

	virtual ~ButcherTable();

	std::vector<std::vector<PhysicalDouble> > a;
	std::vector<PhysicalDouble> b;
	void runge_kutta_step(System& initial) const;
	bool usable=false;

};

std::ostream& operator<<(std::ostream &out, ButcherTable t);

extern const ButcherTable euler, rk4, ssp_rk3, ssp_rk4;

#endif /* RUNGEKUTTA_H_ */
