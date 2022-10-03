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
	ButcherTable(size_t steps, std::vector<PhysicalDouble> a_data, std::vector<PhysicalDouble> b_data,
		std::vector<PhysicalDouble> b_err_data={});

	virtual ~ButcherTable();

	bool embedded=false;
	std::vector<std::vector<PhysicalDouble> > a;
	std::vector<PhysicalDouble> b;
	std::vector<PhysicalDouble> b_err;
	std::vector<PhysicalDouble> b_trunc;
	void runge_kutta_step(System &initial) const;
	bool usable = false;
	PhysicalDouble error_tolerance=1e-16, delta_cutoff = 1e-6;

};

std::ostream& operator<<(std::ostream &out, ButcherTable t);

extern const ButcherTable euler, rk4, ssp_rk3, ssp_rk4, rkf;

#endif /* RUNGEKUTTA_H_ */
