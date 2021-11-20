/*
 * Integrator.h
 *
 *  Created on: Sep 9, 2021
 *      Author: andrzej
 */

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_
#include <vector>
#include <array>
#include <queue>
#include "system.h"
#include "stepfunction.h"
#include "plot.h"
#include "gaussquadrature.h"
#include "rungekutta.h"
#include "realvector.h"

class Integrator {
public:
	Integrator(std::vector<std::string> arg);
	virtual ~Integrator();

	std::vector<std::string> system_configuration;
	System system;

	std::vector<System> snapshots;

	const ButcherTable runge_kutta_method = ssp_rk4;

	PhysicalDouble kappa_min = 0;
	PhysicalDouble precision = 1e-4;

	void restart_system(PhysicalDouble kappa = -1);
	int integrate();
	void find_criticality();
	void save_snapshots(std::string file);

	const GaussQuadrature GLIntegrator;

};

#endif /* INTEGRATOR_H_ */
