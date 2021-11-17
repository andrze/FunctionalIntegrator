/*
 * Integrator.cpp
 *
 *  Created on: Sep 9, 2021
 *      Author: andrzej
 */
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "integrator.h"
#include "rungekutta.h"

Integrator::Integrator(std::vector<std::string> arg) :
		system_configuration(arg) {

	system = System(this, arg);
}

Integrator::~Integrator() {
}

void Integrator::restart_system(double kappa) {
	snapshots.clear();
	backup.clear();
	system = System(this, system_configuration, kappa);
}

int Integrator::integrate() {
	size_t max_steps = 10e+7;
	double max_time = 100;
	snapshots.push_back(system);

	for (size_t i = 0; i < max_steps && system.time < max_time; i++) {
		try {
			runge_kutta_method.runge_kutta_step(system);
			system.zoom_in();
			system.rescale();

		} catch (const std::runtime_error &err) {
			std::cout << err.what() << "\n";
			std::cout << system << "\n";
			if (backup.empty()) {
				std::cout << "Faza nieokreślona. Brak kopii zapasowej do przywrócenia.\n";
				return 0;
			}

			system = backup.back();

			std::cout << system.print_phase();

			return system.find_phase();
		}
		int phase = system.find_phase();
		if (phase != 0) {
			system.time_after_phase_diagnosis -= system.delta_t;
			if (system.time_after_phase_diagnosis <= 0) {
				std::cout << '\n';
				std::cout << system.print_phase(phase);
				return phase;
			}
		}

		if (snapshots.size() == 0 || time_eta_distance(system, snapshots[snapshots.size() - 1]) > 0.2) {
			snapshots.push_back(system);
		}

		if (system.step % 500 == 0) {
			if (backup.size() == 10) {
				backup.pop_front();
			}
			backup.push_back(system);
			std::cout << system.time << ", " << std::flush;
		}
	}
	std::cout << '\n';
	std::cout << "Przekroczono maksymalną liczbę kroków.\n";
	return system.find_phase();
}

void Integrator::find_criticality() {
	double kappa_max = system.kappa, kappa_min = this->kappa_min;

	size_t max_iter = 50;

	for (size_t i = 0; i < max_iter && kappa_max - kappa_min > precision; i++) {
		restart_system((kappa_min + kappa_max) / 2);
		std::cout << "Kappa = " << system.kappa << '\n';
		int phase = integrate();
		if (phase == 1) {
			kappa_max = system.kappa;
		} else if (phase == -1) {
			kappa_min = system.kappa;
		} else {
			auto V = system.V();
			size_t num_of_negative = 0;
			for (size_t j = 0; j < V.num_points; j++) {
				if (V[j] < 0) {
					num_of_negative++;
				}
			}

			if (system.time < 1 || (num_of_negative < V.num_points / 3 && V[0] > -system.a * 0.7)) {
				kappa_min = system.kappa;
				std::cout << "Zakończono w fazie nieuporządkowanej (prawdopodobnie)\n";
			} else {
				kappa_max = system.kappa;
				std::cout << "Zakończono w fazie algebraicznej (prawdopodobnie)\n";
			}
		}
	}

	restart_system(kappa_max);
	system.delta_t = 5e-5;
	integrate();

}

void Integrator::save_snapshots(std::string file) {
	if (snapshots.size() == 0) {
		std::cout << "Snapshots are empty \n";
		return;
	}

	std::ofstream out(file.c_str());

	out << "t,eta,Z\n";
	for (size_t i = 0; i < snapshots.size(); i++) {
		System s = snapshots[i];
		out << s.time << ',' << s.eta << ',' << s.z_dim << '\n';
	}
	out << "\n";

	for (auto s : snapshots) {
		out << "Rho,V,Zs,Zp\n";
		for (size_t i = 0; i < s.V().num_points; i++) {
			out << s.V().xs()[i] << ',' << s.V()[i] << ',' << s.Zs()[i] << ',' << s.Zp()[i] << '\n';
		}
		out << "\n";
	}
	out.close();
	return;
}
