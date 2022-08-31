/*
 * Integrator.cpp
 *
 *  Created on: Sep 9, 2021
 *      Author: andrzej
 */
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include "integrator.h"
#include "rungekutta.h"
#include "terminalplot.h"

void integrate(Task *task, Integrator *integrator) {
	*(task->result) = integrator->GLIntegrator.integrate(*(task->integrand));
}

void executeTasks(Integrator *integrator) {
	while (true) {
		std::unique_lock<std::mutex> lock(integrator->tasks_mutex);
		while (integrator->tasks.empty()) {
			integrator->tasks_not_empty.wait(lock);
		}
		Task *task = integrator->tasks.front();
		integrator->tasks.pop();
		lock.unlock();
		if (task->shutdown) {
			return;
		}
		integrate(task, integrator);
		lock.lock();
		integrator->active_tasks_count--;
		if (integrator->active_tasks_count == 0) {
			integrator->all_tasks_done.notify_one();
		}
		lock.unlock();
	}
}

Task::Task(Integrator *integrator, std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> *integrand,
		PhysicalDouble *result) :
		integrator(integrator), integrand(integrand), result(result) {
}

Integrator::Integrator(std::vector<std::string> arg, size_t num_threads) :
		system_configuration(arg), num_threads(num_threads), tasks_not_empty() {
	if (num_threads <= 0) {
		throw std::runtime_error("Number of threads should be larger than 0.");
	}

	system = System(this, arg);
	kappa_max = system.kappa;

	for (size_t i = 0; i < arg.size() - 1; i++) {
		std::string opt = arg[i];
		if (opt == "-kappa_min") {
			kappa_min = std::strtod(arg[i + 1].c_str(), nullptr);
		} else if (opt == "-precision") {
			precision = std::strtod(arg[i + 1].c_str(), nullptr);
		}
	}

	for (size_t i = 0; i < this->num_threads; i++) {
		threads.emplace_back(executeTasks, this);
	}

	GLIntegrator = GaussQuadrature(system.d);
}

//Integrator::Integrator(std::vector<std::string> arg, size_t num_threads) :
//	system_configuration(arg), num_threads(num_threads) {
//	if (num_threads <= 0) {
//		throw std::runtime_error("Number of threads should be larger than 0.");
//	}
//
//	system = System(this, arg);
//
//	for (size_t i = 0; i < arg.size() - 1; i++) {
//		std::string opt = arg[i];
//		if (opt == "-kappa_min") {
//			kappa_min = std::strtod(arg[i + 1].c_str(), nullptr);
//		} else if (opt == "-precision") {
//			precision = std::strtod(arg[i + 1].c_str(), nullptr);
//		}
//	}
//
//	GLIntegrator = GaussQuadrature(system.d);
//}

Integrator::~Integrator() {
	Task task(this, nullptr, nullptr);
	task.shutdown = true;
	std::unique_lock<std::mutex> lock(tasks_mutex);
	for (size_t t = 0; t < num_threads; t++) {
		tasks.push(&task);
		tasks_not_empty.notify_one();
	}
	lock.unlock();
	for (auto &&t : threads) {
		t.join();
	}
}

void Integrator::restart_system(PhysicalDouble kappa) {
	snapshots.clear();
	system = System(this, system_configuration, kappa);
}

int Integrator::integrate(bool verbose) {
	size_t max_steps = 1e+7;
	PhysicalDouble max_time = 50;
	system.time_derivative();
	snapshots.push_back(system);
	TerminalPlot plot;
//	PhysicalDouble base_delta_t = system.delta_t;

	for (size_t i = 0; i < max_steps && system.time < max_time; i++) {
		if (i % 1000 == 0) {
			std::cout << std::setprecision(5) << system.time << ", " << system.eta << '\n';
			if (verbose) {
				plot.plot(system.parameters);
			}
		}
		try {
			runge_kutta_method.runge_kutta_step(system);
//			system.cut_domain();
			//system.zoom_in();
			//system.rescale();

		} catch (const std::runtime_error &err) {
			std::cout << err.what() << "\n";

			std::cout << system.print_phase();

			return system.find_phase();
		}
		int phase = system.find_phase();
		if (phase != 0) {
			system.time_after_phase_diagnosis -= system.delta_t;
			if (system.time_after_phase_diagnosis <= 0) {
				std::cout << '\n';
				std::cout << system.print_phase();
				return phase;
			}
		}

		if (snapshots.size() == 0 || time_eta_distance(system, snapshots[snapshots.size() - 1]) > 0.05) {
			snapshots.push_back(system);
		}
//
//		PhysicalDouble distance_from_pole = system.a + system.V()[0];
//		PhysicalDouble new_delta = base_delta_t * std::pow(10, int(std::log10(distance_from_pole)));
//		if(std::abs(new_delta - system.delta_t) > 1e-15){
//			if(new_delta < 1e-9){
//				break;
//			}
//			std::cout << std::setprecision(16);
//			std::cout << system.V()<<'\n';
//			std::cout << distance_from_pole << ' ' << int(-std::log10(distance_from_pole)) << ' ' << std::pow(10, 1+int(-std::log10(distance_from_pole))) <<  '\n';
//			std::cout << "Changing delta_t to " << new_delta << '\n';
//			system.delta_t = new_delta;
//		}

	}
	std::cout << '\n';
	std::cout << "Przekroczono maksymalną liczbę kroków.\n";
	return system.find_phase();
}

void Integrator::find_criticality() {
	kappa_max = system.kappa;

	size_t max_iter = 50;

	for (size_t i = 0; i < max_iter && kappa_max - kappa_min > precision; i++) {
		restart_system((kappa_min + kappa_max) / 2);
		std::cout << std::setprecision(-int(std::log10(precision))+2);
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

	out << "t,eta,Z,beta\n";
	for (size_t i = 0; i < snapshots.size(); i++) {
		System s = snapshots[i];
		out << s.time << ',' << s.eta << ',' << s.z_dim << ',' << s.beta_squared << '\n';
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

void Integrator::push_integrand_function(std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> f) {
	integrand_functions.push_back(f);
}

void Integrator::reset_integrals(size_t new_size) {
	integrand_functions.clear();
	if (new_size > 0) {
		integrand_functions.reserve(new_size);
	}
}

std::vector<PhysicalDouble>* Integrator::evaluate_integrals() {
	size_t task_count = integrand_functions.size();
	integral_values = std::vector<PhysicalDouble>(task_count);
//
//	for (size_t i = 0; i < task_count; i++) {
//		integral_values[i] = GLIntegrator.integrate(integrand_functions[i]);
//	}

	std::vector<Task> task_vector;
	task_vector.reserve(task_count);
	std::unique_lock<std::mutex> lock(tasks_mutex);

	for (size_t i = 0; i < task_count; i++) {
		task_vector.emplace_back(this, &(integrand_functions[i]), &(integral_values[i]));
		tasks.push(&task_vector.back());
		active_tasks_count++;
		tasks_not_empty.notify_one();
	}

	while (active_tasks_count != 0) {
		all_tasks_done.wait(lock);
	}

	return &integral_values;
}
