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
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <queue>
#include "system.h"
#include "stepfunction.h"
#include "gaussquadrature.h"
#include "rungekutta.h"
#include "realvector.h"

struct Task {
	Task(Integrator *integrator, std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> *integrand,
		PhysicalDouble *result);
	class Integrator *integrator;
	std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> *integrand;
	PhysicalDouble *result;
	bool shutdown = false;
};

class Integrator {
public:
	Integrator(std::vector<std::string> arg, size_t num_threads = 1);
	virtual ~Integrator();

	std::vector<std::string> system_configuration;
	System system;

	std::vector<System> snapshots;

	const ButcherTable runge_kutta_method = rk4; //ssp_rk4;

	PhysicalDouble kappa_min = 0, kappa_max;
	PhysicalDouble precision = 1e-4;
	PhysicalDouble max_time = 20;

	void restart_system(PhysicalDouble kappa = -1);
	int integrate(bool verbose = false);
	void find_criticality();
	void save_snapshots(std::string file);

	size_t num_threads = 24; //number of threads on which integration is evaluated
	std::queue<Task*> tasks; //queue of tasks to be executed by threads
	std::atomic<int> active_tasks_count { 0 }; //number of currently active tasks
	std::atomic<size_t> current_task { 0 }; //number of currently active tasks
	std::mutex tasks_mutex; //mutex locking tasks queue and tasksActive counter
	std::condition_variable all_tasks_done; //notified when tasksActive is set to 0
	std::condition_variable tasks_not_empty; //notified when task is being pushed to tasks queue

	GaussQuadrature GLIntegrator;

	void push_integrand_function(std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> f);
	void reset_integrals(size_t new_size = 0);
	std::vector<PhysicalDouble>* evaluate_integrals();

private:
	std::vector<std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> > integrand_functions; //vector of integrals to be carried out
	std::vector<PhysicalDouble> integral_values; //vector of results of the integrals
	std::vector<std::thread> threads; //vector of threads working on assigned tasks
	std::vector<std::unique_ptr<std::thread> > workers;

};

void integrate(Task *task, Integrator *integrator);
void executeTasks(Integrator *integrator);

#endif /* INTEGRATOR_H_ */
