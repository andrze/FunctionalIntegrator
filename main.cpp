//#define _GNU_SOURCE
#include <fenv.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include "realvector.h"
#include "stepfunction.h"
#include "gaussquadrature.h"
#include "regulator.h"
#include "system.h"
#include "integrator.h"

int main(int argc, char *argv[]) {
	//feenableexcept(FE_INVALID | FE_OVERFLOW);

	feraiseexcept(FE_INVALID | FE_OVERFLOW);

	if (argc < 2) {
		std::cerr << "Program requires a supported task name as the first argument\n";
		return 0;
	}

	auto start = std::chrono::steady_clock::now();
	std::vector<std::string> arg(argv, argv + argc);

	std::string outfile = "outfile.csv";
	size_t num_threads = 1;
	bool verbose = false;
	for (size_t i = 0; i < arg.size(); i++) {
		std::string opt = arg[i];
		if (opt == "-out" && i < arg.size() - 1) {
			outfile = arg[i + 1];
		} else if (opt == "-threads" && i < arg.size() - 1) {
			num_threads = size_t(std::atoi(arg[i + 1].c_str()));
		} else if (opt == "verbose") {
			verbose = true;
		}

	}
	std::cout << "Running calculations on " << num_threads << " threads.\n";

	Integrator integrator(arg, num_threads);

	std::string task = std::string(argv[1]); // First command line argument reserved for task specification
	if (task == "single") {
		std::cout << "Running single flow integration\n";
		std::cout << "Flow configuration " << integrator.system.print_configuration() << '\n';
		integrator.integrate(verbose);
		auto snaps = integrator.snapshots;
		//std::cout << std::setprecision(10);
		//std::cout << snaps[snaps.size() - 6].time_derivative() << '\n';
	} else if (task == "critical") {
		std::cout << "Searching for the critical temperature.\n";
		integrator.find_criticality();
	} else {
		std::cerr << "Program requires a supported task name as the first argument\n";
		return 0;
	}

	integrator.save_snapshots(outfile);

	auto end = std::chrono::steady_clock::now();
	int seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
	std::cout << "Calculations lasted for " << seconds / 3600 << ":" << (seconds % 3600) / 60 << ":" << (seconds % 60)
		<< ".\n";

	return 0;
}
