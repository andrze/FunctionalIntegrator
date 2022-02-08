//#define _GNU_SOURCE
#include <fenv.h>
#include <iostream>
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
		std::cerr << "Nie podano argumentów wejściowych\n";
		return 0;
	}

	auto start = std::chrono::steady_clock::now();
	std::vector<std::string> arg(argv, argv + argc);

	std::string outfile = "outfile.csv";
	size_t num_threads = 1;
	for (size_t i = 0; i < arg.size() - 1; i++) {
		std::string opt = arg[i];
		if (opt == "-out") {
			outfile = arg[i + 1];
		} else if (opt == "-threads") {
			num_threads = size_t(std::atoi(arg[i + 1].c_str()));
		}
	}
	std::cout << "Running calculations on "<< num_threads << " threads.\n";

	Integrator integrator(arg, num_threads);

	std::string task = std::string(argv[1]); // First command line argument reserved for task specification
	if (task == "single") {
		std::cout << integrator.system.print_configuration();
		integrator.integrate();
	} else if (task == "critical") {
		integrator.find_criticality();
	} else {
		std::cerr << "Podano niewłaściwe zadanie dla programu \n";
		return 0;
	}

	integrator.save_snapshots(outfile);

	auto end = std::chrono::steady_clock::now();
	std::cout << "Calculations lasted for " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
			<< " seconds.\n";

	return 0;
}
