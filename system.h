#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H
#include <vector>
#include <functional>
#include "stepfunction.h"
#include "plot.h"
#include <limits>

class Integrator;

class System {
public:
	System();
	System(Integrator* integrator, StepFunction V, double delta_t, double d, double n, bool sigma_normalization);
	System(Integrator* integrator, std::vector<std::string> configuration, double kappa = -1);

	Integrator* integrator=nullptr;
	std::vector<StepFunction> parameters;
	int phase = 0;

	double time = 0.;
	double delta_t = 1e-3;
	size_t step = 0;
	double time_after_phase_diagnosis = 2.;

	double eta = 0.;
	double z_dim = 1.;

	double d = 3.;
	double vd = 1.;
	double n = 1.;
	bool sigma_normalization = true;
	size_t num_points = 60;
	double rho0_to_rhomax = 1.5;
	double kappa = 1.;
	double u = 0.1;
	double a = 2.;

	double z_violation = 1.;
	bool zoomed = false;

	void reparametrize();
	size_t num_functions();
	StepFunction& operator[](size_t i);
	StepFunction& V();
	StepFunction& Zs();
	StepFunction& Zp();

	void find_eta();
	System time_derivative();

	void rescale();
	void zoom_in();

	double last_val();
	double first_val();
	int find_phase();
	std::string print_phase(int phase = 0);

	std::array<double, 3> kappa_u_z();
	std::string print_configuration();

	System& operator+=(System rhs);
	System& operator*=(double rhs);

	double r(double y);
	double rp(double y);
	double rp2(double y);
	double prefactor(double y);
	double G(double m, double Z, double y);

	double gauss_legendre_integrate(std::function<double(double)> f);
};

System operator+(System lhs, System rhs);
System operator*(double lhs, System rhs);
double time_eta_distance(System lhs, System rhs);

std::ostream& operator<<(std::ostream &out, System s);

#endif // PARAMETRIZATION_H
