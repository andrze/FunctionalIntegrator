#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H
#include <vector>
#include <functional>
#include "stepfunction.h"
#include "plot.h"
#include <limits>

class System {
public:
	System();
	System(StepFunction V, double delta_t, double d, double n, size_t norm_point, bool sigma_normalization);
	System(std::vector<std::string> configuration, double kappa = -1);

	std::vector<StepFunction> parameters;
	int phase = 0;

	double time = 0.;
	double delta_t = 1e-3;
	double min_delta_t = 1e-6;
	size_t step = 0;
	double time_after_phase_diagnosis = 2.;

	double eta = 0.;
	double z_dim = 1.;

	double d = 3.;
	double vd = 1.;
	double n = 1.;
	size_t norm_point = 0;
	bool sigma_normalization = true;
	size_t num_points = 60;
	double rho0_to_rhomax = 1.5;
	double kappa = 1.;
	double u = 0.1;
	double a = 2.;

	void reparametrize();
	size_t num_functions();
	StepFunction& operator[](size_t i);
	StepFunction& V();
	StepFunction& Zs();
	StepFunction& Zp();

	void find_eta();
	System time_derivative(std::vector<Plot> *integrands = nullptr);

	void rk4_step(StepFunction *y_der = nullptr); // RK4
	void ssp_rk3_step(StepFunction *y_der = nullptr); // Strong Stability Preserving RK3
	void euler_step(StepFunction *y_der = nullptr); // Euler method

	void rescale();
	void fix_v();

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
	double prefactor(double y, double eta);
	double G(double m, double Z, double y);
};

System operator+(System lhs, System rhs);
System operator*(double lhs, System rhs);
double time_eta_distance(System lhs, System rhs);

std::ostream& operator<<(std::ostream &out, System s);

#endif // PARAMETRIZATION_H
