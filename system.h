#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H
#include <vector>
#include <functional>
#include <limits>
#include "stepfunction.h"
#include "realvector.h"

class Integrator;

class SystemDelta {
public:
	SystemDelta();
	SystemDelta(std::vector<StepFunction> parameters);
	size_t num_functions();
	RealVector vector_representation();

	StepFunction& operator[](size_t i);
	SystemDelta& operator+=(SystemDelta rhs);
	SystemDelta& operator*=(PhysicalDouble rhs);

	void plot_parameters();

private:
	std::vector<StepFunction> parameters;
};

SystemDelta operator+(SystemDelta lhs, SystemDelta rhs);
SystemDelta operator-(SystemDelta lhs, PhysicalDouble rhs);
SystemDelta operator*(PhysicalDouble lhs, SystemDelta rhs);
SystemDelta operator/(SystemDelta lhs, SystemDelta rhs);

class System {
public:
	System();
	System(Integrator *integrator, StepFunction V, PhysicalDouble delta_t, PhysicalDouble d, PhysicalDouble n,
		bool sigma_normalization);
	System(Integrator *integrator, std::vector<std::string> configuration, PhysicalDouble kappa = -1);
	System(Integrator *integrator, std::string filename);

	Integrator *integrator = nullptr;
	int phase = 0;

	PhysicalDouble time = 0.;
	PhysicalDouble delta_t = 1e-3;
	size_t step = 0;
	PhysicalDouble time_after_phase_diagnosis = 0.5;

	PhysicalDouble eta = 0.;
	PhysicalDouble z_dim = 1.;
	PhysicalDouble beta_squared = std::numeric_limits<PhysicalDouble>::max();

	PhysicalDouble d = 3., d_inv = 1. / 3;
	PhysicalDouble vd = 1.;
	PhysicalDouble n = 2;
	bool sigma_normalization = true;
	size_t norm_point = 0;
	size_t num_points = 60;
	PhysicalDouble rho0_to_rhomax = 1.5;
	PhysicalDouble kappa = 1.;
	PhysicalDouble u = 0.1;
	PhysicalDouble a = 2.;
	PhysicalDouble step_size = 0;

	PhysicalDouble z_correction = 1.;
	bool zoomed = false;

	void reparametrize();
	size_t num_functions();
	StepFunction operator[](size_t i);
	StepFunction V();
	StepFunction Zs();
	StepFunction Zp();
	RealVector full_vector_representation();
	StepFunction rho_func, V1, V2, Zs1, Zs2, Zp1, Zp2;
	void precalculate_rho_derivatives();
	bool rho_derivatives_calculated = false;

	void find_eta();
	void push_time_derivative_integrals(size_t i);
	PhysicalDouble time_derivative(size_t func_i, size_t rho_i, size_t integrals_position);
	SystemDelta time_derivative();
	SystemDelta beta_function, scaling_coefficients;

	void rescale();
	void zoom_in();
	void cut_domain();

	PhysicalDouble last_val();
	PhysicalDouble first_val();
	int find_phase();
	std::string print_phase();

	std::array<PhysicalDouble, 3> kappa_u_z();
	std::string print_configuration();

	PhysicalDouble y_step = 0.05, y_max = 4;
	void cache_regulator();
	std::vector<PhysicalDouble> cached_regulator_vals;

	void plot_parameters();

	System& operator+=(SystemDelta rhs);

private:
	std::vector<StepFunction> parameters;
};

PhysicalDouble time_eta_distance(System lhs, System rhs);

std::ostream& operator<<(std::ostream &out, System s);

#endif // PARAMETRIZATION_H
