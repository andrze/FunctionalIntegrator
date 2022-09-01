#ifndef PARAMETRIZATION_H
#define PARAMETRIZATION_H
#include <vector>
#include <functional>
#include <limits>
#include "stepfunction.h"
#include "realvector.h"

template<typename T, size_t N> bool ArrayLowerEqualComp(std::array<T, N> lhs, std::array<T, N> rhs, size_t pos = 0) {
	if (pos == N || N == 0) {
		return true;
	}
	if (lhs[pos] < rhs[pos]) {
		return true;
	} else if (lhs[pos] > rhs[pos]) {
		return false;
	} else {
		return ArrayLowerEqualComp<T, N>(lhs, rhs, pos + 1);
	}
}

class Integrator;

class System {
public:
	System();
	System(Integrator *integrator, StepFunction V, PhysicalDouble delta_t, PhysicalDouble d, PhysicalDouble n,
		bool sigma_normalization);
	System(Integrator *integrator, std::vector<std::string> configuration, PhysicalDouble kappa = -1);

	Integrator *integrator = nullptr;
	std::vector<StepFunction> parameters;
	int phase = 0;

	PhysicalDouble time = 0.;
	PhysicalDouble delta_t = 1e-3;
	size_t step = 0;
	PhysicalDouble time_after_phase_diagnosis = 2.;

	PhysicalDouble eta = 0.;
	PhysicalDouble z_dim = 1.;
	PhysicalDouble beta_squared = std::numeric_limits<PhysicalDouble>::max();

	PhysicalDouble d = 3., d_inv = 1. / 3;
	PhysicalDouble vd = 1.;
	PhysicalDouble n = 1.;
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
	StepFunction& operator[](size_t i);
	StepFunction& V();
	StepFunction& Zs();
	StepFunction& Zp();
	RealVector full_vector_representation();
	StepFunction rho_func, V2, V3, Zs1, Zs2, Zp1, Zp2;
	void precalculate_rho_derivatives();

	void find_eta();
	void find_eta_minimum();
	void find_eta_norm_point();
	System time_derivative();

	void rescale();
	void zoom_in();
	void cut_domain();

	PhysicalDouble last_val();
	PhysicalDouble first_val();
	int find_phase();
	std::string print_phase();

	std::array<PhysicalDouble, 3> kappa_u_z();
	std::string print_configuration();

	System& operator+=(System rhs);
	System& operator*=(PhysicalDouble rhs);

	PhysicalDouble gauss_legendre_integrate(std::function<PhysicalDouble(PhysicalDouble)> f);
	PhysicalDouble gauss_legendre_integrate(std::function<PhysicalDouble(std::array<PhysicalDouble, 6>)> f);
};

System operator+(System lhs, System rhs);
System operator*(PhysicalDouble lhs, System rhs);
PhysicalDouble time_eta_distance(System lhs, System rhs);

std::ostream& operator<<(std::ostream &out, System s);

#endif // PARAMETRIZATION_H
