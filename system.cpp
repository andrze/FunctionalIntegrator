#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <sstream>
#include "system.h"
#include "regulator.h"
#include "gaussquadrature.h"

double Power(double x, double y) {
	return std::pow(x, y);
}

System::System() {
}

System::System(StepFunction V, double delta_t, double d, double n, bool sigma_normalization) :
		delta_t(delta_t), d(d), n(n), sigma_normalization(sigma_normalization) {
	parameters.push_back(V);
	parameters.push_back(StepFunction(V.step_size, V.num_points, [](double) {
		return 1;
	}));
	parameters.push_back(StepFunction(V.step_size, V.num_points, [](double) {
		return 1;
	}));
	vd = std::pow(2., -1 - d) * std::pow(M_PI, -d / 2) / std::tgamma(d / 2);
}

System::System(std::vector<std::string> configuration, double kappa) {

	for (size_t i = 0; i < configuration.size() - 1; i++) {
		std::string opt = configuration[i], val = configuration[i + 1];
		if (opt == "-kappa") {
			this->kappa = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-u") {
			u = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-delta") {
			delta_t = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-min_delta") {
			min_delta_t = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-dim") {
			d = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-N") {
			n = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-num_points") {
			num_points = size_t(std::atoi(val.c_str()));
		} else if (opt == "-rhomax") {
			rho0_to_rhomax = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-a") {
			a = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-sigma_normalization") {
			sigma_normalization = (val == "true");
		}
	}

	if (kappa > 0) {
		this->kappa = kappa;
	}

	reparametrize();

	vd = std::pow(2., -1 - d) * std::pow(M_PI, -d / 2) / std::tgamma(d / 2);

}

void System::reparametrize() {

	double step_size = kappa * rho0_to_rhomax / num_points;

	StepFunction V(step_size, num_points, [&](double x) {
		return u * (x - kappa);
	});

	parameters.clear();
	parameters.push_back(V);
	parameters.push_back(StepFunction(V.step_size, V.num_points, [=](double) {
		return 1;
	}));
	parameters.push_back(StepFunction(V.step_size, V.num_points, [=](double) {
		return 1;
	}));
}

size_t System::num_functions() {
	return parameters.size();
}

StepFunction& System::operator[](size_t i) {
	return parameters[i];
}

StepFunction& System::V() {
	return parameters[0];
}

StepFunction& System::Zs() {
	return parameters[1];
}

StepFunction& System::Zp() {
	return parameters[2];
}

void System::find_eta() {
	double potential_minimum = V().kappa_u().first;
	if (potential_minimum <= std::numeric_limits<double>::epsilon()) {
		throw std::runtime_error("Potential minimum is out of the grid");
	}
	size_t norm = size_t((potential_minimum - V().domain_begin) / V().step_size);

	double v1 = V()(potential_minimum);
	double v2 = V().derivative(1)(potential_minimum);
	double v3 = V().derivative(2)(potential_minimum);
	double zs = Zs()(potential_minimum);
	double zs1 = Zs().derivative(1)(potential_minimum);
	double zs2 = Zs().derivative(2)(potential_minimum);
	double zp = Zp()(potential_minimum);
	double zp1 = Zp().derivative(1)(potential_minimum);
	double zp2 = Zp().derivative(2)(potential_minimum);
	double rho = potential_minimum;

	double z, z1;
	if (sigma_normalization) {
		z = zs;
		z1 = zs1;
	} else {
		z = zp;
		z1 = zp1;
	}

	auto v_common_part = [=](double y2) {
		double gp = std::pow(v1 + y2 * zp + r(y2), -1);
		double gs = std::pow(v1 + 2 * rho * v2 + y2 * zs + r(y2), -1);
		double res = -2 * Power(gp, 2) * (-1 + n) * (v2 + y2 * zp1)
				- 2 * Power(gs, 2) * (3 * v2 + 2 * rho * v3 + y2 * zs1);
		return res;
	};

	auto v_free_integrand = [=](double y) {
		double y2 = y * y;
		return prefactor(y2) *  v_common_part(y2) * std::pow(y, d - 1.);
	};

	auto v_eta_integrand = [=](double y) {
		double y2 = y * y;
		return r(y2) * v_common_part(y2) * std::pow(y, d - 1.);
	};

	auto z_common_part =
			[=](double y2) {
				if (norm == 0) {
					double g = std::pow(v1 + y2 * zp + r(y2), -1);

					return -2 * Power(g, 2) * ((-1 + n) * zp1 + zs1);
				} else if (sigma_normalization) {

					double gp = std::pow(v1 + y2 * zp + r(y2), -1);
					double gs = std::pow(v1 + 2 * rho * v2 + y2 * zs + r(y2), -1);

					return 8 * Power(gp, 3) * (-1 + n) * ((rho * y2 * Power(zp1, 2)) / d + (v2 + y2 * zp1) * (-zp + zs))
							- (2 * Power(gp, 2) * (-1 + n) * (zp - zs + rho * zs1)) / rho
							+ (8 * Power(gs, 3) * rho * zs1 * (y2 * zs1 + 2 * d * (3 * v2 + 2 * rho * v3 + y2 * zs1)))
									/ d - 2 * Power(gs, 2) * (zs1 + 2 * rho * zs2)
							+ (32 * Power(gs, 5) * rho * y2 * Power(3 * v2 + 2 * rho * v3 + y2 * zs1, 2)
									* Power(zs + rp(y2), 2)) / d
							- (8 * Power(gp, 5) * (-1 + n) * rho * Power(v2 + y2 * zp1, 2) * (zp + rp(y2))
									* (d * v1 + (-4 + d) * y2 * zp + d * r(y2) - 4 * y2 * rp(y2))) / d
							- (16 * Power(gp, 4) * (-1 + n) * rho * y2 * (v2 + y2 * zp1)
									* (2 * zp1 * (zp + rp(y2)) + (v2 + y2 * zp1) * rp2(y2))) / d
							- (8 * Power(gs, 4) * rho * (3 * v2 + 2 * rho * v3 + y2 * zs1)
									* ((3 * d * v2 + 2 * d * rho * v3 + (4 + d) * y2 * zs1) * (zs + rp(y2))
											+ 2 * y2 * (3 * v2 + 2 * rho * v3 + y2 * zs1) * rp2(y2))) / d;
				} else {

					double gp = std::pow(v1 + y2 * zp + r(y2), -1);
					double gs = std::pow(v1 + 2 * rho * v2 + y2 * zs + r(y2), -1);

					return (2
							* (d * Power(gp, 2) * (zp + rho * (zp1 - n * zp1) - zs)
									- 2 * Power(gp, 2) * gs
											* (-2 * y2 * Power(zp + rho * zp1 - zs, 2)
													+ d * (zp - zs) * (2 * rho * v2 + y2 * (-zp + zs)))
									+ 4 * Power(gp, 2) * Power(gs, 3) * y2 * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
											* Power(zs + rp(y2), 2)
									- Power(gs, 2)
											* (d * rho * (zp1 + 2 * rho * zp2)
													- 4 * gp * rho * zp1
															* (2 * d * rho * v2 + rho * y2 * zp1 + d * y2 * (-zp + zs))
													+ Power(gp, 3) * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
															* (zp + rp(y2))
															* (d * v1 + (-4 + d) * y2 * zp + d * r(y2) - 4 * y2 * rp(y2))
													+ Power(gp, 2)
															* (d * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
																	* (zs + rp(y2))
																	+ 8 * y2 * (-2 * rho * v2 + y2 * (zp - zs))
																			* (zp - zs) * (-(rho * zp1) + zs + rp(y2))
																	+ 4 * y2 * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
																			* rp2(y2))))) / (d * rho);
				}

			};

	auto z_free_integrand = [=](double y) {
		double y2 = y * y;
		return prefactor(y2) * z_common_part(y2) * std::pow(y, d - 1.);
	};

	auto z_eta_integrand = [=](double y) {
		double y2 = y * y;
		return r(y2) * z_common_part(y2) * std::pow(y, d - 1.);
	};

	double v_free = gauss_legendre_integrate(v_free_integrand);
	double v_eta = gauss_legendre_integrate(v_eta_integrand);
	double z_free = gauss_legendre_integrate(z_free_integrand);
	double z_eta = gauss_legendre_integrate(z_eta_integrand);

	double free_component, eta_component;

	free_component = vd * (z_free * v2 - v_free * z1);
	eta_component = -v2 * z + vd * (v_eta * z1 - z_eta * v2);

	eta = free_component / eta_component;
}

System System::time_derivative(std::vector<Plot> *integrands) {
	find_eta();
	StepFunction V = this->V();
	StepFunction V2 = V.derivative(1);
	StepFunction V3 = V.derivative(2);
	StepFunction Zs = this->Zs();
	StepFunction Zs1 = Zs.derivative(1);
	StepFunction Zs2 = Zs.derivative(2);
	StepFunction Zp = this->Zp();
	StepFunction Zp1 = Zp.derivative(1);
	StepFunction Zp2 = Zp.derivative(2);

	size_t num_points = V.num_points;
	double step = V.step_size;
	std::vector<double> v_vals, zs_vals, zp_vals;
	StepFunction rho_func = V.x_func();

	for (size_t i = 0; i < num_points; i++) {
		double rho = rho_func[i];
		double v1 = V[i], v2 = V2[i], v3 = V3[i];
		double zs = Zs[i], zs1 = Zs1[i], zs2 = Zs2[i];
		double zp = Zp[i], zp1 = Zp1[i], zp2 = Zp2[i];

		bool singularity = false;
		std::string error;
		if (v1 < -a) {
			error = "V(" + std::to_string(i) + ") = " + std::to_string(v1);
			singularity = true;
		} else if ((v1 + 2 * rho * v2) < -a) {
			error = "V + 2 r V'(" + std::to_string(i) + ") = " + std::to_string((v1 + 2 * rho * v2));
			singularity = true;
		} else if (zs < 0) {
			error = "Zs(" + std::to_string(i) + ") = " + std::to_string(zs);
			singularity = true;
		} else if (zp < 0) {
			error = "Zp(" + std::to_string(i) + ") = " + std::to_string(zp);
			singularity = true;
		}

		if (singularity) {
			throw std::runtime_error("Equations contain a singular propagator " + error);
		}

		// V flow equation
		auto v_integrand = [=](double y) {
			double y2 = y * y;
			double gp = std::pow(v1 + y2 * zp + r(y2), -1);
			double gs = std::pow(v1 + 2 * rho * v2 + y2 * zs + r(y2), -1);
			double pref = prefactor(y2, eta);
			double res = -2 * Power(gp, 2) * (-1 + n) * (v2 + y2 * zp1)
					- 2 * Power(gs, 2) * (3 * v2 + 2 * rho * v3 + y2 * zs1);
			return std::pow(y, d - 1) * res * pref * std::pow(y, d - 1);
		};
		v_vals.push_back(gauss_legendre_integrate(v_integrand));

		// Zs flow equation
		if (rho <= std::numeric_limits<double>::epsilon()) {
			auto zs_integrand = [=](double y) {
				double y2 = y * y;
				double gp = std::pow(v1 + y2 * zp + r(y2), -1);
				double pref = prefactor(y2, eta);

				double res = -2 * Power(gp, 2) * ((-1 + n) * zp1 + zs1);
				return std::pow(y, d - 1) * res * pref;
			};
			zs_vals.push_back(gauss_legendre_integrate(zs_integrand));
		} else {

			auto zs_integrand = [=](double y) {
				double y2 = y * y;
				double gp = std::pow(v1 + y2 * zp + r(y2), -1);
				double gs = std::pow(v1 + 2 * rho * v2 + y2 * zs + r(y2), -1);
				double pref = prefactor(y2, eta);

				double res = 8 * Power(gp, 3) * (-1 + n)
						* ((rho * y2 * Power(zp1, 2)) / d + (v2 + y2 * zp1) * (-zp + zs))
						- (2 * Power(gp, 2) * (-1 + n) * (zp - zs + rho * zs1)) / rho
						+ (8 * Power(gs, 3) * rho * zs1 * (y2 * zs1 + 2 * d * (3 * v2 + 2 * rho * v3 + y2 * zs1))) / d
						- 2 * Power(gs, 2) * (zs1 + 2 * rho * zs2)
						+ (32 * Power(gs, 5) * rho * y2 * Power(3 * v2 + 2 * rho * v3 + y2 * zs1, 2)
								* Power(zs + rp(y2), 2)) / d
						- (8 * Power(gp, 5) * (-1 + n) * rho * Power(v2 + y2 * zp1, 2) * (zp + rp(y2))
								* (d * v1 + (-4 + d) * y2 * zp + d * r(y2) - 4 * y2 * rp(y2))) / d
						- (16 * Power(gp, 4) * (-1 + n) * rho * y2 * (v2 + y2 * zp1)
								* (2 * zp1 * (zp + rp(y2)) + (v2 + y2 * zp1) * rp2(y2))) / d
						- (8 * Power(gs, 4) * rho * (3 * v2 + 2 * rho * v3 + y2 * zs1)
								* ((3 * d * v2 + 2 * d * rho * v3 + (4 + d) * y2 * zs1) * (zs + rp(y2))
										+ 2 * y2 * (3 * v2 + 2 * rho * v3 + y2 * zs1) * rp2(y2))) / d;
				return std::pow(y, d - 1) * res * pref;
			};

			zs_vals.push_back(gauss_legendre_integrate(zs_integrand));
		}

		//Zp flow equation
		if (rho <= std::numeric_limits<double>::epsilon()) {
			auto zp_integrand = [=](double y) {
				double y2 = y * y;
				double gp = std::pow(v1 + y2 * zp + r(y2), -1);
				double pref = prefactor(y2, eta);

				double res = -2 * Power(gp, 2) * ((-1 + n) * zp1 + zs1);
				return std::pow(y, d - 1) * res * pref;
			};

			zp_vals.push_back(gauss_legendre_integrate(zp_integrand));

		} else {
			auto zp_integrand = [=](double y) {
				double y2 = y * y;
				double gp = std::pow(v1 + y2 * zp + r(y2), -1);
				double gs = std::pow(v1 + 2 * rho * v2 + y2 * zs + r(y2), -1);
				double pref = prefactor(y2, eta);

				double res = (2
						* (d * Power(gp, 2) * (zp + rho * (zp1 - n * zp1) - zs)
								- 2 * Power(gp, 2) * gs
										* (-2 * y2 * Power(zp + rho * zp1 - zs, 2)
												+ d * (zp - zs) * (2 * rho * v2 + y2 * (-zp + zs)))
								+ 4 * Power(gp, 2) * Power(gs, 3) * y2 * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
										* Power(zs + rp(y2), 2)
								- Power(gs, 2)
										* (d * rho * (zp1 + 2 * rho * zp2)
												- 4 * gp * rho * zp1
														* (2 * d * rho * v2 + rho * y2 * zp1 + d * y2 * (-zp + zs))
												+ Power(gp, 3) * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
														* (zp + rp(y2))
														* (d * v1 + (-4 + d) * y2 * zp + d * r(y2) - 4 * y2 * rp(y2))
												+ Power(gp, 2)
														* (d * Power(2 * rho * v2 + y2 * (-zp + zs), 2) * (zs + rp(y2))
																+ 8 * y2 * (-2 * rho * v2 + y2 * (zp - zs)) * (zp - zs)
																		* (-(rho * zp1) + zs + rp(y2))
																+ 4 * y2 * Power(2 * rho * v2 + y2 * (-zp + zs), 2)
																		* rp2(y2))))) / (d * rho);
				return std::pow(y, d - 1) * res * pref;
			};
			zp_vals.push_back(gauss_legendre_integrate(zp_integrand));
		}
	}

	StepFunction v_der = (2 - eta) * V - (d - 2 + eta) * rho_func * V2
			- vd * StepFunction(step, v_vals, V.domain_begin);
	StepFunction zs_der = -eta * Zs - (d - 2 + eta) * rho_func * Zs1 - vd * StepFunction(step, zs_vals, V.domain_begin);
	StepFunction zp_der = -eta * Zp - (d - 2 + eta) * rho_func * Zp1 - vd * StepFunction(step, zp_vals, V.domain_begin);

	System derivative = *this;
	derivative.parameters = std::vector<StepFunction>();
	derivative.parameters.push_back(v_der);
	derivative.parameters.push_back(zs_der);
	derivative.parameters.push_back(zp_der);
	derivative.eta = eta;

	return derivative;
}

void System::rescale() {
	double potential_minimum = V().kappa_u().first;
	if (potential_minimum <= std::numeric_limits<double>::epsilon()) {
		throw std::runtime_error("Potential minimum is out of the grid");
	}
	double z_norm;
	if (sigma_normalization) {
		z_norm = Zs()(potential_minimum);
	} else {
		z_norm = Zp()(potential_minimum);
	}

	for (auto &&func : this->parameters) {
		func = func / z_norm;
		func.step_size = func.step_size / z_norm;
		func.domain_begin = func.domain_begin / z_norm;
	}
	if (z_norm > 1.01 || z_norm < .99) {
		std::cout << "Rescaling by a factor of " << z_norm << ", old potential minimum at " << potential_minimum
				<< '\n';
	}

	z_dim = z_dim * z_norm;
	z_violation *= z_norm;
}

void System::zoom_in() {
	auto v_vals = V().vals, xs = V().xs();
	const int DELTA_MIN_TO_DELTA_MAX = 10;
	const double V_THRESHOLD = -a + .1;
	const double V_GOAL = -a + .15;

	if (v_vals[0] > V_THRESHOLD && v_vals[0] < v_vals[1]) {
		return;
	}

	size_t rescale_begin_index = 0, rescale_end_index = xs.size() - 1;

	for (size_t i = 1; i < v_vals.size(); i++) {
		if (v_vals[i] > V_GOAL && v_vals[i] > v_vals[i - 1]) {
			rescale_begin_index = i;
			break;
		}
	}

	for (size_t i = rescale_begin_index + 1; i < v_vals.size(); i++) {
		if (v_vals[i] > 0) {
			rescale_end_index = std::min(rescale_begin_index + (i - rescale_begin_index) * (DELTA_MIN_TO_DELTA_MAX + 1),
					xs.size() - 1);
			break;
		}
	}

	double rescale_begin = xs[rescale_begin_index];
	double rescale_end = xs[rescale_end_index];
	for (auto &&func : this->parameters) {
		func = func.zoom_in(rescale_begin, rescale_end);
	}
	std::cout << "Zooming in on the interval [" << rescale_begin << ',' << rescale_end << "].\n";
	zoomed = true;
}

double System::last_val() {
	return V().vals.back();
}

double System::first_val() {
	return V().vals.front();
}

int System::find_phase() {
	phase = 0;
	auto mini_maxi = V().minmax();
	if (mini_maxi.first > 0) { // Disordered phase
		phase = -1;
	} else if (last_val() < 0) { // Ordered phase
		phase = 1;
	}
	return phase;
}

std::string System::print_phase(int phase) {
	if (phase == 0) {
		if (this->phase != 0) {
			phase = this->phase;
		} else {
			phase = find_phase();
		}
	}

	std::string message = "Zakończono w fazie ";
	if (phase == 1) {
		message += "uporządkowanej.\n";
	} else if (phase == -1) {
		message += "nieuporządkowanej.\n";
	} else {
		message += "nieokreślonej.\n";
	}

	return message;
}

std::array<double, 3> System::kappa_u_z() {
	auto kappa_u = V().kappa_u();
	double kappa = kappa_u.first, u = kappa_u.second;
	double z = Zs()(kappa);

	return std::array<double, 3> { { kappa, u, z } };
}

std::string System::print_configuration() {
	std::stringstream configuration;
	configuration << "-delta," << delta_t << ',';
	configuration << "-min_delta," << min_delta_t << ',';
	configuration << "-dim," << d << ',';
	configuration << "-N," << n << ',';
	configuration << "-num_points," << num_points << ',';
	configuration << "-a," << a << ',';
	configuration << "-sigma_normalization," << sigma_normalization << ',';
	/*if (precision < 0) {
	 configuration << "-prec," << precision << ',';
	 }*/
	configuration << "-rho_max," << rho0_to_rhomax << '\n';

	return configuration.str();
}

System& System::operator+=(System rhs) {
	if (num_functions() != rhs.num_functions()) {
		throw std::invalid_argument("Parametrizations with different numbers of functions provided");
	}

	for (size_t i = 0; i < parameters.size(); i++) {
		parameters[i] = parameters[i] + rhs[i];
	}
	return *this;
}

System& System::operator*=(double rhs) {
	for (size_t i = 0; i < parameters.size(); i++) {
		parameters[i] = parameters[i] * rhs;
	}
	return *this;
}

double System::r(double y) {
	return a / 2 * R(y);
}

double System::rp(double y) {
	return a / 2 * Rp(y);
}

double System::rp2(double y) {
	return a / 2 * Rp2(y);
}

double System::prefactor(double y) {
	return a / 2 * Prefactor(y);
}

double System::prefactor(double y, double eta) {
	return a / 2 * Prefactor(y, eta);
}

double System::G(double m, double Z, double y) {
	return 1 / (m + Z * y + r(y));
}

System operator+(System lhs, System rhs) {
	lhs += rhs;
	return lhs;
}

System operator*(double lhs, System rhs) {
	rhs *= lhs;
	return rhs;
}

std::ostream& operator<<(std::ostream &out, System s) {
	out << "Rho,V,Zs,Zp\n";
	for (size_t i = 0; i < s.V().num_points; i++) {
		out << s.V().x_func()[i] << ',' << s.V()[i] << ',' << s.Zs()[i] << ',' << s.Zp()[i] << '\n';
	}

	return out;
}

double time_eta_distance(System lhs, System rhs) {
	double t_dist = lhs.time - rhs.time;
	double eta_dist = lhs.eta - rhs.eta;

	return std::sqrt(t_dist * t_dist + eta_dist * eta_dist);
}
