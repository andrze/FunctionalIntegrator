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

System::System(StepFunction V, double delta_t, double d, double n, size_t norm_point, bool sigma_normalization) :
		delta_t(delta_t), d(d), n(n), norm_point(norm_point), sigma_normalization(sigma_normalization) {
	parameters.push_back(V);
	parameters.push_back(StepFunction(V.step_size, V.num_points, [](double) {
		return 1;
	}));
	parameters.push_back(StepFunction(V.step_size, V.num_points, [](double) {
		return 1;
	}));
	vd = std::pow(2., -1 - d) * std::pow(M_PI, -d / 2) / std::tgamma(d / 2);
	if (this->norm_point >= V.num_points) {
		this->norm_point = V.num_points - 1;
	}
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
		} else if (opt == "-norm") {
			norm_point = size_t(std::atoi(val.c_str()));
		} else if (opt == "-sigma_normalization") {
			sigma_normalization = (val == "true");
		}
	}

	if (kappa > 0) {
		this->kappa = kappa;
	}

	reparametrize();

	vd = std::pow(2., -1 - d) * std::pow(M_PI, -d / 2) / std::tgamma(d / 2);

	if (this->norm_point >= num_points) {
		this->norm_point = num_points - 1;
	}
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
	size_t norm = norm_point;

	double v1 = V()[norm];
	double v2 = V().derivative(1, norm);
	double v3 = V().derivative(2, norm);
	double zs = Zs()[norm];
	double zs1 = Zs().derivative(1, norm);
	double zs2 = Zs().derivative(2, norm);
	double zp = Zp()[norm];
	double zp1 = Zp().derivative(1, norm);
	double zp2 = Zp().derivative(2, norm);
	double rho = V().step_size * double(norm);

	auto common_part =
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

	auto free_integrand = [=](double y) {
		double y2 = y * y;
		return prefactor(y2) * common_part(y2) * std::pow(y, d - 1.);
	};

	auto eta_integrand = [=](double y) {
		double y2 = y * y;
		return r(y2) * common_part(y2) * std::pow(y, d - 1.);
	};

	double i1 = gauss_legendre_integrate(free_integrand);
	double i2 = gauss_legendre_integrate(eta_integrand);

	double free_component, eta_component;
	if (sigma_normalization) {
		free_component = (d - 2) * rho * zs1 + vd * i1;
		eta_component = zs + rho * zs1 - vd * i2;
	} else {
		free_component = (d - 2) * rho * zp1 + vd * i1;
		eta_component = zp + rho * zp1 - vd * i2;
	}

	eta = -free_component / eta_component;

}

System System::time_derivative(std::vector<Plot> *integrands) {
	find_eta();
	StepFunction V = this->V();
	StepFunction Zs = this->Zs();
	StepFunction V2 = V.derivative(1);
	StepFunction V3 = V.derivative(2);
	StepFunction Zs1 = Zs.derivative(1);
	StepFunction Zs2 = Zs.derivative(2);
	StepFunction Zp = this->Zp();
	StepFunction Zp1 = Zp.derivative(1);
	StepFunction Zp2 = Zp.derivative(2);

	size_t num_points = V.num_points;
	double step = V.step_size;
	std::vector<double> v_vals, zs_vals, zp_vals;

	for (size_t i = 0; i < num_points; i++) {
		bool is_stationary = false;
		double rho = step * double(i);
		double v1 = V[i], v2 = V2[i], v3 = V3[i];
		double zs = Zs[i], zs1 = Zs1[i], zs2 = Zs2[i];
		double zp = Zp[i], zp1 = Zp1[i], zp2 = Zp2[i];
		if (v1 < 0 && V[i + 1] > 0) {
			is_stationary = true;
		}

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
		double threshold = 1e+5;
		if (std::abs(v1) > threshold) {
			error = "V(" + std::to_string(i) + ") = " + std::to_string(v1);
			singularity = true;
		} else if (std::abs(zs1) > threshold) {
			error = "Zs'(" + std::to_string(i) + ") = " + std::to_string(zs1);
			singularity = true;
		} else if (std::abs(zp1) > threshold) {
			error = "Zp'(" + std::to_string(i) + ") = " + std::to_string(zp1);
			singularity = true;
		}
		if (singularity) {
			throw std::runtime_error("System parametrization is singular " + error);
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

		// Plots of v integrands
		if (integrands && (i == 0 || is_stationary || i == num_points - 1)) {
			Plot plot;
			for (double s = 0.; s <= 7; s += 0.05) {
				plot.push(s, v_integrand(s));
			}
			integrands->push_back(plot);
		}

		// Zs flow equation
		if (/*i==norm_point ||*/i == 0) {
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
		if (i == 0) {
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

		if (is_stationary) {
			//std::cout<<der_vals.back()*4/d*vd<<' '<<scaling_vals.back()<<'\n';
		}
	}

	StepFunction rho(step, num_points, [](double x) {
		return x;
	});
	StepFunction v_der = (2 - eta) * V - (d - 2 + eta) * rho * V2 - vd * StepFunction(step, v_vals);
	StepFunction zs_der = -eta * Zs - (d - 2 + eta) * rho * Zs1 - vd * StepFunction(step, zs_vals);
	StepFunction zp_der = -eta * Zp - (d - 2 + eta) * rho * Zp1 - vd * StepFunction(step, zp_vals);

	//double normalization;
	/*if (sigma_normalization) {
	 zs_der.vals[norm_point] = 0.;
	 if (norm_point == 0) {
	 zp_der.vals[0] = 0.;
	 }
	 //normalization = zs_der.vals[norm_point];
	 } else {
	 zp_der.vals[norm_point] = 0.;
	 //normalization = zp_der.vals[norm_point];
	 }*/

	zs_der.vals[0] = zp_der.vals[0];
	/*zs_der = zs_der - normalization;
	 zp_der = zp_der - normalization;*/

	System derivative = *this;
	derivative.parameters = std::vector<StepFunction>();
	derivative.parameters.push_back(v_der);
	derivative.parameters.push_back(zs_der);
	derivative.parameters.push_back(zp_der);
	return derivative;
}

void System::rk4_step(StepFunction *y_der) {

	System second, third, fourth;
	System first_eval, second_eval, third_eval;

	first_eval = this->time_derivative();

	if (y_der != nullptr) {
		*y_der = first_eval.Zp();
	}

	second = *this + delta_t / 2. * first_eval;
	second_eval = second.time_derivative();

	third = *this + delta_t / 2 * second_eval;
	third_eval = third.time_derivative();

	fourth = *this + delta_t * third_eval;

	*this += (delta_t / 6.) * (first_eval + 2. * second_eval + 2. * third_eval + fourth.time_derivative());

	z_dim *= 1 + (delta_t / 6.) * (this->eta + 2 * second.eta + 2 * third.eta + fourth.eta);
	time += delta_t;
	step++;
}

void System::ssp_rk3_step(StepFunction *y_der) {

	System second, third;
	System first_eval, second_eval;

	first_eval = this->time_derivative();

	if (y_der != nullptr) {
		*y_der = first_eval.Zp();
	}

	second = *this + delta_t * first_eval;
	second_eval = second.time_derivative();

	third = *this + delta_t / 4. * (first_eval + second_eval);

	*this += (delta_t / 6.) * (first_eval + second_eval + 4. * third.time_derivative());

	rescale();

	z_dim *= 1 + (delta_t / 6.) * (this->eta + second.eta + 4 * third.eta);
	time += delta_t;
	step++;
}

void System::euler_step(StepFunction *y_der) {

	System first_eval;

	first_eval = this->time_derivative();

	if (y_der != nullptr) {
		*y_der = first_eval.Zp();
	}

	*this += delta_t * first_eval;

	z_dim *= 1 + delta_t * this->eta;
	time += delta_t;
	step++;
}

void System::rescale() {
	double scaling_factor = 1.;
	if (sigma_normalization) {
		scaling_factor = 1 / Zs()[norm_point];
	} else {
		scaling_factor = 1 / Zp()[norm_point];
	}
	for (auto &&f : parameters) {
		f = f * scaling_factor;
	}
}

void System::fix_v() {
	StepFunction &V = this->V();
	for (size_t i = 0; i < V.num_points; i++) {
		/*if (V[i] < -1.95) {
		 V[i] = -1.95;
		 }*/
		if (i > 0 && V[i] < V[i - 1]) {
			V[i] = V[i - 1];
		}
	}
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
	configuration << "-norm," << norm_point << ',';
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
		out << i * s.V().step_size << ',' << s.V()[i] << ',' << s.Zs()[i] << ',' << s.Zp()[i] << '\n';
	}

	return out;
}

double time_eta_distance(System lhs, System rhs) {
	double t_dist = lhs.time - rhs.time;
	double eta_dist = lhs.eta - rhs.eta;

	return std::sqrt(t_dist * t_dist + eta_dist * eta_dist);
}
