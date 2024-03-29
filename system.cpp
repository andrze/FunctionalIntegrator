#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include "system.h"
#include "regulator.h"
#include "gaussquadrature.h"
#include "integrator.h"
#include "realvector.h"
#include "numericalmethods.h"
#include "terminalplot.h"

SystemDelta::SystemDelta() {

}

SystemDelta::SystemDelta(std::vector<StepFunction> parameters) :
	parameters(parameters) {

}

size_t SystemDelta::num_functions() {
	return parameters.size();
}

RealVector SystemDelta::vector_representation() {
	std::vector<PhysicalDouble> params;
	for (size_t i = 0; i < num_functions(); i++) {
		params.insert(params.end(), parameters[i].begin(), parameters[i].end());
	}
	return RealVector(params);
}

StepFunction& SystemDelta::operator[](size_t i) {
	return parameters[i];
}

void SystemDelta::plot_parameters() {
	TerminalPlot plot;
	plot.min_val = -.1;
	plot.max_val = .1;
	plot.fixed_range = true;
	plot.plot(parameters);
}

SystemDelta& SystemDelta::operator+=(SystemDelta rhs) {
	if (num_functions() != rhs.num_functions()) {
		throw std::invalid_argument("System: Adding systems with different numbers of functions");
	}

	for (size_t i = 0; i < num_functions(); i++) {
		parameters[i] = parameters[i] + rhs[i];
	}
	return *this;
}

SystemDelta& SystemDelta::operator*=(PhysicalDouble rhs) {
	for (size_t i = 0; i < num_functions(); i++) {
		parameters[i] = parameters[i] * rhs;
	}
	return *this;
}

SystemDelta operator+(SystemDelta lhs, SystemDelta rhs) {
	lhs += rhs;
	return lhs;
}

SystemDelta operator-(SystemDelta lhs, PhysicalDouble rhs) {
	for (size_t i = 0; i < lhs.num_functions(); i++) {
		lhs[i] = lhs[i] - rhs;
	}
	return lhs;
}

SystemDelta operator*(PhysicalDouble lhs, SystemDelta rhs) {
	rhs *= lhs;
	return rhs;
}

SystemDelta operator/(SystemDelta lhs, SystemDelta rhs) {
	for (size_t i = 0; i < lhs.num_functions(); i++) {
		lhs[i] = lhs[i] / rhs[i];
	}
	return lhs;
}

System::System() {
}

System::System(Integrator *integrator, StepFunction V, PhysicalDouble delta_t, PhysicalDouble d, PhysicalDouble n,
	bool sigma_normalization) :
	integrator(integrator), delta_t(delta_t), d(d), n(n), sigma_normalization(sigma_normalization), num_points(
		V.num_points), step_size(V.step_size) {
	parameters.push_back(V);
	parameters.push_back(StepFunction(V.step_size, V.num_points, [](PhysicalDouble) {
		return 1;
	}));
	parameters.push_back(StepFunction(V.step_size, V.num_points, [](PhysicalDouble) {
		return 1;
	}));

}

System::System(Integrator *integrator, std::vector<std::string> configuration, PhysicalDouble kappa) :
	integrator(integrator) {

	for (size_t i = 0; i < configuration.size() - 1; i++) {
		std::string opt = configuration[i], val = configuration[i + 1];
		if (opt == "-kappa") {
			this->kappa = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-u") {
			u = std::strtod(val.c_str(), nullptr);
		} else if (opt == "-delta") {
			delta_t = std::strtod(val.c_str(), nullptr);
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
		} else if (opt == "-norm_point") {
			norm_point = size_t(std::atoi(val.c_str()));
		}
	}

	if (kappa > 0) {
		this->kappa = kappa;
	}

	step_size = this->kappa * rho0_to_rhomax / num_points;

	reparametrize();
}

System::System(Integrator *integrator, std::string filename) :
	integrator(integrator) {
	std::ifstream file(filename.c_str(), std::ios::in);
	//std::ifstream file("/home/andrzej/Documents/Uczelnia/Anizotropie/IsingFunctional/Release/fixed_point_2.00.csv", std::ios::in);

	if (!file.is_open()) {
		throw std::runtime_error("Fixed point file " + filename + " does not exist");
	}
	std::vector<std::string> config;

	std::string line, line2;
	std::getline(file, line);

	std::istringstream lineStream(line);
	for (std::string cell; std::getline(lineStream, cell, ',');) {
		config.push_back(cell);
	}

	*this = System(integrator, config);

	std::vector<std::vector<PhysicalDouble> > functions;
	for (std::string line; std::getline(file, line, '\n');) {

		std::stringstream lineStream(line);
		std::string cell;
		for (size_t i = 0; std::getline(lineStream, cell, ','); i++) {
			PhysicalDouble val;
			try {
				val = std::strtod(cell.c_str(), nullptr);
			} catch (const std::invalid_argument&) {
				continue;
			}
			if (i >= functions.size()) {
				functions.emplace_back();
			}
			functions[i].push_back(val);
		}
	}

	parameters.clear();
	for (size_t j = 0; j < functions.size(); j++) {
		parameters.emplace_back(kappa * rho0_to_rhomax / num_points, functions[j]);
	}

	file.close();
}

void System::reparametrize() {

	StepFunction V(step_size, num_points, [&](PhysicalDouble x) {
		return u * (x - kappa);
	});

	parameters.clear();
	parameters.push_back(V);
	parameters.push_back(StepFunction(step_size, num_points, [=](PhysicalDouble) {
		return 1;
	}));
	parameters.push_back(StepFunction(step_size, num_points, [=](PhysicalDouble) {
		return 1;
	}));
}

size_t System::num_functions() {
	return parameters.size();
}

StepFunction System::operator[](size_t i) {
	return parameters[i];
}

StepFunction System::V() {
	return parameters[0];
}

StepFunction System::Zs() {
	return parameters[1];
}

StepFunction System::Zp() {
	return parameters[2];
}

RealVector System::full_vector_representation() {
	std::vector<PhysicalDouble> params;
	for (size_t i = 0; i < num_functions(); i++) {
//		if (i != 2) {
//			params.insert(params.end(), parameters[i].begin(), parameters[i].end());
//		} else {
//			params.insert(params.end(), parameters[i].begin() + 1, parameters[i].end());
//		}
		params.insert(params.end(), parameters[i].begin(), parameters[i].end());
	}
	return RealVector(params);
}

void System::precalculate_rho_derivatives() {

	if (!rho_derivatives_calculated) {
		V1 = V().derivative(1);
		V2 = V().derivative(2);
		Zs1 = Zs().derivative(1);
		Zs2 = Zs().derivative(2);
		Zp1 = Zp().derivative(1);
		Zp2 = Zp().derivative(2);

		rho_func = V().x_func();
		d_inv = 1.l / d;

		vd = std::pow(2., -1 - d) * std::pow(M_PI, -d / 2) / std::tgamma(d / 2);

	}
	rho_derivatives_calculated = true;
}

void System::push_time_derivative_integrals(size_t i) {

	PhysicalDouble rho = rho_func[i];
	PhysicalDouble v = V()[i], v1 = V1[i], v2 = V2[i];
	PhysicalDouble zs = Zs()[i], zs1 = Zs1[i], zs2 = Zs2[i];
	PhysicalDouble zp = Zp()[i], zp1 = Zp1[i], zp2 = Zp2[i];

	if (cached_regulator_vals.size() < y_max / y_step) {
		cache_regulator();
	}

	for (size_t k = 0; k <= y_max / y_step; k++) {
		PhysicalDouble y = k * y_step, r = cached_regulator_vals[k];
		if (v + zp * y * y + r < 0.005) {
			throw std::runtime_error(
				"Equations contain a singular propagator at rho=" + std::to_string(rho) + ", y=" + std::to_string(y));
		}
	}

	PhysicalDouble eta_copy = eta;
	if (i == 0) {
		PhysicalDouble f5 = ((-1 + n) * zp1) + zs1;

		auto integrands = [=](IntegrandArgument args) {
			PhysicalDouble y2 = args[0], yd = args[1], ry = args[2], prefactor = args[5];
			PhysicalDouble pref = yd * (prefactor - eta_copy * ry);
			PhysicalDouble g = (v + ry + y2 * zs);
			PhysicalDouble f0 = g * g;

			PhysicalDouble v_int = (-2 * ((2 + n) * v1 + f5 * y2)) / f0;
			PhysicalDouble zs_int = (-2 * f5) / f0;
//			zs_int = 0; //LPA test
			return std::array<PhysicalDouble, 3> { pref * v_int, pref * zs_int, pref * zs_int };
		};
		integrator->push_integrand_function(integrands);
	} else {
		PhysicalDouble f3 = 2 * rho;
		PhysicalDouble f6 = -1 + n;
		PhysicalDouble f19 = 2 * d;
		PhysicalDouble f21 = 4 + d;
		PhysicalDouble f22 = 8 * f6;
		PhysicalDouble f23 = 8 * rho;
		PhysicalDouble f29 = d * v1;
		PhysicalDouble f31 = zp1 * zp1;
		PhysicalDouble f38 = 4 * d;
		PhysicalDouble f117 = rho * zp1;
		PhysicalDouble f119 = rho * rho;

		auto integrands = [=](IntegrandArgument args) {
			PhysicalDouble y2 = args[0], yd = args[1], ry = args[2], rpy = args[3], rpy2 = args[4], prefactor = args[5];
			PhysicalDouble pref = yd * (prefactor - eta_copy * ry);
			PhysicalDouble gp = v + ry + y2 * zp;
			PhysicalDouble gs = v + 2 * rho * v1 + ry + y2 * zs;

			PhysicalDouble f0 = y2 * zp1;
			PhysicalDouble f1 = gp * gp;
			PhysicalDouble f4 = y2 * zs1;
			PhysicalDouble f5 = gs * gs;
			PhysicalDouble f9 = v1 + f0;
			PhysicalDouble f13 = ((3 * v1) + (f3 * v2)) + f4;
			PhysicalDouble f17 = f1 * gp;
			PhysicalDouble f18 = f5 * gs;
			PhysicalDouble f20 = rpy + zp;
			PhysicalDouble f25 = rpy + zs;
			PhysicalDouble f30 = rho * y2;
			PhysicalDouble f34 = f1 * rho;
			PhysicalDouble f47 = f20 * f20;
			PhysicalDouble f56 = f25 * f25;
			PhysicalDouble f61 = d * f18;
			PhysicalDouble f92 = (((3 * d) * v1) + ((f19 * rho) * v2)) + (f21 * f4);
			PhysicalDouble f116 = 8 * y2;
			PhysicalDouble f118 = d * f34;
			PhysicalDouble f134 = (d * gp) * rho;
			PhysicalDouble f135 = f20 + f117;

			PhysicalDouble v_int = (-2 * f13) / f5 - (2 * f6 * f9) / f1;
			PhysicalDouble zs_int = (32 * (f13 * f13) * f30 * f56) / (d * f18 * f5)
				+ (32 * f30 * f47 * f6 * (f9 * f9)) / (d * f1 * f17)
				+ (f22 * f9 * rho * (-(f20 * (f0 * f21 + f29)) - 2 * f9 * rpy2 * y2)) / (d * f17 * gp)
				- (f13 * f23 * (f92 * rpy + 2 * f13 * rpy2 * y2 + f92 * zs)) / (d * f18 * gs)
				+ (f22 * ((f30 * f31) / d + f9 * (-zp + zs))) / f17
				+ (f23 * (6 * f29 + (1 + f19) * f4 + f38 * rho * v2) * zs1) / f61
				- (2 * f6 * (zp - zs + rho * zs1)) / f34 - (2 * (zs1 + f3 * zs2)) / f5;
			PhysicalDouble zp_int = (f116 * f47) / (d * f17 * rho) + (f116 * f56) / (f61 * rho)
				- (f20 * f38 + f117 * f19 * f6 + 8 * rpy2 * y2) / f118
				+ ((f116 * (f119 * f31 - f47)) / f118 + (8 * d * f135 + 16 * rpy2 * y2) / f134) / gs
				+ ((f116 * f135 * (f117 - rpy + zp - 2 * zs)) / f134
					- (2 * (4 * rpy2 * y2 + d * (5 * f117 + 2 * rpy + zp + 2 * f119 * zp2 + zs))) / (d * rho)) / f5;
//			zs_int = 0; //LPA test
//			zp_int = 0; //LPA test
			return std::array<PhysicalDouble, 3> { pref * v_int, pref * zs_int, pref * zp_int };
		};
		integrator->push_integrand_function(integrands);
	}

}

void System::find_eta() {
	PhysicalDouble rho = rho_func[norm_point];

	integrator->reset_integrals(2);

	eta = 0;
	push_time_derivative_integrals(norm_point);
	eta = 1;
	push_time_derivative_integrals(norm_point);

	auto results = integrator->evaluate_integrals();

	IntegrandValue i_free = (*results)[0];
	IntegrandValue i_total = (*results)[1];

	PhysicalDouble z, z1;
	PhysicalDouble i1, i2;
	if (sigma_normalization) {
		i1 = i_free[1];
		i2 = i_total[1] - i_free[1];
		z = Zs()[norm_point];
		z1 = Zs1[norm_point];
	} else {
		i1 = i_free[2];
		i2 = i_total[2] - i_free[2];
		z = Zp()[norm_point];
		z1 = Zp1[norm_point];
	}

	PhysicalDouble free_component = (d - 2) * rho * z1 + vd * i1;
	PhysicalDouble eta_component = z + rho * z1 + vd * i2;

	eta = -free_component / eta_component;
//	eta = 0; //LPA test
}

PhysicalDouble System::time_derivative(size_t func_i, size_t rho_i, size_t integrals_position) {

	PhysicalDouble integral = -vd * integrator->integral_result(integrals_position)[func_i];

	PhysicalDouble scaling;
	PhysicalDouble rho = rho_func[rho_i];
	if (func_i == 0) {
		scaling = (2 - eta) * V()[rho_i] - (d - 2 + eta) * rho * V1[rho_i];
	} else if (func_i == 1) {
		scaling = -eta * Zs()[rho_i] - (d - 2 + eta) * rho * Zs1[rho_i];
	} else if (func_i == 2) {
		scaling = -eta * Zp()[rho_i] - (d - 2 + eta) * rho * Zp1[rho_i];
	} else {
		throw std::runtime_error(
			"System::time_derivative: Scaling for function " + std::to_string(func_i) + " not implemented");
	}

	return integral + scaling;
}

SystemDelta System::time_derivative() {
	precalculate_rho_derivatives();
	find_eta();

	size_t num_points = V().num_points;
	PhysicalDouble step = V().step_size;

	integrator->reset_integrals(num_points);

	for (size_t i = 0; i < num_points; i++) {
		push_time_derivative_integrals(i);
	}
	auto integrals = integrator->evaluate_integrals();

	std::vector<PhysicalDouble> v_vals(num_points), zs_vals(num_points), zp_vals(num_points);
	for (size_t i = 0; i < num_points; i++) {
		auto &&res = (*integrals)[i];
		v_vals[i] = res[0];
		zs_vals[i] = res[1];
		zp_vals[i] = res[2];
	}

	StepFunction v_der = (2 - eta) * V() - (d - 2 + eta) * rho_func * V1
		- vd * StepFunction(step, v_vals, V().domain_begin);
	StepFunction zs_der = -eta * Zs() - (d - 2 + eta) * rho_func * Zs1
		- vd * StepFunction(step, zs_vals, V().domain_begin);
	StepFunction zp_der = -eta * Zp() - (d - 2 + eta) * rho_func * Zp1
		- vd * StepFunction(step, zp_vals, V().domain_begin);

	beta_function = SystemDelta( { v_der, zs_der, zp_der });
	beta_squared = beta_function.vector_representation().norm();

	return beta_function;
}

void System::rescale() {
	PhysicalDouble potential_minimum = V().kappa_u().first;
	if (potential_minimum <= std::numeric_limits<PhysicalDouble>::epsilon()) {
		throw std::runtime_error("System: Potential minimum is out of the grid");
	}
	PhysicalDouble z_norm;
	if (sigma_normalization) {
		z_norm = Zs()[norm_point];
	} else {
		z_norm = Zp()[norm_point];
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

	z_dim *= z_norm;
	z_correction *= z_norm;
}

void System::zoom_in() {
	std::vector<PhysicalDouble> v_vals(V().begin(), V().end()), xs = V().xs();
	const PhysicalDouble V_LOWER_THRESHOLD = -a;
	const PhysicalDouble V_UPPER_THRESHOLD = 2e+6;
	const PhysicalDouble V_LOWER_GOAL = -a;
	const PhysicalDouble V_UPPER_GOAL = 2e+6;

	if (v_vals[0] > V_LOWER_THRESHOLD && v_vals[0] < v_vals[1] && v_vals[v_vals.size() - 1] < V_UPPER_THRESHOLD) {
		return;
	}

	size_t rescale_begin_index = 0, rescale_end_index = xs.size() - 1;

	for (size_t i = 1; i < v_vals.size(); i++) {
		if (v_vals[i] > V_LOWER_GOAL && v_vals[i] > v_vals[i - 1]) {
			rescale_begin_index = i;
			break;
		}
	}

	for (size_t i = rescale_begin_index + 1; i < v_vals.size(); i++) {
		if (v_vals[i] > V_UPPER_GOAL) {
			rescale_end_index = i - 1;
			break;
		}
	}

	PhysicalDouble rescale_begin = xs[rescale_begin_index];
	PhysicalDouble rescale_end = xs[rescale_end_index];
	for (auto &&func : this->parameters) {
		func = func.zoom_in(rescale_begin, rescale_end);
	}
	std::cout << "Zooming in on the interval [" << rescale_begin << ',' << rescale_end << "].\n";
	zoomed = true;
}

void System::cut_domain() {
	const double LOWER_THRESHOLD = -a * 2;
	const double UPPER_THRESHOLD = 5e+4;

	if (*V().begin() > LOWER_THRESHOLD && *std::next(V().end(), -1) < UPPER_THRESHOLD) {
		return;
	}

	std::vector<PhysicalDouble> v_vals = V().vals_copy(), xs = V().xs();

	size_t rescale_begin = 0, rescale_end = v_vals.size();

	if (v_vals[0] < LOWER_THRESHOLD) {
		for (size_t i = 1; i < v_vals.size(); i++) {
			if (v_vals[i] > LOWER_THRESHOLD && v_vals[i] > v_vals[i - 1]) {
				rescale_begin = i;
				break;
			}
		}
	}

	if (v_vals.back() > UPPER_THRESHOLD) {
		for (size_t i = rescale_begin + 1; i < v_vals.size(); i++) {
			if (v_vals[i] > UPPER_THRESHOLD) {
				rescale_end = i;
				break;
			}
		}
	}

	for (auto &&func : this->parameters) {
		func = func.cut_domain(rescale_begin, rescale_end);
	}
	std::cout << "Zooming in on the interval [" << xs[rescale_begin] << ',' << xs[rescale_end] << "].\n";
	zoomed = true;
	rho_derivatives_calculated = false;
}

PhysicalDouble System::last_val() {
	return *std::next(V().end(), -1);
}

PhysicalDouble System::first_val() {
	return *V().begin();
}

int System::find_phase() {
	phase = 0;
	PhysicalDouble v_at_point_1 = V()[1];
	PhysicalDouble v_at_point_m1 = V()[V().num_points - 2];
	if (v_at_point_1 > 0) { // Disordered phase
		phase = -1;
	} else if (v_at_point_m1 < 0) { // Ordered phase
		phase = 1;
	}
	return phase;
}

std::string System::print_phase() {
	if (phase == 0) {
		if (this->phase != 0) {
			phase = this->phase;
		} else {
			phase = find_phase();
		}
	}

	std::string message = "Flow ended in the ";
	if (phase == 1) {
		message += "ordered";
	} else if (phase == -1) {
		message += "disordered";
	} else {
		message += "unidentified";
	}
	message += " phase.\n";

	return message;
}

std::array<PhysicalDouble, 3> System::kappa_u_z() {
	precalculate_rho_derivatives();
	auto kappa_u = V().kappa_u(V1);

	if (std::abs(kappa_u.first - (-1)) < std::numeric_limits<PhysicalDouble>::epsilon()) {
		return {-1, -1, -1};
	}

	PhysicalDouble z;
	if (sigma_normalization) {
		z = Zs()(kappa_u.first);
	} else {
		z = Zp()(kappa_u.first);
	}

	return std::array<PhysicalDouble, 3> { { kappa_u.first, kappa_u.second, z } };
}

std::string System::print_configuration() {
	std::stringstream configuration;
	configuration << "-delta," << delta_t << ',';
	configuration << "-dim," << d << ',';
	configuration << "-N," << n << ',';
	configuration << "-num_points," << num_points << ',';
	configuration << "-a," << a << ',';
	configuration << "-sigma_normalization," << sigma_normalization << ',';
	configuration << "-norm_point," << norm_point << ',';
	configuration << "-rho_max," << rho0_to_rhomax << '\n';

	return configuration.str();
}

System& System::operator+=(SystemDelta rhs) {
	if (num_functions() != rhs.num_functions()) {
		throw std::invalid_argument("System: Adding systems with different numbers of functions");
	}

	for (size_t i = 0; i < num_functions(); i++) {
		parameters[i] = parameters[i] + rhs[i];
	}
	rho_derivatives_calculated = false;
	return *this;
}

void System::cache_regulator() {
	cached_regulator_vals.clear();
	for (size_t i = 0; i <= y_max / y_step; i++) {
		PhysicalDouble y = i * y_step;
		cached_regulator_vals.push_back(a / 2 * R(y * y));
	}
}

void System::plot_parameters() {
	TerminalPlot plot;
	plot.plot(parameters);
}

std::ostream& operator<<(std::ostream &out, System s) {
	out << "Rho,V,Zs,Zp,V1,Zs1,Zp1,V2,Zs2,Zp2\n";
	s.precalculate_rho_derivatives();
	for (size_t i = 0; i < s.V().num_points; i++) {
		out << s.V().x_func()[i] << ',' << s.V()[i] << ',' << s.Zs()[i] << ',' << s.Zp()[i] << ',' << s.V1[i] << ','
			<< s.Zs1[i] << ',' << s.Zp1[i] << ',' << s.V2[i] << ',' << s.Zs2[i] << ',' << s.Zp2[i] << '\n';
	}

	return out;
}

PhysicalDouble time_eta_distance(System lhs, System rhs) {
	PhysicalDouble t_dist = lhs.time - rhs.time;
	PhysicalDouble eta_dist = lhs.eta - rhs.eta;

	return std::sqrt(t_dist * t_dist + eta_dist * eta_dist);
}
