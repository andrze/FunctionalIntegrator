#include <cmath>
#include "regulator.h"

// Regulator function and derivatives divided by Z (variable regulating fluctuactions)

const double long_cutoff = 30.;
const double short_cutoff = 1e-2;

double R(double y) {
	if (y < short_cutoff) {
		double y2 = y * y, y4 = y2 * y2;
		return 2. - y + y2 / 6. - y4 / 360. + y2 * y4 / 15120.;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	return 2 * y / expm1(y);
}

double Rp(double y) {
	if (y < short_cutoff) {
		double y3 = y * y * y, y5 = y3 * y * y;
		return -1. + y / 3 - y3 / 90 + y5 / 2520;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	double expm = expm1(y);
	double expy = expm + 1;
	double denominator = std::pow(expm, -2.);

	return 2 * (expm - y * expy) * denominator;
}

double Rp2(double y) {
	if (y < short_cutoff) {
		double y2 = y * y, y4 = y2 * y2;
		return 1. / 3 - y2 / 30 + y4 / 504 + y4 * y2 / 10800;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	double expm = expm1(y);
	double expy = expm + 1;
	double denominator = std::pow(expm, -3.);

	return 2 * expy * (-2 * expm + y * (expy + 1)) * denominator;

}

double Prefactor(double y) {
	if (y < short_cutoff) {
		double y2 = y * y, y4 = y2 * y2;
		return 4 - y2/3. + y4/60 - y2*y4/1512 + y4*y4/43200.;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	double expm = expm1(y);
	double expy = expm + 1;
	double denominator = std::pow(expm, -2);

	return 4 * y * y * expy * denominator;
}

double Prefactor(double y, double eta) {
	if (y < short_cutoff) {
		double y2 = y * y, y3 = y2 * y, y4 = y2 * y2, y5 = y4 * y, y6 = y4 * y2;
		return -2 * (-2 + eta) + eta * y + ((-2 - eta) * y2) / 6.
				+ 2 * (-0.6666666666666666 + (5 * (2 - eta / 2.)) / 12. + (-2 + eta) / 12. + eta / 8.) * y3
				+ 2 * (-0.25 + (2 - eta) / 240. + (5 * (1 - eta / 6.)) / 12. + (-2 + eta / 2.) / 12. + eta / 30.) * y4
				+ 2
						* (-0.06666666666666667 + (2 - eta) / 720. + (2 - eta / 2.) / 240.
								+ (5 * (0.3333333333333333 - eta / 24.)) / 12. + (-1 + eta / 6.) / 12. + eta / 144.)
						* y5
				+ 2
						* (-0.013888888888888888 + (2 - eta / 2.) / 720. + (1 - eta / 6.) / 240.
								+ (5 * (0.08333333333333333 - eta / 120.)) / 12.
								+ (-0.3333333333333333 + eta / 24.) / 12. + (-2 + eta) / 6048. + eta / 840.) * y6;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	double expm = expm1(y);
	double expy = expm + 1;
	double denominator = std::pow(expm, -2);

	return (y * 2 * (expy * (2 * y - eta) + eta)) * denominator;
}
