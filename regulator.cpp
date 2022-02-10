#include <cmath>
#include "regulator.h"

// Regulator function and derivatives divided by Z (variable regulating fluctuactions)

const PhysicalDouble long_cutoff = 30.;
const PhysicalDouble short_cutoff = 1e-2;

PhysicalDouble R(PhysicalDouble y) {
	if (y < short_cutoff) {
		PhysicalDouble y2 = y * y, y4 = y2 * y2;
		return 2. - y + y2 / 6. - y4 / 360. + y2 * y4 / 15120.;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	return 2 * y / expm1(y);
}

PhysicalDouble R(PhysicalDouble y, PhysicalDouble expm) {
	if (y < short_cutoff) {
		PhysicalDouble y2 = y * y, y4 = y2 * y2;
		return 2. - y + y2 / 6. - y4 / 360. + y2 * y4 / 15120.;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	return 2 * y / expm;
}

PhysicalDouble Rp(PhysicalDouble y) {
	if (y < short_cutoff) {
		PhysicalDouble y3 = y * y * y, y5 = y3 * y * y;
		return -1. + y / 3 - y3 / 90 + y5 / 2520;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	PhysicalDouble expm = expm1(y);
	PhysicalDouble expy = expm + 1;
	PhysicalDouble denominator = std::pow(expm, -2);

	return 2 * (expm - y * expy) * denominator;
}

PhysicalDouble Rp(PhysicalDouble y, PhysicalDouble expm) {
	if (y < short_cutoff) {
		PhysicalDouble y3 = y * y * y, y5 = y3 * y * y;
		return -1. + y / 3 - y3 / 90 + y5 / 2520;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	PhysicalDouble expy = expm + 1;
	PhysicalDouble denominator = std::pow(expm, -2);

	return 2 * (expm - y * expy) * denominator;
}

PhysicalDouble Rp2(PhysicalDouble y) {
	if (y < short_cutoff) {
		PhysicalDouble y2 = y * y, y4 = y2 * y2;
		return 1. / 3 - y2 / 30 + y4 / 504 + y4 * y2 / 10800;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	PhysicalDouble expm = expm1(y);
	PhysicalDouble expy = expm + 1;
	PhysicalDouble denominator = std::pow(expm, -3);

	return 2 * expy * (-2 * expm + y * (expy + 1)) * denominator;
}

PhysicalDouble Rp2(PhysicalDouble y, PhysicalDouble expm) {
	if (y < short_cutoff) {
		PhysicalDouble y2 = y * y, y4 = y2 * y2;
		return 1. / 3 - y2 / 30 + y4 / 504 + y4 * y2 / 10800;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	PhysicalDouble expy = expm + 1;
	PhysicalDouble denominator = std::pow(expm, -3);

	return 2 * expy * (-2 * expm + y * (expy + 1)) * denominator;
}

PhysicalDouble Prefactor(PhysicalDouble y) {
	if (y < short_cutoff) {
		PhysicalDouble y2 = y * y, y4 = y2 * y2;
		return 4 - y2/3. + y4/60 - y2*y4/1512 + y4*y4/43200.;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	PhysicalDouble expm = expm1(y);
	PhysicalDouble expy = expm + 1;
	PhysicalDouble denominator = std::pow(expm, -2);

	return 4 * y * y * expy * denominator;
}

PhysicalDouble Prefactor(PhysicalDouble y, PhysicalDouble expm) {
	if (y < short_cutoff) {
		PhysicalDouble y2 = y * y, y4 = y2 * y2;
		return 4 - y2/3. + y4/60 - y2*y4/1512 + y4*y4/43200.;
	}
	if (y > long_cutoff) {
		return 0.;
	}
	PhysicalDouble expy = expm + 1;
	PhysicalDouble denominator = std::pow(expm, -2);

	return 4 * y * y * expy * denominator;
}
