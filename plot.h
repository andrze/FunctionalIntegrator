/*
 * plot.h
 *
 *  Created on: 15 gru 2018
 *      Author: andrzej
 */

#ifndef PLOT_H_
#define PLOT_H_
#include <vector>
#include <ostream>
#include <array>
#include <cmath>

class Plot {
public:
	Plot();
	~Plot();

	std::vector<double> values;
	std::vector<double> derivatives;
	std::vector<double> times;

	size_t size();
	double exp_time_log_der(size_t k);
	std::array<double, 3> operator[](size_t k);
	void push(double t, double v, double d = std::nan(""));
};

std::ostream& operator<<(std::ostream &out, Plot p);

#endif /* PLOT_H_ */
