/*
 * TerminalPlot.cpp
 *
 *  Created on: Feb 10, 2022
 *      Author: andrzej
 */

#include "terminalplot.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include "realvector.h"

TerminalPlot::TerminalPlot() {

}

TerminalPlot::~TerminalPlot() {
}

void TerminalPlot::plot(std::vector<StepFunction> functions) {
	if (functions.size() > markers.size()) {
		throw std::invalid_argument("TerminalPlot: Asked to plot more functions than supported");
	}

	PhysicalDouble min = std::numeric_limits<PhysicalDouble>::max(), max = -std::numeric_limits<PhysicalDouble>::min();

	for (auto &&f : functions) {
		auto minmax = f.minmax();

		min = std::min(minmax.first, min);
		max = std::max(minmax.second, max);
	}

	min = std::max(min, min_val);
	max = std::min(max, max_val);

	std::vector<PhysicalDouble> x_vals = functions.front().xs();

	double y_unit = (max - min) / height;
	double x_start = x_vals.front();
	double x_unit = (x_vals.back() - x_start) / width;

	int zero_pos = -min / y_unit;

	std::vector<char> empty_row(size_t(width), ' ');
	std::vector<std::vector<char> > plot_symbols(size_t(height), empty_row);
	for (int j = 0; j < width; j++) {
		double x = x_start + x_unit * j;
		for (int i = 0; i < int(functions.size()); i++) {
			double val = functions[size_t(i)](x);
			int pos = static_cast<size_t>((val - min) / y_unit);
			if (pos < 0 || pos >= height) {
				continue;
			}

			if (plot_symbols[size_t(pos)][size_t(j)] == ' ') {
				plot_symbols[size_t(pos)][size_t(j)] = markers[size_t(i)];
			} else {
				plot_symbols[size_t(pos)][size_t(j)] = '*';
			}
		}
		if (zero_pos >= 0 && zero_pos < height && plot_symbols[size_t(zero_pos)][size_t(j)] == ' ') {
			plot_symbols[size_t(zero_pos)][size_t(j)] = '-';
		}
	}

	int prec = 2;
	std::cout << std::setprecision(prec) << std::fixed;
	for (int j = 0; j < width + 2; j++) {
		std::cout << '#';
	}
	std::cout << '\n';
	for (int i = height - 1; i >= 0; i--) {
		std::cout << "#";
		for (int j = 0; j < width; j++) {
			std::cout << plot_symbols[size_t(i)][size_t(j)];
		}
		std::cout << "#";
		if (i == 0) {
			std::cout << ' ' << min;
		}
		if (i == zero_pos) {
			std::cout << ' ' << 0;
		}
		if (i == height - 1) {
			std::cout << ' ' << max;
		}
		std::cout << "\n";
	}
	for (int j = 0; j < width + 2; j++) {
		std::cout << '#';
	}
	std::cout << '\n';
	for (int j = 0; j < width + 2; j++) {
		if (j == 1) {
			std::cout << x_vals.front();
			j += prec + 2;
		} else if (j == width - prec) {
			std::cout << x_vals.back();
			j += prec + 2;
		} else {
			std::cout << ' ';
		}
	}
	std::cout << '\n';
}

