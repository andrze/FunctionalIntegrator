/*
 * TerminalPlot.h
 *
 *  Created on: Feb 10, 2022
 *      Author: andrzej
 */

#ifndef TERMINALPLOT_H_
#define TERMINALPLOT_H_

#include "stepfunction.h"
#include <cstdlib>
#include <vector>

class TerminalPlot {
public:
	TerminalPlot();
	virtual ~TerminalPlot();

	void plot(std::vector<StepFunction> functions);

	PhysicalDouble max_val = 2, min_val = -2;
	bool fixed_range = false;
private:
	int width = 80, height = 25;
	std::vector<char> markers { '1', '2', '3', '4' };
};

#endif /* TERMINALPLOT_H_ */
