#ifndef PLOTSET_H_
#define PLOTSET_H_

#include "realvector.h"
#include "plot.h"

class PlotSet {
public:
	PlotSet(size_t i=0);
	virtual ~PlotSet();

	void push_to_each(RealVector values, RealVector ders, double t);
	void pop_from_each();
	RealVector back_vals();
	RealVector back_ders();
	double back_time();

	size_t plot_size();
	Plot& operator[](size_t i);
	size_t plot_number();

	double eta(size_t k);
    Plot rescaled(size_t k, double d);
    int phase_diagnosis(double d=0);

private:
    std::vector<Plot> plots;

};

std::ostream& operator<<(std::ostream& out, PlotSet plots);

#endif /* PLOTSET_H_ */
