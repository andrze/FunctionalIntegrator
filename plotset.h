#ifndef PLOTSET_H_
#define PLOTSET_H_

#include "realvector.h"
#include "plot.h"

class PlotSet {
public:
	PlotSet(size_t i=0);
	virtual ~PlotSet();

	void push_to_each(RealVector values, RealVector ders, PhysicalDouble t);
	void pop_from_each();
	RealVector back_vals();
	RealVector back_ders();
	PhysicalDouble back_time();

	size_t plot_size();
	Plot& operator[](size_t i);
	size_t plot_number();

	PhysicalDouble eta(size_t k);
    Plot rescaled(size_t k, PhysicalDouble d);
    int phase_diagnosis(PhysicalDouble d=0);

private:
    std::vector<Plot> plots;

};

std::ostream& operator<<(std::ostream& out, PlotSet plots);

#endif /* PLOTSET_H_ */
