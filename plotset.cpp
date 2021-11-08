/*
 * PlotSet.cpp
 *
 *  Created on: 15 gru 2018
 *      Author: andrzejchlebicki
 */

#include "plotset.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <iostream>

PlotSet::PlotSet(size_t i) {
	for(size_t j=0; j<i; j++){
		plots.push_back(Plot());
	}
}

PlotSet::~PlotSet() {
}

void PlotSet::push_to_each(RealVector values, RealVector ders, double t){
	for(size_t i=0; i<plots.size(); i++){
		plots[i].values.push_back(values[i]);
		plots[i].derivatives.push_back(ders[i]);
		plots[i].times.push_back(t);
	}
}

void PlotSet::pop_from_each(){
	for(size_t i=0; i<plots.size(); i++){
		plots[i].values.pop_back();
		plots[i].derivatives.pop_back();
		plots[i].times.pop_back();
	}
}

RealVector PlotSet::back_vals(){
	std::vector<double> coords;
	for(size_t i=0; i<plots.size(); i++){
		coords.push_back(plots[i].values.back());
	}
	return RealVector(coords);
}

RealVector PlotSet::back_ders(){
	std::vector<double> coords;
	for(size_t i=0; i<plots.size(); i++){
		coords.push_back(plots[i].derivatives.back());
	}
	return RealVector(coords);
}

double PlotSet::back_time(){
	double t=plots[0].times.back();
	double precision = t*std::numeric_limits<double>::epsilon();
	for(size_t i=1; i<plots.size(); i++){
		if(std::abs(t - plots[i].times.back()) > precision){
			throw std::invalid_argument("Simulation time different in plots on same level");
		}
	}
	return t;
}

size_t PlotSet::plot_size(){
	size_t size = plots[0].size();
	for(size_t i=1; i<plots.size(); i++){
		if(plots[i].size() != size){
			throw std::invalid_argument("Sizes of plots in PlotSet differ");
		}
	}
	return size;
}

Plot& PlotSet::operator[](size_t i){

	return plots[i];
}

size_t PlotSet::plot_number(){
	return plots.size();
}

double PlotSet::eta(size_t k){
    return plots[4].exp_time_log_der(k);
}

Plot PlotSet::rescaled(size_t k, double d){
    if(k > this->plot_number()){
        throw std::invalid_argument("Argument larger than number of plots");
    }

    if(this->plot_number() < 5){
        throw std::invalid_argument("Not enough plots to allow rescaling");
    }

    Plot rescaled;
    rescaled.times = this->plots[k].times;
    for(size_t i=0; i<this->plot_size(); i++){
        double Z = this->plots[4].values[i];
        double T = this->plots[5].values[i];
        double scaling = 1;
        if(k==0){
            scaling = std::sqrt(Z)*std::pow(rescaled.times[i],(2-d)/2)/sqrt(T);
        } else if(k==1 || k==2){
            scaling = T/(Z*rescaled.times[i]*rescaled.times[i]);
        }

        rescaled.values.push_back(this->plots[k].values[i]*scaling);
        rescaled.derivatives.push_back(this->plots[k].derivatives[i]*scaling);
    }

    return rescaled;
}

int PlotSet::phase_diagnosis(double d){
	size_t plot_size = this->plot_size();
	size_t start = plot_size-10;
	if(plot_size < 11){
		start = 1;
	}
    for(size_t i=start; i<plot_size; i++){
        if(std::abs(plots[0].exp_time_log_der(i)) < 1e-04){
        //if(plots[0].values[i] > 1e-0{
            return 2;
        }
    }

    auto alpha = this->rescaled(0,d);
    auto mpi = this->rescaled(2,d);

    if(d!=0 && mpi.values[1]/(alpha.values[1]*alpha.values[1]) < mpi.values[0]/(alpha.values[0]*alpha.values[0])){
        return 1;
    }
    /*
    for(size_t i=0; i<plot_size; i++){
        double k = std::pow(plots[0].values[i],2)*plots[4].values[i];
        if(k > 0.02 && eta(i) > 0.265){
            return 1;
        }
    }*/

    return 0;
}

std::ostream& operator<<(std::ostream& out, PlotSet plots){
	for(size_t i=0; i<plots.plot_size(); i++){
		for(size_t j=0; j<plots.plot_number(); j++){
			auto slice = plots[j][i];
			if(j == 0){
				out << -log(slice[0]) << ", ";
			}
			out << slice[1] << ", " << slice[2] << ", ";
		}
		out << "\n";
	}
	return out;
}
