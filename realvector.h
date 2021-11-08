#ifndef REALVECTOR_H
#define REALVECTOR_H

#include <utility>
#include <vector>
#include <ostream>
#include <functional>

class RealVector {
public:
	RealVector();
	RealVector(std::vector<double> coords);

	std::vector<double> coords;

	size_t size();
	RealVector& operator +=(RealVector rhs);
	RealVector& operator *=(double rhs);
	RealVector& operator -=(RealVector rhs);
	RealVector& operator /=(double rhs);
	double& operator [](size_t i);
};

RealVector operator +(RealVector lhs, RealVector rhs);
RealVector operator -(RealVector lhs, RealVector rhs);
RealVector operator *(RealVector lhs, double rhs);
RealVector operator *(double lhs, RealVector rhs);
RealVector operator /(RealVector lhs, double rhs);
double operator *(RealVector lhs, RealVector rhs);
std::ostream& operator <<(std::ostream &out, RealVector v);

RealVector filter(RealVector v, std::vector<bool>);

#endif // REALVECTOR_H
