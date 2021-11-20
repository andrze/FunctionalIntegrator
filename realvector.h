#ifndef REALVECTOR_H
#define REALVECTOR_H

#include <utility>
#include <vector>
#include <ostream>
#include <functional>

typedef long double PhysicalDouble;

class RealVector {
public:
	RealVector();
	RealVector(std::vector<PhysicalDouble> coords);

	std::vector<PhysicalDouble> coords;

	size_t size();
	RealVector& operator +=(RealVector rhs);
	RealVector& operator *=(PhysicalDouble rhs);
	RealVector& operator -=(RealVector rhs);
	RealVector& operator /=(PhysicalDouble rhs);
	PhysicalDouble& operator [](size_t i);
};

RealVector operator +(RealVector lhs, RealVector rhs);
RealVector operator -(RealVector lhs, RealVector rhs);
RealVector operator *(RealVector lhs, PhysicalDouble rhs);
RealVector operator *(PhysicalDouble lhs, RealVector rhs);
RealVector operator /(RealVector lhs, PhysicalDouble rhs);
PhysicalDouble operator *(RealVector lhs, RealVector rhs);
std::ostream& operator <<(std::ostream &out, RealVector v);

RealVector filter(RealVector v, std::vector<bool>);

#endif // REALVECTOR_H
