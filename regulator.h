#ifndef REGULATOR_H
#define REGULATOR_H
#include "realvector.h"

PhysicalDouble R(PhysicalDouble y); //  r(y))
PhysicalDouble Rp(PhysicalDouble y); // r'(y)
PhysicalDouble Rp2(PhysicalDouble y); // r''(y)
PhysicalDouble Prefactor(PhysicalDouble y); // 2 (r(y) - y r'(y))


#endif // REGULATOR_H
