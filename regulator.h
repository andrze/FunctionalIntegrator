#ifndef REGULATOR_H
#define REGULATOR_H
#include "realvector.h"

PhysicalDouble R(PhysicalDouble y); //  r(y))
PhysicalDouble R(PhysicalDouble y, PhysicalDouble expm1); //  r(y))
PhysicalDouble Rp(PhysicalDouble y); // r'(y)
PhysicalDouble Rp(PhysicalDouble y, PhysicalDouble expm1); //  r(y))
PhysicalDouble Rp2(PhysicalDouble y); // r''(y)
PhysicalDouble Rp2(PhysicalDouble y, PhysicalDouble expm1); //  r(y))
PhysicalDouble Prefactor(PhysicalDouble y); // 2 (r(y) - y r'(y))
PhysicalDouble Prefactor(PhysicalDouble y, PhysicalDouble expm1); //  r(y))


#endif // REGULATOR_H
