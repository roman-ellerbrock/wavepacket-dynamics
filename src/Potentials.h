//
// Created by Roman Ellerbrock on 8/22/21.
//

#ifndef POTENTIALS_H
#define POTENTIALS_H
#include "Core/Vector.h"

double V_HO(const Vectord& x);

double V_tullyA_diabatic11(const Vectord& x);
double V_tullyA_diabatic12(const Vectord& x);
double V_tullyA_diabatic22(const Vectord& x);

#endif //POTENTIALS_H
