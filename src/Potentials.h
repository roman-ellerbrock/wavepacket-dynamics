//
// Created by Roman Ellerbrock on 8/22/21.
//

#ifndef POTENTIALS_H
#define POTENTIALS_H
#include "Core/Vector.h"

/// Harmonic Oscillator potential as PES
double V_HO(const Vectord& x);

/// Eckhard Barrier
double V_Eckhard(const Vectord& x);

/// Model Potential from Sec. IV A in Tully 1990
double V_tullyA_diabatic11(const Vectord& x);
double V_tullyA_diabatic12(const Vectord& x);
double V_tullyA_diabatic22(const Vectord& x);

double V_tullyA_adiabatic1(const Vectord& x);
#endif //POTENTIALS_H
