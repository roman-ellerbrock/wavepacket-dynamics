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
double tullyA_diabatic11(const Vectord& q);
double tullyA_diabatic12(const Vectord& q);
double tullyA_diabatic22(const Vectord& q);


double V_tullyA_adiabatic1(const Vectord& q);
double V_tullyA_adiabatic2(const Vectord& q);
#endif //POTENTIALS_H
