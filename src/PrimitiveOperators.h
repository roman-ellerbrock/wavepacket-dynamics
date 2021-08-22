//
// Created by Roman Ellerbrock on 8/19/21.
//

#ifndef PRIMITIVEOPERATORS_H
#define PRIMITIVEOPERATORS_H
#include "Core/Matrix.h"

/// Harmonic Oscillator
Matrixcd x_HO(size_t dim, double freq);
Matrixcd p_HO(size_t dim, double freq);
Matrixcd kin_HO(size_t dim, double freq);
pair<Matrixcd, Vectord> dvr_HO(size_t dim, double freq, double x0);

/// Equidistant Grid / FFT
Matrixcd x_FFT(size_t dim, double x0, double x1);
Matrixcd p_FFT(size_t dim, double x0, double x1);
Matrixcd kin_FFT(size_t dim, double x0, double x1);
pair<Matrixcd, Vectord> dvr_FFT(size_t dim, double x0, double x1);

/// Electronic-/Number-basis


#endif //PRIMITIVEOPERATORS_H
