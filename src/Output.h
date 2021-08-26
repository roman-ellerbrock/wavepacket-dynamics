//
// Created by Roman Ellerbrock on 8/25/21.
//

#ifndef OUTPUT_H
#define OUTPUT_H
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "Basis.h"
#include "PrimitiveOperators.h"

void output(const Wavefunction& Psi, const Hamiltonian& H, const Basis& basis);

#endif //OUTPUT_H
