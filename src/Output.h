//
// Created by Roman Ellerbrock on 8/25/21.
//

#ifndef OUTPUT_H
#define OUTPUT_H
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "Basis.h"
#include "PrimitiveOperators.h"

/// Print information about wavefunction to cout and plot the wavefunction to os_psi if provided
void output(const Wavefunction& Psi, const Hamiltonian& H, const Basis& basis, ostream* os_psi = nullptr);

#endif //OUTPUT_H
