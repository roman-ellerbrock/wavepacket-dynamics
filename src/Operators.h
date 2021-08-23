//
// Created by Roman Ellerbrock on 8/19/21.
//

#ifndef OPERATORS_H
#define OPERATORS_H
#include "Hamiltonian.h"
#include "Basis.h"

Hamiltonian harmonic_osciallator(const Basis& basis);

Hamiltonian kinetic_energy(const Basis& basis);

Hamiltonian harmonic_osciallator_pes(const Basis& basis);

Hamiltonian tully_A(const Basis& basis);

#endif //OPERATORS_H
