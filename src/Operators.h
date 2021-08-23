//
// Created by Roman Ellerbrock on 8/19/21.
//

#ifndef OPERATORS_H
#define OPERATORS_H
#include "Hamiltonian.h"
#include "Basis.h"

/**
 * Rationale:
 * - This file contains Hamiltoinians for different model systems.
 * - Feel free to add more and play around.
 */

Hamiltonian harmonic_osciallator(const Basis& basis);

Hamiltonian kinetic_energy(const Basis& basis);

Hamiltonian harmonic_osciallator_pes(const Basis& basis);

Hamiltonian tully_A(const Basis& basis);

#endif //OPERATORS_H
