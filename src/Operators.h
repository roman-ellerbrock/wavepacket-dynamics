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

/// Simple n-dimensional harmonic oscillator
Hamiltonian harmonic_osciallator(const Basis& basis);

/// Scatter a wave packet at an Eckhard Barrier
Hamiltonian eckhard_pes(const Basis& basis);

/// Free particle (can be used in other Hamiltonians)
Hamiltonian kinetic_energy(const Basis& basis);

/// n-dimensional HO but V is evaluated using a PES
Hamiltonian harmonic_osciallator_pes(const Basis& basis);

/// Model system for diabatic coupling (Tully 1990, Sect. IV A)
Hamiltonian tully_A(const Basis& basis);

#endif //OPERATORS_H
