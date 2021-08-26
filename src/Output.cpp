//
// Created by Roman Ellerbrock on 8/25/21.
//

#include "Output.h"

/**
 * \brief Calculate the expectation value <Psi|P|Psi>
 * @param Psi
 * @param P
 * @param basis
 * @return expectation value <P>
 */
double expect(const Wavefunction& Psi, const ProductOperator& P, const Basis& basis) {
	Wavefunction PPsi = P.apply(Psi, basis);
	Matrixcd Prep = Psi.dotProduct(PPsi);
	return real(Prep(0,0));
}

/**
 * \brief Calculate energy expectation value <Psi|H|Psi>
 * @param Psi
 * @param H
 * @param basis
 * @return expected energy E = <H>
 */
double energy_expect(const Wavefunction& Psi, const Hamiltonian& H, const Basis& basis) {
	Wavefunction HPsi = H.apply(Psi, basis);
	Matrixcd Hrep = Psi.dotProduct(HPsi);
	return real(Hrep(0,0));
}

/**
 * Calculate expectation values of wavepacket and print
 * @param Psi
 * @param H
 * @param basis
 */
void output(const Wavefunction& Psi, const Hamiltonian& H, const Basis& basis) {
	/// Note: This is a very simple output routine.
	/// For number Basis, x is the particle number operator

	/// Calculate total Energy
	double energy = energy_expect(Psi, H, basis);
	cout << "<H> = " << energy << endl;
	for (const PrimitiveBasis& prim : basis) {
		/// Calculate kinetic Energy
		ProductOperator T(prim.kin_, prim.coord_);
		double kin = expect(Psi, T, basis);
		ProductOperator X(prim.x_, prim.coord_);
		double x = expect(Psi, X, basis);
		X.emplace_back(prim.x_, prim.coord_);
		double x2 = expect(Psi, X, basis);
		double dx = sqrt(abs(x*x - x2));

		/// Make Output:
		cout << "coord: " << prim.coord_ << ", <T> = " << kin << ",\t <x> = " << x << ",\t <dx> = " << dx << endl;
	}
}

