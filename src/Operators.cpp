//
// Created by Roman Ellerbrock on 8/19/21.
//

#include "Operators.h"
#include "Potentials.h"

Hamiltonian harmonic_osciallator(const Basis& basis) {
	Hamiltonian H;
	for (const PrimitiveBasis& prim : basis) {
		size_t coord = prim.coord_;

		/// Kinetic Energy
		{
			ProductOperator P;
			P.emplace_back(prim.kin_, coord);
			H.emplace_back(1., P);
		}

		/// HO Potential
		{
			ProductOperator P;
			P.emplace_back(prim.x_ * prim.x_, coord);
			H.emplace_back(0.5, P);
		}

		/// quartic Potential
		{
			ProductOperator P;
			auto x2 = prim.x_ * prim.x_;
			auto x4 = x2 * x2;
			P.emplace_back(x4 , coord);
			double coeff = 0.01;
			H.emplace_back(coeff, P);
		}
	}
	return H;
}

Hamiltonian kinetic_energy(const Basis& basis) {
	Hamiltonian H;
	for (const PrimitiveBasis& prim : basis) {
		size_t coord = prim.coord_;

		/// Kinetic Energy
		{
			ProductOperator P;
			P.emplace_back(prim.kin_, coord);
			H.emplace_back(1., P);
		}
	}

	return H;
}

Hamiltonian abs_pes(const Basis& basis, double omega) {
	/// T
	Hamiltonian H = kinetic_energy(basis);

	/// V_HO-PES
	{
		ProductOperator P;
		auto f = [omega](const Vectord& x) {
			double V = 0.;
			for (size_t k = 0; k < x.dim(); ++k) {
				V += omega * abs(x(k));
			}
			return V;
		};

		P.V_ = f;
		H.emplace_back(1., P);
	}
	return H;
}

Hamiltonian harmonic_osciallator_pes(const Basis& basis) {
	/// T
	Hamiltonian H = kinetic_energy(basis);

	/// V_HO-PES
	{
		ProductOperator P;
		P.V_ = V_HO;
		H.emplace_back(1., P);
	}
	return H;
}

Hamiltonian eckhard_pes(const Basis& basis) {
	/// T
	Hamiltonian H = kinetic_energy(basis);

	/// V_HO-PES
	{
		ProductOperator P;
		P.V_ = V_Eckhard;
		H.emplace_back(1., P);
	}
	return H;
}

Hamiltonian tully_A_adiabatic(const Basis& basis) {
	/**
	 * Rationale:
	 * - from Sec. IV.) A of Ref: http://dx.doi.org/10.1063/1.459170
	 * - use input: Coordinate 1: FFT grid, Coordinate 2: number basis
	 * - written in mass-weighted coordinates V is V(q), not V(x), where q=sqrt(m)*x
	 * - adiabatic model, i.e. no state interaction
	 */

	/// H = T + |0><0|*V_adiabatic1 + |1><1|*V_adiabatic2
	Hamiltonian H;
	{ /// kinetic energy for FFT
		ProductOperator P;
		P.emplace_back(basis[0].kin_, 0);
		H.emplace_back(1., P);
	}

	size_t nstates = 2;
	{
		ProductOperator P;
		P.emplace_back(element(nstates, 0, 0), 1);
		P.V_ = V_tullyA_adiabatic1;
		H.emplace_back(1., P);
	}

	{
		ProductOperator P;
		P.emplace_back(element(nstates, 1, 1), 1);
		P.V_ = V_tullyA_adiabatic2;
		H.emplace_back(1., P);
	}
	return H;
}

Hamiltonian tully_A(const Basis& basis) {
	/**
	 * Rationale:
	 * - from Sec. IV.) A of Ref: http://dx.doi.org/10.1063/1.459170
	 * - use input: Coordinate 1: FFT grid, Coordinate 2: number basis
	 * - written in mass-weighted coordinates V is V(q), not V(x), where q=sqrt(m)*x
	 */

	/// H = T + |0><0|*V_11 + (|0><1|+|1><0|)*V12 + |1><1|*V_22
	Hamiltonian H;
	{ /// kinetic energy for FFT
		ProductOperator P;
		P.emplace_back(basis[0].kin_, 0);
		H.emplace_back(1., P);
	}

	size_t nstates = 2;
	{
		/// |0><0|*V_11
		ProductOperator P;
		P.emplace_back(element(nstates, 0, 0), 1);
		P.V_ = tullyA_diabatic11;
		H.emplace_back(1., P);
	}

	{
		/// (|0><1|+|1><0|)*V12
		Matrixcd proj = element(nstates, 0, 1) + element(nstates, 1, 0);
		ProductOperator P;
		P.emplace_back(proj, 1);
		P.V_ = tullyA_diabatic12;
		H.emplace_back(1., P);
	}

	{
		/// |1><1|*V_22
		ProductOperator P;
		P.emplace_back(element(nstates, 1, 1), 1);
		P.V_ = tullyA_diabatic22;
		H.emplace_back(1., P);
	}

	return H;
}
