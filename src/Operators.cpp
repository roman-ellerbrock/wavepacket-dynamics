//
// Created by Roman Ellerbrock on 8/19/21.
//

#include "Operators.h"

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
			double coeff = 0.5 * prim.freq_ * prim.freq_;
			H.emplace_back(coeff, P);
		}

		/// quartic Potential
		{
			ProductOperator P;
			auto x2 = prim.x_ * prim.x_;
			auto x4 = x2 * x2;
			P.emplace_back(x4 , coord);
			double coeff = 0.05;
			H.emplace_back(coeff, P);
		}
	}
	return H;
}
