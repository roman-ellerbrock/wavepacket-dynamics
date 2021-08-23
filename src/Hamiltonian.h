//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include "ProductOperator.h"
#include <vector>

using namespace std;

/**
 * \class Hamiltonian
 * \brief This is a Hamiltonian with sum-of-product structure. ProductOperators can contain a PES.
 */
class Hamiltonian : public vector<pair<double, ProductOperator>> {
public:
	using vector::emplace_back;

	Hamiltonian() = default;
	~Hamiltonian() = default;

	void emplace_back(double coeff, const ProductOperator& P) {
		pair<double, ProductOperator> po({coeff, P});
		emplace_back(po);
	}

	/**
	 * \brief apply Hamiltonian to a wavefunction
	 * @param Psi The wavefunction
	 * @param basis The basis corresponding to Psi
	 * @return HPsi
	 */
	Wavefunction apply(const Wavefunction& Psi, const Basis& basis) const{
		Wavefunction HPsi(Psi.shape());
		/// Loop over all summands in the wavefunction and apply them
		for (const pair<double, ProductOperator>& pa : *this) {
			const double c = pa.first;
			const ProductOperator& P = pa.second;

			/// Apply product operator and multiply with coefficient
			HPsi += c * P.apply(Psi, basis);
		}
		return HPsi;
	}

};


#endif //HAMILTONIAN_H
