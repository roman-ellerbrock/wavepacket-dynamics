//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include "ProductOperator.h"
#include <vector>

using namespace std;

class Hamiltonian : public vector<pair<double, ProductOperator>> {
public:
	using vector::emplace_back;

	Hamiltonian() = default;
	~Hamiltonian() = default;

	void emplace_back(double coeff, const ProductOperator& P) {
		pair<double, ProductOperator> po({coeff, P});
		emplace_back(po);
	}

	Wavefunction apply(const Wavefunction& Psi) const{
		Wavefunction HPsi(Psi.shape_);
		/// Loop over all summands in the wavefunction and apply them
		for (const pair<double, ProductOperator>& pa : *this) {
			const double c = pa.first;
			const ProductOperator& P = pa.second;

			HPsi += c * P.apply(Psi);
		}
		return HPsi;
	}

};


#endif //HAMILTONIAN_H
