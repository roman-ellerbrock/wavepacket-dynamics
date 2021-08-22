//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef PRODUCTOPERATOR_H
#define PRODUCTOPERATOR_H
#include "Core/Matrix.h"
#include "Wavefunction.h"

class ProductOperator : public vector<pair<Matrixcd, int>>{
public:
	using vector::emplace_back;

	ProductOperator() = default;
	~ProductOperator() = default;

	ProductOperator(const Matrixcd& A, int coord) {
		emplace_back(A, coord);
	}

	void emplace_back(const Matrixcd& A, int coord) {
		pair<Matrixcd, int> op(A, coord);
		emplace_back(op);
	}

	[[nodiscard]] Wavefunction apply(const Wavefunction& Psi) const {
		Wavefunction pPsi(Psi);
		for (const pair<Matrixcd, int>& op : (*this)) {
			const Matrixcd& m = op.first;
			const int k = op.second;
			pPsi = matrixTensor(m, pPsi, k);
		}
		return pPsi;
	}

	void print() const {

		for (const pair<Matrixcd, int>& op : (*this)) {
			const Matrixcd& m = op.first;
			const int k = op.second;
			cout << "k=" << k << endl;
			m.print();
		}
	}

	/// Potential energy surface (replace by class to allow memory initialization)
	function<double(const Vectord& x)> V_;
};


#endif //PRODUCTOPERATOR_H
