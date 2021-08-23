//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef PRODUCTOPERATOR_H
#define PRODUCTOPERATOR_H
#include "Core/Matrix.h"
#include "Wavefunction.h"

/**
 * \class Product Operator
 * \brief This is a product of one-dimensional operators, i.e. P = h_1 * h_2 * ... * h_size
 */
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

	/**
	 * \brief Translate idx (i1, i2, ...) into a grid point (x_i1, x_i2, ...)
	 * @param idx compoint index as vector of many indices
	 * @param basis
	 * @return grid point (x_i1, x_i2, ...)
	 */
	Vectord gridPoint(const vector<size_t> idx, const Basis& basis) const {
		Vectord x(basis.size());
		for (const PrimitiveBasis& prim : basis) {
			size_t k = prim.coord_;
			const Vectord& grid = prim.grid_;
			x(k) = grid(idx[k]);
		}
		return x;
	}

	/**
	 * \brief Apply PES using discrete variable representation (DVR)
	 * @param Psi Wavefunction (in and output)
	 * @param basis is the basis of the Psi
	 */
	void applyPotential(Wavefunction& Psi, const Basis& basis) const {
		if (!V_) { return; }
		const TensorShape& shape = Psi.shape();
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			auto idx = indexMapping(I, shape);
			auto x = gridPoint(idx, basis);
			Psi(I) *= V_(x);
		}
	}

	[[nodiscard]] Wavefunction apply(const Wavefunction& Psi, const Basis& basis) const {
		Wavefunction pPsi(Psi);
		/// Swipe through all operators (op) in product operator (this)
		for (const pair<Matrixcd, int>& op : (*this)) {
			const Matrixcd& m = op.first; /// matrix representation of operator
			const int k = op.second; /// coordinate the operator is applied to

			/// Check matrix Dimension to give error message
			if (m.dim1() != pPsi.shape()[k]) {
				cerr << "Error: matrix dimension does not fit tensors dimension.\n";
				exit(1);
			}

			/// Perform matrix-tensor product
			pPsi = matrixTensor(m, pPsi, k);
		}

		/// If a function pointer is stored, apply PES using DVR
		if (V_) {
			applyPotential(pPsi, basis);
		}
		return pPsi;
	}

	/**
	 * \brief print the matrix representaitons of all operators
	 */
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
