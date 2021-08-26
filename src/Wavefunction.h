//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "Core/Tensor.h"
#include "Basis.h"

/**
 * \brief This class represents the direct product wavefunction. It is a decorated tensor.
 */
class Wavefunction: public Tensorcd {
public:

	Wavefunction() = default;

	/// \brief Create wavefunction and allocate its memory. Elements are set to zero
	explicit Wavefunction(const Basis& basis)
		: Tensorcd(basis.shape_) {}

	/// \brief Occupy the wavefunction with direct product of 1D wavefunctions from Primitive basis.
	void occupy(const Basis& basis) {
		Tensorcd& A = *this;
		const TensorShape& shape = A.shape();

		/// Create even superposition
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			A[I] = 1.;
		}

		for (const auto& prim : basis) {
			/// Get 1D wavefunction and create a product
			auto psi1D = prim.initial1DWavefunction();
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				auto idx = indexMapping(I, shape);
				A[I] *= psi1D(idx[prim.coord_]);
			}
		}
		normalize();
	}

	/**
	 * \brief Plot wavefunction without dynamic dimensionality.
	 * @param basis Basis of the wavefunction
	 * @param coords Coordinates that shall be plottet (up to two). If empty, it will plot the first one(s).
	 * @param os stream that the output is written to
	 */
	void plot(const Basis& basis, const vector<size_t>& coords, ostream& os = cout) const {
		if (coords.empty()) {
			if (basis.empty()) {
				cerr << "Error: empty basis!\n";
				exit(1);
			}
			vector<size_t> cs;
			for (const PrimitiveBasis& prim : basis) {
				cs.emplace_back(prim.coord_);
			}
			plot(basis, cs, os);
			return;
		}
		if (coords.size() == 1) {
			plot1D(basis, coords.front(), os);
		} else {
			plot2D(basis, coords[0], coords[1], os);
		}
	}

	/// \brief: Transform from finite basis rep. (FBR) to discrete variable rep. (DVR)
	void fbrToDVR(const Basis& basis) {
		Tensorcd& A = *this;
		for (const auto& coord : basis) {
			A = matrixTensor(coord.trafo_, A, coord.coord_);
		}
	}

	/// \brief: Transform from discrete variable rep. (DVR) to finite basis rep. (FBR)
	void dvrToFBR(const Basis& basis) {
		Tensorcd& A = *this;
		for (const auto& coord : basis) {
			A = matrixTensor(coord.trafo_.adjoint(), A, coord.coord_);
		}
	}

	/// \brief normalize the wavefunction and return norm.
	double normalize() {
		Tensorcd& A = *this;
		double norm = A.dotProduct(A).frobeniusNorm();
		norm = sqrt(norm);
		A /= norm;
		return norm;
	}

	/// Fancy template stuff, just necessary C++ details... ==>
	using Tensorcd::Tensor;
	using Tensorcd::operator=;
	using Tensorcd::operator*=;

	Wavefunction(const Tensorcd& A)
		: Tensorcd(A) {}
	/// <===
private:

	/// \brief Plot wavefunction in 1D
	void plot1D(const Basis& basis, size_t k1 = 0, ostream& os = cout) const {
		const PrimitiveBasis& prim = basis.front();
		const Tensorcd& A = *this;
		for (size_t i = 0; i < prim.dim_; ++i) {
			os << prim.grid_(i) << " " << pow(abs(A[i]), 2) << endl;
		}
	}

	void plot2D(const Basis& basis, size_t k1, size_t k2, ostream& os = cout) const {
		/**
		 * Rationale:
		 * - Plot wavefunction after integrating out all but two coordinates.
		 */
		/// Calculate 2D density (integrate out all other coords)
		size_t dim1 = basis[k1].dim_;
		size_t dim2 = basis[k2].dim_;
		Matrixd d(dim1, dim2);
		const Tensorcd& A = *this;
		for (size_t I = 0; I < basis.shape_.totalDimension(); ++I) {
			auto idx = indexMapping(I, basis.shape_);
			d(idx[k1], idx[k2]) += pow(abs(A[I]), 2);
		}

		/// Plot 2D density
		const Vectord& x1 = basis[k1].grid_;
		const Vectord& x2 = basis[k2].grid_;
		for (size_t i1 = 0; i1 < dim1; ++i1) {
			for (size_t i2 = 0; i2 < dim2; ++i2) {
				os << x1[i1] << " " << x2[i2] << " " << d(i1, i2) << endl;
			}
			os << endl;
		}
	}

};


#endif //WAVEFUNCTION_H
