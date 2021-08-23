//
// Created by Roman Ellerbrock on 8/19/21.
//

#include "PrimitiveOperators.h"
#include "Util/QMConstants.h"

/// Harmonic Oscillator
/**
 * Rationale:
 * - primitive operators for harmonic oscillator basis in FBR.
 * - x- and p- representation are obtained from creation-
 *   and annihilation-operators
 */

/// \brief x-operator representation for Harmonic Osciallator
Matrixcd x_HO(size_t dim, double freq) {
	Matrixcd x(dim, dim);
	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			if (m == n + 1) {
				x(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				x(m, n) = sqrt(n * 1.);
			}
			x(m, n) *= sqrt(1. / (2. * freq));
		}
	}
	return x;
}

/// \brief p-operator representation for Harmonic Osciallator
Matrixcd p_HO(size_t dim, double freq) {
	Matrixcd p(dim, dim);
	complex<double> imag(0., 1.);
	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			if (m == n + 1) {
				p(m, n) = sqrt(n + 1.);
			}
			if (m == n - 1) {
				p(m, n) = -sqrt(n * 1.);
			}
			p(m, n) *= imag * sqrt(freq / 2.);
		}
	}
	return p;
}

/// \brief T-operator representation for Harmonic Osciallator
Matrixcd kin_HO(size_t dim, double freq) {
	Matrixcd Kin(dim, dim);

	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i == j + 2) {
				Kin(j, i) = sqrt((j + 1.) * (j + 2.));
			}
			if (i == j) {
				Kin(j, i) = -(2 * j + 1);
			}
			if (i == j - 2) {
				Kin(j, i) = sqrt(j * (j - 1.));
			}
			// factor from p_-representation and from kinetic energy kin_=0.5*p_�
			Kin(j, i) *= (-freq / 2.) * (1 / 2.);
		}
	}
	return Kin;
}

/// \brief Build x-eigenvalues and FBR-DVR transformation matrix
pair<Matrixcd, Vectord> dvr_HO(size_t dim, double freq, double x0) {
	/// Build DVR from diagonalizing X matrix
	/// x = U * x_ev * U^dagger
	Matrixcd x = x_HO(dim, freq);
	auto spec = diagonalize(x);
	Vectord& grid = spec.second;

	for (int i = 0; i < grid.dim(); i++)
		grid(i) += x0;

	return spec;
}

/// FFT
/**
 * Rationale:
 * - This is a prinitive implementation, not optimized for performance.
 *   In an efficient implementation, one would keep x- and p- operator as
 *   a vector (i.e. as a grid in both x and p) and apply a fast Fourier Transform
 *   to transform the wavefunction from x- to momentum representation.
 *   Here, we build x- and p- as matrices. Hence, the FFT basis has the same
 *   implementation as the orthogonal basis functions.
 */

/// \brief build grid in x-space
Vectord xgrid_FFT(size_t dim, double x0, double x1) {
	auto N = (double) dim;
	double deltaX = (x1 - x0) / (N - 1);

	Vectord x(dim);
	for (size_t i = 0; i < dim; ++i) {
		x(i) = x0 + i * deltaX;
	}
	return x;
}

/// \brief build grid in p-space
Vectord pgrid_FFT(size_t dim, double x0, double x1) {
	auto N = (double) dim;
	double deltaX = (x1 - x0) / (N - 1);
	double deltaP = QM::two_pi / (deltaX * N);
	double lenP = (N - 1) * deltaP;
	double p0 = -lenP / 2.;

	Vectord p(dim);
	for (size_t i = 0; i < dim; ++i) {
		p(i) = p0 + i * deltaP;
	}
	return p;
}

/// \brief build dense x-operator in x-space (i.e. a full matrix)
Matrixcd x_FFT(size_t dim, double x0, double x1) {
	Vectord xgrid = xgrid_FFT(dim, x0, x1);

	Matrixcd x(dim, dim);
	for (int i = 0; i < dim; ++i) {
		x(i, i) = xgrid(i);
	}
	return x;
}

/// \brief build dense p-operator in p-space (i.e. a full matrix)
Matrixcd p_FFT(size_t dim, double x0, double x1) {
	/**
	 * Rationale:
	 * - For most use-cases the p-grid should be centered around 0.
	 *   deltaX determines width of the p-grid and (x1-x0) & dim determines
	 *   deltaP.
	 */
	Vectord pgrid = pgrid_FFT(dim, x0, x1);
	Matrixcd p(dim, dim);
	for (int i = 0; i < dim; ++i) {
		p(i, i) = pgrid(i);
	}
	return p;
}

/// \brief build dense T-operator in p-space (i.e. a full matrix)
Matrixcd kin_FFT(size_t dim, double x0, double x1) {
	Vectord pgrid = pgrid_FFT(dim, x0, x1);

	Matrixcd kin(dim, dim);
	for (int i = 0; i < dim; ++i) {
		kin(i, i) = 0.5 * pgrid(i) * pgrid(i);
	}
	return kin;
}

/// \brief build x-grid and transformation matrix
pair<Matrixcd, Vectord> dvr_FFT(size_t dim, double x0, double x1) {
	Vectord x = xgrid_FFT(dim, x0, x1);
	Vectord p = pgrid_FFT(dim, x0, x1);


	Matrixcd U(dim, dim);
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			U(j, i) = exp(-QM::im * (x(i) - x0) * p(j)) / sqrt((double) dim);
		}
	}
	return {U, x};
}

/// Number Basis
/**
 * Rationale:
 * - this is a logic basis set for creation- annihilation basis |n>.
 * - no representation for operator x and p.
 * - can be used for electronic state basis, spin basis, and many more.
 */
Matrixcd element(size_t dim, size_t n, size_t m) {
	/// Matrix-representation for |m><n|
	Matrixcd P(dim, dim);
	P(n, m) = 1.;
	return P;
}

/**
 * One could add Pauli x-, y-, z- operators, creation and annihilation operators.
 * This would allow solving Ising model Hamiltonians, Spin-Boson models, etc.
 * If you are interested, let me know. I can add content if someone is interested or
 * help you implementing it yourself.
 * Check out Wikipedia for more information on these systems, I can recommend implementing
 * transversal field Ising models in 1D and 2D and solving these both classically and
 * quantum mechanically.
 */
