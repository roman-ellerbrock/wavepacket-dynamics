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
Matrixcd x_FFT(size_t dim, double x0, double x1) {
	auto N = (double) dim;
	double deltaX = (x1 - x0) / (N - 1);

	Matrixcd x(dim, dim);
	for (int i = 0; i < dim; ++i) {
		x(i, i) = x0 + i * deltaX;
	}
	return x;
}

Matrixcd p_FFT(size_t dim, double x0, double x1) {
	/**
	 * Rationale:
	 * - For most use-cases the p-grid should be centered around 0.
	 *   deltaX determines width of the p-grid and (x1-x0) & dim determines
	 *   deltaP.
	 */
	auto N = (double) dim;
	double deltaX = (x1 - x0) / (N - 1);
	double lenP = QM::two_pi / deltaX;
	double p0 = - 0.5 * lenP;
	double deltaP = lenP / (N - 1);

	Matrixcd p(dim, dim);
	for (int i = 0; i < dim; ++i) {
		p(i, i) = p0 + i * deltaP;
	}
	p *= QM::im;
	return p;
}

Matrixcd kin_FFT(size_t dim, double x0, double x1) {
	auto N = (double) dim;
	double deltaX = (x1 - x0) / (N - 1);
	double lenP = QM::two_pi / deltaX;
	double p0 = - 0.5 * lenP;
	double deltaP = lenP / (N - 1);


	Matrixcd kin(dim, dim);
	for (int i = 0; i < dim; ++i) {
		double p = p0 + i * deltaP;
		kin(i, i) = - 0.5 * p * p;
	}
	return kin;
}

pair<Matrixcd, Vectord> dvr_FFT(size_t dim, double x0, double x1) {
	auto N = (double) dim;
	double deltaX = (x1 - x0) / (N - 1);
	double lenP = QM::two_pi / deltaX;
	double p0 = - 0.5 * lenP;
	double deltaP = lenP / (N - 1);

	Matrixcd U(dim, dim);
	Vectord grid(dim);
	for (int i = 0; i < dim; ++i) {
		double x = x0 + i * deltaX;
		grid(i) = x;
		for (int j = 0; j < dim; ++j) {
			U(j, i) = exp(-QM::im * QM::two_pi * (double) i * (double) j / N) / sqrt(N);
		}
	}
	return {U, grid};
}

/// Equidistant Grid
