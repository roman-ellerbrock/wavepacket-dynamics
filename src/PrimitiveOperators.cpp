//
// Created by Roman Ellerbrock on 8/19/21.
//

#include "PrimitiveOperators.h"

/// Harmonic Oscillator
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


/// Equidistant Grid
