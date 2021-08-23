//
// Created by Roman Ellerbrock on 8/22/21.
//
#include "Potentials.h"

double V_HO(const Vectord& x) {
	double V = 0.;
	for (size_t k = 0; k < x.dim(); ++k) {
		V += 0.5 * x(k) * x(k);
	}
	return V;
}

double V_Eckhard(const Vectord& x) {
	double omega = 1.;
	double D = 0.1;
	return D * exp(-omega*(x(0)*x(0)));
}

/// Tully Model Potential from Sec. IV A in 1990 paper

double V_tullyA_diabatic11(const Vectord& x) {
	double A = 0.01;
	double B = 1.6;

	double V = 0;
	if (x(0) >= 0) {
		V = A * (1. - exp(-B * x(0)));
	} else {
		V = -A * (1. - exp(-B * x(0)));
	}
	return V;
}

double V_tullyA_diabatic12(const Vectord& x) {
	double C = 0.005;
	double D = 1.0;

	double V = C * exp(-D*x(0)*x(0));
	return V;
}

double V_tullyA_diabatic22(const Vectord& x) {
	return -V_tullyA_diabatic11(x);
}
