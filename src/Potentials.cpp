//
// Created by Roman Ellerbrock on 8/22/21.
//
#include "Potentials.h"
#include "Core/Matrix.h"

double V_HO(const Vectord& x) {
	double V = 0.;
	for (size_t k = 0; k < x.dim(); ++k) {
		V += 0.5 * x(k) * x(k);
	}
	return V;
}

double V_Eckhard(const Vectord& x) {
	double omega = 1.;
	double D = 1.;
	return D * exp(-omega*(x(0)*x(0)));
}

/// Tully Model Potential from Sec. IV A in 1990 paper

Vectord unweight(const Vectord& q) {
	double m = 1000.;
	Vectord x(1);
	x(0) = q(0) / sqrt(m);
	return x;
}

double V_tullyA_diabatic11(const Vectord& x) {
	double A = 0.01;
	double B = 1.6;

	double V = 0;
	if (x(0) >= 0) {
		V = A * (1. - exp(-B * x(0)));
	} else {
		V = -A * (1. - exp(B * x(0)));
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

double V_tullyA_adiabatic1(const Vectord& q) {
	auto x = unweight(q);
	Matrixcd V(2, 2);
	V(0, 0) = V_tullyA_diabatic11(x);
	V(1, 0) = V_tullyA_diabatic12(x);
	V(0, 1) = V_tullyA_diabatic12(x);
	V(1, 1) = V_tullyA_diabatic22(x);
	auto eigen = diagonalize(V);
	Vectord ev = eigen.second;
	return ev(0);
}
