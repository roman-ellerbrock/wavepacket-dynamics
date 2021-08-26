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

/**
 * Tully Model systems
 * Rationale:
 * - Tully Model Potential from Sec. IV A in 1990 paper
 * - The wave packet code does not explicitely treat mass, instead everything
 *   is expressed in terms of mass-weighted coordiantes
 */

/**
 * \brief transform from massweighted coordinate q = sqrt(m) * x to x
 * @param q mass weighted coordinate
 * @return x (mass unweighted coordinate)
 */
double unweight(const double q) {
	double m = 1000.;
	double x = q / sqrt(m);
	return x;
}

double V_tullyA_diabatic11(const double q) {
	auto x = unweight(q);
	double A = 0.01;
	double B = 1.6;

	double V = 0;
	if (x >= 0) {
		V = A * (1. - exp(-B * x));
	} else {
		V = -A * (1. - exp(B * x));
	}
	return V;
}

double V_tullyA_diabatic12(const double q) {
	auto x = unweight(q);
	double C = 0.005;
	double D = 1.0;

	double V = C * exp(-D*x*x);
	return V;
}

double V_tullyA_diabatic22(const double q) {
	auto x = unweight(q);
	return -V_tullyA_diabatic11(x);
}

/**
 * \brief build 2x2 matrix of V in Tully Model A
 * @return V_diabatic
 */
Matrixcd tullyA_diabaticMatrix(const double q) {
	///
	Matrixcd Vdiab(2, 2);
	Vdiab(0, 0) = V_tullyA_diabatic11(q);
	Vdiab(1, 0) = V_tullyA_diabatic12(q);
	Vdiab(0, 1) = conj(Vdiab(1, 0));
	Vdiab(1, 1) = V_tullyA_diabatic22(q);
	return Vdiab;
}

/**
 * @return transformation matrix that transforms diabatic rep to adiabatic at x=-\infty
 */
Matrixcd U_left() {
	/// x at -infinity
	double q = -350;
	Matrixcd Vdiab = tullyA_diabaticMatrix(q);
	auto eigen = diagonalize(Vdiab);
	return eigen.first;
}

Matrixcd tullyA_diabatixMatrixTransformed(const double q) {
	Matrixcd Vdiab = tullyA_diabaticMatrix(q);
	Matrixcd U = U_left();
	return unitarySimilarityTrafo(Vdiab, U);
}

double tullyA_diabatic11(const Vectord& q) {
	Matrixcd Vdiab = tullyA_diabatixMatrixTransformed(q(0));
	return real(Vdiab(0, 0));
}

double tullyA_diabatic12(const Vectord& q) {
	Matrixcd Vdiab = tullyA_diabatixMatrixTransformed(q(0));
	return real(Vdiab(0, 1));
}

double tullyA_diabatic22(const Vectord& q) {
	Matrixcd Vdiab = tullyA_diabatixMatrixTransformed(q(0));
	return real(Vdiab(1, 1));
}

/**
 * /// Tully model A, adiabatic V-elements
 */

Vector<double> V_tullyA_adiabatic(const Vectord& q) {
	Matrixcd V(2, 2);
	V(0, 0) = V_tullyA_diabatic11(q(0));
	V(1, 0) = V_tullyA_diabatic12(q(0));
	V(0, 1) = V_tullyA_diabatic12(q(0));
	V(1, 1) = V_tullyA_diabatic22(q(0));
	auto eigen = diagonalize(V);
	return eigen.second;
}

double V_tullyA_adiabatic1(const Vectord& q) {
	Vectord ev = V_tullyA_adiabatic(q);
	return ev(0);
}

double V_tullyA_adiabatic2(const Vectord& q) {
	Vectord ev = V_tullyA_adiabatic(q);
	return ev(1);
}
