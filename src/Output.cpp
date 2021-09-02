//
// Created by Roman Ellerbrock on 8/25/21.
//

#include "Output.h"

/**
 * \brief Calculate the expectation value <Psi|P|Psi>
 * @param Psi
 * @param P
 * @param basis
 * @return expectation value <P>
 */
double expect(const Wavefunction& Psi, const ProductOperator& P, const Basis& basis) {
	Wavefunction PPsi = P.apply(Psi, basis);
	Matrixcd Prep = Psi.dotProduct(PPsi);
	return real(Prep(0, 0));
}

std::complex<double> autocorrelation(const Wavefunction &Psi, const Wavefunction &Psi0, const Basis &basis) {
  Matrixcd autoc = Psi.dotProduct(Psi0);
  return autoc(0, 0);
}

/**
 * \brief Calculate energy expectation value <Psi|H|Psi>
 * @param Psi
 * @param H
 * @param basis
 * @return expected energy E = <H>
 */
double energy_expect(const Wavefunction& Psi, const Hamiltonian& H, const Basis& basis) {
	Wavefunction HPsi = H.apply(Psi, basis);
	Matrixcd Hrep = Psi.dotProduct(HPsi);
	return real(Hrep(0, 0));
}

/**
 * Calculate expectation values of wavepacket and print them. Print wavefunction if ostream is given
 * @param Psi Wavefunction
 * @param H Hamiltonian of the system
 * @param basis Basis of the system
 * @param os_psi Output stream for plotting wavefunction
 */
void output(const Wavefunction& Psi, const Wavefunction &Psi0, const Hamiltonian& H, const Basis& basis, ostream* os_psi) {

	/// Calculate total Energy
	double norm = abs(Psi.dotProduct(Psi)(0, 0));
	double energy = energy_expect(Psi, H, basis);
  std::complex<double> ac = autocorrelation(Psi, Psi0, basis);
	cout << "<H> = " << energy << ", |Psi|^2 = " << norm << endl;
  cout << "<Psi(0)|Psi(t)> = " << real(ac) << " + " << imag(ac) << "i" << endl;
	/// for each coordinate, give output depending on the basis type
	for (const PrimitiveBasis& prim : basis) {
		if (prim.type_ == "HO" || prim.type_ == "FFT") {
			/// Calculate kinetic Energy
			ProductOperator T(prim.kin_, prim.coord_);
			double kin = expect(Psi, T, basis);
			ProductOperator X(prim.x_, prim.coord_);
			double x = expect(Psi, X, basis);
			ProductOperator P(prim.p_, prim.coord_);
			double p = expect(Psi, P, basis);
			X.emplace_back(prim.x_, prim.coord_);
			double x2 = expect(Psi, X, basis);
			double dx = sqrt(abs(x * x - x2));

			/// Make Output:
			cout << "coord: " << prim.coord_ << ", <T> = " << kin << ", <p> = " << p << ", <x> = " << x << ", <dx> = "
				 << dx << endl;
		} else if (prim.type_ == "NumberBasis") {
			cout << "coord: " << prim.coord_;
			double av = 0;
			for (size_t i = 0; i < prim.dim_; ++i) {
				ProductOperator proj(element(prim.dim_, i, i), prim.coord_);
				double el = expect(Psi, proj, basis);
				cout << ", <P_" << i << "> = " << el;
				av += i * el;
			}
			cout << ", <N> = " << av << endl;
		}
	}
	cout << endl;

	/// Plot wavefunction if ostream is provided
	if (os_psi) {
		ostream& os = *os_psi;
		Psi.plot(basis, {}, os);
	}
}

