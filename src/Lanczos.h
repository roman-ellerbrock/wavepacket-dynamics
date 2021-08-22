//
// Created by Roman Ellerbrock on 8/19/21.
//

#ifndef LANCZOS_H
#define LANCZOS_H
#include "Wavefunction.h"
#include "Hamiltonian.h"

class Lanczos {
public:
	Lanczos() = default;
	~Lanczos() = default;

	void calculate(Wavefunction Psi, const Hamiltonian& H,
		const Basis& basis, size_t krylov_size) {
		/// { |Psi>, H|Psi>, H^2|Psi>, ..., H^n|Psi> }

		Psi.occupy(basis);
		ofstream os("plot.dat");
		Psi.plot2D(basis, 0, 1, os);

		krylov_space = {Psi};
		vector<double> beta;
		vector<double> alpha;

		/// First step
		Wavefunction HPsi = H.apply(Psi);

		Matrixcd alphamat = Psi.dotProduct(HPsi);
		alpha.push_back(real(alphamat(0, 0)));

		Psi = HPsi - alpha.back() * Psi;

		/// remaining steps
		for (size_t l = 1; l < krylov_size+1; ++l) {
			beta.push_back(Psi.normalize());
			krylov_space.push_back(Psi);

			HPsi = H.apply(Psi);
			alphamat = Psi.dotProduct(HPsi);
			alpha.push_back(real(alphamat(0, 0)));

			const Wavefunction& v = krylov_space[l];
			const Wavefunction& vlast = krylov_space[l-1];
			Psi = HPsi - alpha.back() * v - beta.back() * vlast;
		}

		/// Build tri-diagonal hamiltonian matrix
		Hrep = Matrixcd(krylov_size, krylov_size);
		for (size_t l = 0; l < krylov_size; ++l) {
			Hrep(l, l) = alpha[l];
			if (l == krylov_size - 1) { break; }
			Hrep(l + 1, l) = beta[l];
			Hrep(l, l + 1) = beta[l];
		}
	}

	void timeStep(Wavefunction& Psi, double& t, double dt, const Hamiltonian& H,
		const Basis& basis, size_t krylov_size) {
		calculate(Psi, H, basis, krylov_size);

	}

	void print() const {
		cout << "Hamiltonian Representation:\n";
		Hrep.print();
		auto x = diagonalize(Hrep);
		cout << "Eigenvalues:\n";
		auto pres = cout.precision();
		cout << setprecision(8) << fixed;
		x.second.print();
		cout << setprecision(pres);
	}

	vector<Wavefunction> krylov_space;
	Matrixcd Hrep;
};

#endif //LANCZOS_H
