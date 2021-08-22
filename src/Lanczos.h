//
// Created by Roman Ellerbrock on 8/19/21.
//

#ifndef LANCZOS_H
#define LANCZOS_H
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "Util/QMConstants.h"

class Lanczos {
public:
	Lanczos() = default;
	Lanczos(const Wavefunction& Psi, size_t krylov_size) {
		krylovSpace_ = vector<Wavefunction>(krylov_size, Wavefunction(Psi.shape_));
	}
	~Lanczos() = default;

	void calculate(Wavefunction Psi, const Hamiltonian& H,
		const Basis& basis, size_t krylov_size) {
		/// { |Psi>, H|Psi>, H^2|Psi>, ..., H^n|Psi> }

		krylovSpace_ = {Psi};
		vector<double> beta;
		vector<double> alpha;

		/// First step
		Wavefunction HPsi = H.apply(Psi);

		Matrixcd alphamat = Psi.dotProduct(HPsi);
		alpha.push_back(real(alphamat(0, 0)));

		Psi = HPsi - alpha.back() * Psi;

		/// remaining steps
		for (size_t l = 1; l < krylov_size; ++l) {
			beta.push_back(Psi.normalize());
			krylovSpace_.push_back(Psi);

			HPsi = H.apply(Psi);
			alphamat = Psi.dotProduct(HPsi);
			alpha.push_back(real(alphamat(0, 0)));

			const Wavefunction& v = krylovSpace_[l];
			const Wavefunction& vlast = krylovSpace_[l - 1];
			Psi = HPsi - alpha.back() * v - beta.back() * vlast;
		}

		/// Build tri-diagonal hamiltonian matrix
		Hrep_ = Matrixcd(krylov_size, krylov_size);
		for (size_t l = 0; l < krylov_size; ++l) {
			Hrep_(l, l) = alpha[l];
			if (l == krylov_size - 1) { break; }
			Hrep_(l + 1, l) = beta[l];
			Hrep_(l, l + 1) = beta[l];
		}
		auto spectrum = diagonalize(Hrep_);
		U_ = spectrum.first;
		E_ = spectrum.second;
	}

	void timeStep(Wavefunction& Psi, double& t, double dt) {

		/// Calculate propagation coefficients
		size_t krylov_size = krylovSpace_.size();
		if (krylov_size == 0) {
			cerr << "Error: no krylov vectors. Calculate space first.\n";
			exit(1);
		}
		Vectorcd c(krylov_size);
		for (size_t l = 0; l < krylov_size; ++l) {
			for (size_t k = 0; k < krylov_size; ++k) {
				complex<double> exp_k = exp(-QM::im * E_(k) * dt);
				c(l) += U_(l, k) * exp_k * conj(U_(0, k));
			}
		}

		const TensorShape& shape = Psi.shape();
		Psi.zero();
		for (size_t l = 0; l < krylov_size; ++l) {
			const Wavefunction& psi_l = krylovSpace_[l];
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				Psi(I) += c(l) * psi_l(I);
			}
		}
		t += dt;
	}

	void integrateNext(Wavefunction& Psi, double& t, double& dt, double tnext, double tmax,
		const Hamiltonian& H, const Basis& basis, size_t krylov_size) {
		while (t + 1e-10 < tnext) {
			dt = min(dt, tnext - t);
			/// re-calculate Krylov space
			calculate(Psi, H, basis, krylov_size);
			/// make timestep
			timeStep(Psi, t, dt);

			/// minimal output
			auto pres = cout.precision();
			cout << setprecision(2) << fixed;
			cout << "time: " << t << " / " << tmax << " |Psi|^2  = " << Psi.normalize() << endl;
			cout << setprecision(pres);
		}
	}

	void integrate(Wavefunction& Psi, double& t, double& dt, double tmax,
		double out, const Hamiltonian& H, const Basis& basis, size_t krylov_size) {
		Wavefunction Psi_last = Psi;
		bool movie = false;

		size_t count = 0;
		while (t + 1e-6 < tmax) {
			double tnext = min(out, tmax - t) + t;
			integrateNext(Psi, t, dt, tnext, tmax, H, basis, krylov_size);

			if (movie) {
				ofstream os("Psi." + to_string(count++) + ".dat");
				os << "# time: " << t << "/" << tmax << endl;
				Psi.plot2D(basis, 0, 1, os);
			}
		}
	}

	void print() const {
		cout << "Hamiltonian Representation:\n";
		Hrep_.print();
		cout << "Eigenvalues:\n";
		auto pres = cout.precision();
		cout << setprecision(8) << fixed;
		E_.print();
		cout << setprecision(pres);
	}

	vector<Wavefunction> krylovSpace_;
	Matrixcd Hrep_;
	Vectord E_;
	Matrixcd U_;
};

#endif //LANCZOS_H
