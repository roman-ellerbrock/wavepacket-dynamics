//
// Created by Roman Ellerbrock on 8/19/21.
//

#ifndef LANCZOS_H
#define LANCZOS_H
#include "Wavefunction.h"
#include "Hamiltonian.h"
#include "Util/QMConstants.h"

/**
 * \class Lanczos
 * \brief Build Krylov space {Psi, HPsi, H^2Psi, ...} and evaluate eigenstates or integrate in time using this space.
 */
class Lanczos {
public:
	Lanczos() = default;
	Lanczos(const Wavefunction& Psi, size_t krylov_size) {
		initialize(Psi, krylov_size);
	}
	~Lanczos() = default;
	void initialize(const Wavefunction& Psi, size_t krylov_size) {
		krylovSpace_ = vector<Wavefunction>(krylov_size, Wavefunction(Psi.shape()));
	}

	/**
	 * \brief Build Krylov Space K={Psi, HPsi, H^2Psi, ...} using Lanczos algorithm.
	 * @param Psi Wavefunction
	 * @param H Hamiltonian
	 * @param basis multidimensional Basis
	 * @param krylov_size size of the Krylov space
	 *
	 * The Lanczos algorithm is a krylov method obtained from a three point recursion.
	 */
	void calculate(Wavefunction Psi, const Hamiltonian& H,
		const Basis& basis, size_t krylov_size) {
		/// { |Psi>, H|Psi>, H^2|Psi>, ..., H^n|Psi> }

		if (krylovSpace_.size() != krylov_size) { initialize(Psi, krylov_size); }

		krylovSpace_[0] = Psi;
		vector<double> beta;
		vector<double> alpha;

		/// First step
		Wavefunction HPsi = H.apply(Psi, basis);

		Matrixcd alphamat = Psi.dotProduct(HPsi);
		alpha.push_back(real(alphamat(0, 0)));

		Psi = HPsi - alpha.back() * Psi;

		/// remaining steps
		for (size_t l = 1; l < krylov_size; ++l) {
			beta.push_back(Psi.normalize());
			krylovSpace_[l] = Psi;

			HPsi = H.apply(Psi, basis);
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
		/// Diagonalize Hamiltonian matrix and save transformation matrix U and eigenvalues E.
		auto spectrum = diagonalize(Hrep_);
		U_ = spectrum.first;
		E_ = spectrum.second;
	}


	/**
	 * \brief Perform Short iterative Lanczos (SIL) integration to solve time-dependent Schrödinger Equation
	 * @param Psi initial Wavefunction Psi(t), will be Psi(tmax) at out
	 * @param t current time
	 * @param dt initial timestep
	 * @param tmax end of integration
	 * @param out output at t = i * out
	 * @param H Hamiltonian for integration -i d/dt Psi = H Psi
	 * @param basis Basis for wavefunciton
	 * @param krylov_size size of krylov space
	 *
	 * Set movie to true to plot wavefunction frames at every output; allows to make gifs
	 */
	void integrate(Wavefunction& Psi, double& t, double& dt, double tmax,
		double out, const Hamiltonian& H, const Basis& basis, size_t krylov_size) {
		bool movie = true;

		size_t count = 0;
		if (movie) {
			ofstream os("tmp/Psi." + to_string(count++) + ".dat");
			os << "# time: " << t << "/" << tmax << endl;
			Psi.plot(basis, {}, os);
		}
		while (t + 1e-6 < tmax) {
			double tnext = min(out, tmax - t) + t;
			integrateNext(Psi, t, dt, tnext, tmax, H, basis, krylov_size);

			if (movie) {
				ofstream os("tmp/Psi." + to_string(count++) + ".dat");
				os << "# time: " << t << "/" << tmax << endl;
				Psi.plot(basis, {}, os);
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

private:
	/**
	 * \brief Perform a step with a short-iterative Lanczos algorithm
	 * @param Psi Wavefunction
	 * @param t current time
	 * @param dt time step
	 */
	void timeStep(Wavefunction& Psi, double& t, double dt) {

		/// Calculate propagation coefficients
		size_t krylov_size = krylovSpace_.size();
		if (krylov_size == 0) {
			cerr << "Error: no krylov vectors. Calculate space first.\n";
			exit(1);
		}
		/// Calculate coefficients of Psi(t+dt) in krylov basis, |Psi(t+dt)>= \sum_l c_l * |psi_l>
		Vectorcd c(krylov_size);
		for (size_t l = 0; l < krylov_size; ++l) {
			for (size_t k = 0; k < krylov_size; ++k) {
				complex<double> exp_k = exp(-QM::im * E_(k) * dt);
				c(l) += U_(l, k) * exp_k * conj(U_(0, k));
			}
		}

		/// Calculate |Psi(t+dt)>= \sum_l c_l * |psi_l> using pre-calculated coefficients and krylov vectors
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

	/**
	 * Integrate in time up to next output.
	 * @param t current time
	 * @param dt time-step
	 * @param tnext time to next output or end of integration
	 * @param tmax maximal integration time
	 */
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
};

#endif //LANCZOS_H
