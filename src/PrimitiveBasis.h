//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef PRIMITIVEBASIS_H
#define PRIMITIVEBASIS_H
#include "yaml-helper.h"
#include "yaml-cpp/yaml.h"
#include "Core/Tensor.h"
#include "Core/Matrix.h"

class PrimitiveBasis {
public:
	explicit PrimitiveBasis(const YAML::Node& node) {
		dim_ = read_key<int>(node, "dim");
		coord_ = read_key<int>(node, "coord");
		type_ = read_key<string>(node, "type");
		freq_ = read_key<double>(node, "freq");
		x0_ = read_key<double>(node, "x0");
		wffreq_ = read_key<double>(node, "wffreq");
		wfx0_ = read_key<double>(node, "wfx0");

		buildX();
		buildP();
		buildKin();

		/// perform DVR last!
		buildDVR();
	}

	~PrimitiveBasis() = default;

	void buildX() {
		Matrixcd x(dim_, dim_);
		for (int n = 0; n < dim_; n++) {
			for (int m = 0; m < dim_; m++) {
				if (m == n + 1) {
					x(m, n) = sqrt(n + 1.);
				}
				if (m == n - 1) {
					x(m, n) = sqrt(n * 1.);
				}
				x(m, n) *= sqrt(1. / (2. * freq_));
			}
		}
		x_ = x;
	}

	void buildP() {
		Matrixcd p(dim_, dim_);
		complex<double> imag(0., 1.);
		for (int n = 0; n < dim_; n++) {
			for (int m = 0; m < dim_; m++) {
				if (m == n + 1) {
					p(m, n) = sqrt(n + 1.);
				}
				if (m == n - 1) {
					p(m, n) = -sqrt(n * 1.);
				}
				p(m, n) *= imag * sqrt(freq_ / 2.);
			}
		}
		p_ = p;
	}

	void buildKin() {
		Matrixcd Kin(dim_, dim_);

		for (int i = 0; i < dim_; i++) {
			for (int j = 0; j < dim_; j++) {
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
				Kin(j, i) *= (-freq_ / 2.) * (1 / 2.);
			}
		}
		kin_ = Kin;
	}

	void buildDVR() {
		/// Build DVR from diagonalizing X matrix
		/// x = U * x_ev * U^dagger
		auto spec = diagonalize(x_);
		grid_ = spec.second;
		trafo_ = spec.first;

		for (int i = 0; i < dim_; i++)
			grid_(i) += x0_;

		kin_ = unitarySimilarityTrafo(kin_, trafo_);
		p_   = unitarySimilarityTrafo(p_, trafo_);
		x_   = unitarySimilarityTrafo(x_, trafo_);
	}

	Matrixcd x_;
	Matrixcd p_;
	Matrixcd kin_;

	/// for DVR
	Vectord grid_;
	Matrixcd trafo_;

	/// basis parameters
	string type_; /// HO, FFT, Legendre, ...
	int dim_;
	int coord_;
	double freq_;
	double x0_;
	/// Initial wavepacket parameters
	double wffreq_;
	double wfx0_;
};


#endif //PRIMITIVEBASIS_H
