//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef PRIMITIVEBASIS_H
#define PRIMITIVEBASIS_H
#include "yaml-helper.h"
#include "yaml-cpp/yaml.h"
#include "Core/Tensor.h"
#include "Core/Matrix.h"
#include "PrimitiveOperators.h"

class PrimitiveBasis {
public:
	explicit PrimitiveBasis(const YAML::Node& node) {
		dim_ = read_key<int>(node, "dim");
		coord_ = read_key<int>(node, "coord");
		type_ = read_key<string>(node, "type");

		if (type_ == "HO") {
			freq_ = read_key<double>(node, "freq");
			x0_ = read_key<double>(node, "x0");
			wffreq_ = read_key<double>(node, "wffreq");
			wfx0_ = read_key<double>(node, "wfx0");
		} else if (type_ == "NumberBasis") {
			wffreq_ = read_key<double>(node, "theta");
		} else if (type_ == "FFT") {
			x0_ = read_key<double>(node, "x0");
			freq_ = read_key<double>(node, "x1");
			wffreq_ = read_key<double>(node, "wffreq");
			wfx0_ = read_key<double>(node, "wfx0");
		} else {
			cerr << "Error: unknown basis type. Type was: " << type_ << endl;
			exit(3);
		}

		initializeOperators();

		/// perform DVR last!
		buildDVR();
	}

	~PrimitiveBasis() = default;

	void initializeOperators() {
		if (type_ == "HO") {
			/// build operators in FBR
			kin_ = kin_HO(dim_, freq_);
			x_ = x_HO(dim_, freq_);
			p_ = p_HO(dim_, freq_);

			/// calculate DVR
			auto spec = dvr_HO(dim_, freq_, x0_);
			trafo_ = spec.first;
			grid_ = spec.second;

			/// transform operators to DVR
			kin_ = unitarySimilarityTrafo(kin_, trafo_);
			p_   = unitarySimilarityTrafo(p_, trafo_);
			x_   = unitarySimilarityTrafo(x_, trafo_);

		} else if (type_ == "FFT") {
			cerr << "Not implemented, yet.\n";
			exit(1);
		} else if (type_ == "NumberBasis") {
			cerr << "Not implemented, yet.\n";
			exit(1);
		} else {
			cerr << "Error: unknown basis type. Type was: " << type_ << endl;
			exit(3);
		}
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
