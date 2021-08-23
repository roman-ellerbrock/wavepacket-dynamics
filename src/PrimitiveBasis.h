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
#include "Util/QMConstants.h"

/**
 * \class PrimitiveBasis
 * \brief This manages the 1D basis information. It creates operator representations and initial wave packets.
 */
class PrimitiveBasis {
public:
	/**
	 * \brief Initialize the parameters in the Primitive basis from a yaml-node.
	 * @param node
	 * Note: not all variables are actively used by all basis types. It would be
	 *       safer to inherit from this class to specify the type of basis.
	 *       Since this would result in various classes, I avoided the inheritance here.
	 */
	explicit PrimitiveBasis(const YAML::Node& node) {
		dim_ = read_key<int>(node, "dim");
		coord_ = read_key<int>(node, "coord");
		type_ = read_key<string>(node, "type");

		if (type_ == "HO") {
			freq_ = read_key<double>(node, "freq");
			x0_ = read_key<double>(node, "x0");
			wffreq_ = read_key<double>(node, "wffreq");
			wfx0_ = read_key<double>(node, "wfx0");
			wfp0_ = read_key<double>(node, "wfp0");

		} else if (type_ == "NumberBasis") {
			wffreq_ = read_key<double>(node, "theta");

		} else if (type_ == "FFT") {
			x0_ = read_key<double>(node, "x0");
			x1_ = read_key<double>(node, "x1");
			wffreq_ = read_key<double>(node, "wffreq");
			wfx0_ = read_key<double>(node, "wfx0");
			wfp0_ = read_key<double>(node, "wfp0");

		} else {
			cerr << "Error: unknown basis type. Type was: " << type_ << endl;
			exit(3);
		}

		initializeOperators();
	}

	~PrimitiveBasis() = default;

	/**
	 * \brief Initialize the operators that are attached to the basis and
	 * transform every operator to grid representation (DVR) if it exists.
	 */
	void initializeOperators() {
		if (type_ == "HO") {
			/// build operators in FBR
			kin_ = kin_HO(dim_, freq_);
			x_ = x_HO(dim_, freq_);
			p_ = p_HO(dim_, freq_);

			/// Your code comes here
			/// a.) build trafo_ and grid_ using dvr_HO ind PrimitiveOperators.cpp
			/// calculate DVR
			auto spec = dvr_HO(dim_, freq_, x0_);
			trafo_ = spec.first;
			grid_ = spec.second;

			/// b.) transform x, p and kin to DVR
			/// transform operators to DVR

			/// c.) build the Hamiltonian and end program

			exit(0);
		} else if (type_ == "FFT") {
			/// Build x, p, kin
			x_ = x_FFT(dim_, x0_, x1_);
			p_ = p_FFT(dim_, x0_, x1_);
			kin_ = kin_FFT(dim_, x0_, x1_);

			/// get the transformation matrix and x-grid
			auto spec = dvr_FFT(dim_, x0_, x1_);
			trafo_ = spec.first;
			grid_ = spec.second;

			/// transform kin and p; x is already in grid rep.
			kin_ = unitarySimilarityTrafo(kin_, trafo_);
			p_ = unitarySimilarityTrafo(p_, trafo_);

			/// freq isn't used for FFT. Set to high value to see errors early. Argh.. ^_^
			freq_ = -1000.;
		} else if (type_ == "NumberBasis") {
			/// Initialize a virtual grid that allows plotting the wavefunction (x is number operator)
			grid_ = Vectord(dim_);
			for (size_t i = 0; i < dim_; ++i) {
				grid_(i) = i;
			}
		} else {
			cerr << "Error: unknown basis type. Type was: " << type_ << endl;
			exit(3);
		}
	}

	/**
	 * \brief Create a 1D initial wavefunction. Is used to create multidimensional wavefunction
	 * @return Vector(dim) with the 1D wavefunction in DVR (HO, FFT) or in FBR (Number basis)
	 */
	[[nodiscard]] Vectorcd initial1DWavefunction() const {
		Vectorcd psi(dim_);
		if (type_ == "NumberBasis") {
			if (dim_ == 1) {
				psi(0) = 1.;
			} else {
				psi(0) = cos(wffreq_);
				psi(1) = sin(wffreq_);
			}
		} else {
			for (size_t i = 0; i < dim_; ++i) {
				/// gauss wavepacket with momentum wfp0 located at wfx0
				psi(i) = exp(QM::im * wfp0_*grid_(i)) * exp(-0.5 * wffreq_ * (pow(grid_(i) - wfx0_, 2)));
			}
			normalize(psi);
		}
		return psi;
	}

	/// Matrix representations of operators
	Matrixcd x_;
	Matrixcd p_;
	Matrixcd kin_;

	/// for DVR
	Vectord grid_;
	Matrixcd trafo_;

	/// basis type
	string type_; /// HO, FFT, Legendre, ...

	/// basis parameters
	/// Note: Here inheritance for different types would be cleaner but I want to avoid creating many classes.
	int dim_;
	int coord_;
	double freq_; /// freq has different meanings depending on type
	double x0_;
	double x1_; /// only exists for FFT

	/// Initial wavepacket parameters
	double wffreq_;
	double wfx0_;
	double wfp0_;
};


#endif //PRIMITIVEBASIS_H
