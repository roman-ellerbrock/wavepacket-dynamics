//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef BASIS_H
#define BASIS_H
#include "yaml-cpp/yaml.h"
#include "PrimitiveBasis.h"
#include <vector>

/**
 * \class Basis
 * \brief This is the multidimensional basis of a wavefunction.
 */
class Basis : public std::vector<PrimitiveBasis>{
public:
	explicit Basis(const YAML::Node& node) {
		/// Save a primitive basis for every coordinate
		for (const YAML::Node& subnode: node) {
			emplace_back(PrimitiveBasis(subnode));
		}

		/// Create the Tensorshape shape = {dim1, dim2, ...} from all dimensions
		vector<size_t> dims;
		for (const auto& prim : *this) {
			dims.push_back(prim.dim_);
		}
		dims.push_back(1);
		shape_ = TensorShape(dims);
	}
	~Basis() = default;

	TensorShape shape_;
};


#endif //BASIS_H
