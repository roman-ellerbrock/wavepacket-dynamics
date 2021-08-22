//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef BASIS_H
#define BASIS_H
#include "yaml-cpp/yaml.h"
#include "PrimitiveBasis.h"
#include <vector>

class Basis : public std::vector<PrimitiveBasis>{
public:
	explicit Basis(const YAML::Node& node) {
		for (const YAML::Node& subnode: node) {
			emplace_back(PrimitiveBasis(subnode));
		}

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
