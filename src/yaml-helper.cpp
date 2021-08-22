//
// Created by Roman Ellerbrock on 8/18/21.
//

#include "yaml-helper.h"
#include <iostream>

/// safe read for YAML Node
template<class T>
T read_key(const YAML::Node& node, const string& key) {
	if (auto par = node[key]) {
		return node[key].as<T>();
	} else {
		cerr << "Cannot find key in yaml-node. Key was: " << key << endl;
		cerr << "Yaml node: " << node << endl;
		exit(3);
	}
}

template double read_key<double>(const YAML::Node& node, const string& key);
template int read_key<int>(const YAML::Node& node, const string& key);
template size_t read_key<size_t>(const YAML::Node& node, const string& key);
template string read_key<string>(const YAML::Node& node, const string& key);
template YAML::Node read_key<YAML::Node>(const YAML::Node& node, const string& key);

