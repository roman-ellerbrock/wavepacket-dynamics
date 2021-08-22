//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef YAML_HELPER_H
#define YAML_HELPER_H
#include "yaml-cpp/yaml.h"
#include <string>
#include <iostream>

using namespace std;

/// read key in YAML Node with safety check
template<class T>
T read_key(const YAML::Node& node, const string& key);

#endif //YAML_HELPER_H
