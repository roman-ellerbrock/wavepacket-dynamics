//
// Created by Roman Ellerbrock on 8/18/21.
//

#ifndef YAML_PARSER_H
#define YAML_PARSER_H
#include <string>
#include "yaml-helper.h"

/**
 * \brief run parser which reads all required information about the system and runs desired jobs
 * @param filename yaml-input file
 */
void run(const string& filename);

#endif //YAML_PARSER_H
