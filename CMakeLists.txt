cmake_minimum_required(VERSION 3.16)
project(WavepacketDynamics C CXX)
enable_language(C)
enable_language(CXX)

cmake_policy(SET "CMP0076" NEW)
cmake_policy(SET "CMP0040" NEW)
cmake_policy(SET "CMP0042" NEW)

set(CMAKE_CXX_STANDARD 14)

### Add submodule for QuTree
# If you have QuTree installed, feel free to use this instead:
#find_package(QuTree)
add_subdirectory(external/QuTree)
include_directories("external/QuTree")

# If you have Yaml installed, feel free to use this instead:
#find_package(yaml-cpp REQUIRED)
#include_directories(${YAML_CPP_INCLUDE_DIR})
#add_library(Yaml INTERFACE IMPORTED)
add_subdirectory(external/yaml-cpp)
include_directories("external")

add_executable(WavepacketDynamics main.cpp)

target_sources(WavepacketDynamics PUBLIC
        src/Basis.h
        src/PrimitiveBasis.h

        src/ProductOperator.h
        src/Hamiltonian.h
        src/Wavefunction.h
        src/Operators.cpp
        src/Operators.h
        src/Potentials.cpp
        src/Potentials.h
        src/PrimitiveOperators.cpp
        src/PrimitiveOperators.h

        src/Lanczos.h

        src/yaml-parser.cpp
        src/yaml-parser.h
        src/yaml-helper.h
        src/yaml-helper.cpp
        )

target_include_directories(WavepacketDynamics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(WavepacketDynamics PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

# link binary with QuTree and Yaml
target_link_libraries(WavepacketDynamics QuTree::QuTree yaml-cpp)
