//
// Created by Roman Ellerbrock on 8/18/21.
//

#include "yaml-cpp/yaml.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "Lanczos.h"
#include "Hamiltonian.h"
#include "Operators.h"
#include "Potentials.h"

using namespace YAML;

void execute_job(Wavefunction& Psi, const Hamiltonian& H,
	const Node& node, const Basis& basis) {

	auto job = read_key<string>(node, "job");

	if (job == "Eigenstates") {
		cout << "Running Eigenstate calculation...\n";
		auto nkrylov = read_key<size_t>(node, "krylov_size");

		Lanczos lan;
		lan.calculate(Psi, H, basis, nkrylov);
		lan.print();
	} else if (job == "Integrate") {
		auto nkrylov = read_key<size_t>(node, "krylov_size");
		auto dt = read_key<double>(node, "dt");
		auto tmax = read_key<double>(node, "tmax");
		auto out = read_key<double>(node, "out");
		auto accuracy = read_key<double>(node, "accuracy");
		Lanczos lan;
		double t = 0;
		lan.integrate(Psi, t, dt, tmax, out, H, basis, nkrylov);
	} else {
		cerr << "Cannot find job named: " << job << endl;
		exit(3);
	}
}

Hamiltonian read_hamiltonian(const Node& input, const Basis& basis) {
	/**
	 * Read Hamiltonian from yaml node.
	 */

	Node Hnode = read_key<Node>(input, "Hamiltonian");
	auto Hname = read_key<string>(Hnode, "name");

	Hamiltonian H;
	if (Hname == "HO") {
		H = harmonic_osciallator(basis);
	} else if (Hname == "KineticEnergy") {
		H = kinetic_energy(basis);
	} else {
		cerr << "wrong Hamiltonian name.\n";
		exit(3);
	}

	auto Vname = read_key<string>(Hnode, "V");
	if (Vname == "HO") {
		ProductOperator V;
		V.V_ = V_HO;
		H.emplace_back(1., V);
	} else if (Vname == "none") {
		/// do nothing
	} else {
		cerr << "Error: wrong PES name.\n";
		exit(3);
	}

	return H;
}

void run(const string& filename) {
	/// Load input
	Node input = LoadFile(filename);

	/// Read basis from corresponding key
	Node basis_node = read_key<Node>(input, "Basis");
	Basis basis(basis_node);

	/// Create Wavefunction
	Wavefunction Psi(basis);
	Psi.occupy(basis);

	/// Create Hamiltonian
	Hamiltonian H = read_hamiltonian(input, basis);

	/// Read what jobs to perform
	Node run_node = read_key<Node>(input, "run");
	execute_job(Psi, H, run_node, basis);
}