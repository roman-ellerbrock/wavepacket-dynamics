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

/**
 * \brief Execute a job (integrate in time or calculate eigenstates)
 * @param Psi Initial wavefunction
 * @param H Hamiltonian of the system
 * @param node yaml node that describes the job
 * @param basis Basis of the wavefunction
 */
void execute_job(Wavefunction& Psi, const Hamiltonian& H,
	const Node& node, const Basis& basis) {

	/// Read the type of job
	auto job = read_key<string>(node, "job");

	if (job == "Eigenstates") {
		/// Calculate Eigenstates by running a single Lanczos
		cout << "Running Eigenstate calculation...\n";
		auto nkrylov = read_key<size_t>(node, "krylov_size");
		Lanczos lan;
		lan.calculate(Psi, H, basis, nkrylov);
		lan.print();

	} else if (job == "Integrate") {
		/// Integrate in time using a short-iterative lanczos
		auto nkrylov = read_key<size_t>(node, "krylov_size");
		auto dt = read_key<double>(node, "dt");
		auto tmax = read_key<double>(node, "tmax");
		auto out = read_key<double>(node, "out");
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
	 * Rationale:
	 * - Read Hamiltonian from yaml node.
	 * - The Hamiltonians are implemented in Operators.cpp with PESs from Potentials.cpp
	 * - If you want to add your own Hamiltonian, just implement it in Operators.cpp/Potentials.cpp
	 *   and add them to the if/else below.
	 * - Feel free to play around by writing new Hamitonians
	 */

	Node Hnode = read_key<Node>(input, "Hamiltonian");
	auto Hname = read_key<string>(Hnode, "name");

	Hamiltonian H;
	if (Hname == "HO") {
		H = harmonic_osciallator(basis);
	} else if (Hname == "KineticEnergy") {
		H = kinetic_energy(basis);
	} else if (Hname == "HO_pes") {
		H = harmonic_osciallator_pes(basis);
	} else if (Hname == "Eckhard") {
		H = eckhard_pes(basis);
	} else if (Hname == "tullyA") {
		H = tully_A(basis);
	} else {
		cerr << "wrong Hamiltonian name.\n";
		exit(3);
	}

	return H;
}

void run(const string& filename) {
	/// Load input
	Node input = LoadFile(filename);

	/// Read Basis from corresponding key
	Node basis_node = read_key<Node>(input, "Basis");
	Basis basis(basis_node);


	/// For this exercise we will not need more than the primitive Basis
	exit(0);
	/// Create Wavefunction and occupy
	Wavefunction Psi(basis);
	Psi.occupy(basis);

	/// Create Hamiltonian
	Hamiltonian H = read_hamiltonian(input, basis);

	/// Read what jobs to perform
	Node run_node = read_key<Node>(input, "run");
	execute_job(Psi, H, run_node, basis);
}