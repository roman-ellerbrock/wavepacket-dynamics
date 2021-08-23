#include <iostream>
#include "yaml-parser.h"
#include <fstream>

using namespace std;

int main(int argc, char *argv[]) {
	cout << "========================================" << endl;
	cout << "====   Wavepacket Dynamics program  ====" << endl;
	cout << "========================================" << endl;

	if (argc != 2) {
		cerr << "Please provide a mctdh input as argument." << endl;
		exit(1);
	}

	string filename(argv[1]);
	ifstream ifs{filename};

	if (ifs.fail()) {
		cerr << "Error opening the input file. Please make sure the Spelling\n"
			 << "is correct." << endl;
		exit(1);
	}

	/// call run, found in src/yaml-parser.cpp
	run(filename);

	return 0;
}
