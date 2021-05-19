#include "kmc_query.h"

int main(int argc, char** argv) {
	uint32 mm = 0;
	if (argc < 3) {
		std::cout << "KMC query" << std::endl;
		std::cout << "---------" << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << " " << std::string(argv[0]) << " <input> [parameters]" << std::endl;
		std::cout << "With:" << std::endl;
		std::cout << " input: path to KMC database" << std::endl;
		std::cout << "And parameters:" << std::endl;
		std::cout << " -q, --query <kmer>: single k-mer query" << std::endl;
		std::cout << " -f, --queryfile <file>: multiple k-mer query from file"
				<< std::endl;
		std::cout << " -mm, --mismatch <number>: number of mismatches (default "
				<< mm << ")" << std::endl;
	} else {
		//open database
		std::string input_file_name;
		input_file_name = std::string(argv[1]);
		KMC_Db kmc_db(input_file_name);
		//define k-mers for query
		std::vector<Kmer> kmers;
		std::string parameter;
		for (int i = 2; i < argc; i++) {
			parameter = std::string(argv[i]);
			if (parameter == "-q" || parameter == "--query") {
				if ((i + 1) < argc) {
					i++;
					parameter = std::string(argv[i]);
					if (!std::regex_match(parameter, std::regex("^[ACGT]+$"))) {
						std::cerr << "Ignored invalid k-mer: '" << parameter
								<< "'" << std::endl;
					} else if (parameter.size() != kmc_db.kmer_length) {
						std::cerr << "Ignored k-mer because of invalid size: '"
								<< parameter << "'" << std::endl;
					} else {
						kmers.push_back(kmc_db.kmc_kmer.create(parameter));
					}
				} else {
					std::cerr
							<< "No k-mer provided in single k-mer query, program aborted"
							<< std::endl;
					exit(EXIT_FAILURE);
				}
			} else if (parameter == "-f" || parameter == "--queryfile") {
				if ((i + 1) < argc) {
					i++;
					parameter = std::string(argv[i]);
					std::ifstream in(parameter);
					std::string line;
					uint32 total = 0, ignored = 0;
					while (std::getline(in, line)) {
						line = std::regex_replace(line, std::regex("([ACGT]*)"),
								std::string("$1"));
						if (line.size() == kmc_db.kmer_length) {
							kmers.push_back(kmc_db.kmc_kmer.create(line));
							total++;
						} else {
							ignored++;
						}
					}
					if (total && ignored) {
						std::cerr << "Ignored " << ignored << " line(s) in "
								<< parameter << std::endl;
					} else if (!total) {
						std::cerr << "No valid k-mers in " << parameter
								<< std::endl;
					}
				} else {
					std::cerr
							<< "No filename provided in multiple k-mer query, program aborted"
							<< std::endl;
					exit(EXIT_FAILURE);
				}
			} else if (parameter == "-mm" || parameter == "--mismatch") {
				if ((i + 1) < argc) {
					i++;
					parameter = std::string(argv[i]);
					if (!std::regex_match(parameter, std::regex("^[0-9]+$"))) {
						std::cerr << "Ignored invalid minimum: '" << parameter
								<< "'" << std::endl;
					} else {
						mm = std::stoi(parameter);
					}
				} else {
					std::cerr << "No number provided after '" << parameter
							<< "' parameter, program aborted" << std::endl;
					exit(EXIT_FAILURE);
				}
			} else {
				std::cerr << "Invalid parameter '" << parameter
						<< "', program aborted" << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		if (kmers.size() > 0) {
			if (mm > 0) {
				std::set<Kmer> list;
				std::map<Kmer, std::string> types;
				for (Kmer kmer : kmers) {
					kmc_db.kmc_kmer.expand_kmer(kmer, mm, list, types);
				}
				kmers.clear();
				for (Kmer item : list) {
					kmers.push_back(item);
				}
			}
			//search them
			kmc_db.search(kmers);
		} else
			std::cerr << "No k-mers to search, program aborted" << std::endl;
		exit(EXIT_FAILURE);
	}
}
