#include "kmc_analysis.h"

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cout << "KMC analysis" << std::endl;
		std::cout << "------------" << std::endl;
		std::cout << "Usage:" << std::endl;
		std::cout << " kmc_analysis <operation> [operation parameters]"
				<< std::endl;
		std::cout << "Available operations:" << std::endl;
		std::cout << " info           - information about KMC database"
				<< std::endl;
		std::cout << " dump           - dump k-mers" << std::endl;
	} else {
		std::string operation = std::string(argv[1]);
		if (operation == "info") {
			operation_info(argc, argv);
		} else if (operation == "dump") {
			operation_dump(argc, argv);
		} else {
			std::cerr << "Unknown operation '" << operation
					<< "', program aborted" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	return 0;
}

void operation_info(int argc, char** argv) {
	std::cout << "KMC analysis - info" << std::endl;
	std::cout << "-------------------" << std::endl;
	if (argc < 3) {
		std::cout << "Usage:" << std::endl;
		std::cout << " kmc_analysis info <input>" << std::endl;
		std::cout << "With:" << std::endl;
		std::cout << " input: path to KMC database" << std::endl;
	} else {
		std::string input_file_name;
		input_file_name = std::string(argv[2]);
		KMC_Db kmc_db(input_file_name);
		kmc_db.info();
	}
}

void operation_dump(int argc, char** argv) {
	uint32 min = 2;
	uint32 max = 255;
	bool rc = false;
	std::cout << "KMC analysis - dump" << std::endl;
	std::cout << "-------------------" << std::endl;
	if (argc < 4) {
		std::cout << "Usage:" << std::endl;
		std::cout << " kmc_analysis dump <input> <output> [parameters]"
				<< std::endl;
		std::cout << "With:" << std::endl;
		std::cout << " input         - path to KMC database" << std::endl;
		std::cout << " output        - path to output file" << std::endl;
		std::cout << "And optional parameters:" << std::endl;
		std::cout
				<< " -min <value>: exclude k-mers occurring less than <value> times (default "
				<< min << ")" << std::endl;
		std::cout
				<< " -max <value>: exclude k-mers occurring more than <value> times (default "
				<< max << ")" << std::endl;
		std::cout
				<< " -rc: include reverse complements" << std::endl;
	} else {
		std::string input_file_name, output_file_name;
		input_file_name = std::string(argv[2]);
		output_file_name = std::string(argv[3]);
		KMC_Db kmc_db(input_file_name);
		//parse parameters
		std::string parameter;
		for (int i = 4; i < argc; i++) {
			parameter = std::string(argv[i]);
			if (parameter == "-rc") {
				rc = true;
			} else if (parameter == "-min" || parameter == "-max") {
				if ((i + 1) < argc) {
					i++;
					std::string value = std::string(argv[i]);
					if (!std::regex_match(value, std::regex("^[0-9]+$"))) {
						std::cerr << "Ignored invalid value: '" << value << "'"
								<< " for '" << parameter << "'" << std::endl;
					} else {
						if (parameter == "-min") {
							min = std::stoi(value);
						} else if (parameter == "-max") {
							max = std::stoi(value);
						}
					}
				} else {
					std::cerr << "No value provided after '" << parameter
							<< "' parameter, program aborted" << std::endl;
					exit(EXIT_FAILURE);
				}
			} else {
				std::cerr << "Invalid parameter '" << parameter
						<< "', program aborted" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		if (min && max && max < min) {
			std::cerr
					<< "Minimum cannot be larger than maximum, program aborted"
					<< std::endl;
			exit(EXIT_FAILURE);
		}
		kmc_db.dump(output_file_name, min, max, rc);
	}
}
