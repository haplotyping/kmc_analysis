#include "KMC_Db.h"

KMC_Db::KMC_Db(const std::string &file_name) :
		file_base(file_name) {
	if (!initialize()) {
		std::cerr << "couldn't initialize!" << std::endl;
	}

}

bool KMC_Db::initialize() {

	char marker[4];
	size_t result;

	//change only after succesfull initialisation
	initialised = false;

	//define filenames
	file_name_pre = file_base + file_name_pre_suffix;
	file_name_suf = file_base + file_name_suf_suffix;

	//open pre-file
	if ((file_pre = fopen(file_name_pre.c_str(), "rb")) == NULL) {
		std::cerr << "Couldn't open '" << file_name_pre << "'" << std::endl;
		return false;
	}

	//check start marker
	result = fread(marker, 1, 4, file_pre);
	if (result != 4 || (strncmp(marker_pre, marker, 4) != 0)) {
		fclose(file_pre);
		file_pre = NULL;
		return false;
	}

	//move to end, register size
	fseek(file_pre, 0, SEEK_END);
	size_pre = ftell(file_pre);

	//check end marker
	fseek(file_pre, -4, SEEK_CUR);
	result = fread(marker, 1, 4, file_pre);
	if (result != 4 || (strncmp(marker_pre, marker, 4) != 0)) {
		fclose(file_pre);
		file_pre = NULL;
		return false;
	}

	//check version, for now only 0x200 supported
	fseek(file_pre, -12, SEEK_END);
	result = fread(&kmc_version, sizeof(uint32), 1, file_pre);
	int64 header_offset;
	if (kmc_version == 0x200) {
		//get header offset
		fseek(file_pre, -8, SEEK_END);
		header_offset = fgetc(file_pre);

		//process header
		fseek(file_pre, (0LL - (header_offset + 8)), SEEK_END);
		result = fread(&kmer_length, 1, sizeof(uint32), file_pre);
		result = fread(&mode, 1, sizeof(uint32), file_pre);
		result = fread(&suffix_counter_size, 1, sizeof(uint32), file_pre);
		result = fread(&prefix_length, 1, sizeof(uint32), file_pre);
		result = fread(&signature_length, 1, sizeof(uint32), file_pre);
		result = fread(&min_count, 1, sizeof(uint32), file_pre);
		result = fread(&max_count, 1, sizeof(uint32), file_pre);
		result = fread(&total_kmers, 1, sizeof(uint64), file_pre);
		result = fread(&both_strands, 1, 1, file_pre);
		both_strands = !both_strands;
	} else if (kmc_version == 0x0) {
		uint64 value;
		//get header offset
		fseek(file_pre, -8, SEEK_END);
		header_offset = fgetc(file_pre);
		//process header
		fseek(file_pre, (0LL - (header_offset + 8)), SEEK_END);
		result = fread(&kmer_length, 1, sizeof(uint32), file_pre);
		result = fread(&mode, 1, sizeof(uint32), file_pre);
		result = fread(&suffix_counter_size, 1, sizeof(uint32), file_pre);
		result = fread(&prefix_length, 1, sizeof(uint32), file_pre);
		result = fread(&min_count, 1, sizeof(uint32), file_pre);
		result = fread(&max_count, 1, sizeof(uint32), file_pre);
		result = fread(&total_kmers, 1, sizeof(uint64), file_pre);
		result = fread(&value, 1, sizeof(uint64), file_pre);
		both_strands = (value & 0x000000000000000F) == 1;
		both_strands = !both_strands;
		max_count += value & 0xFFFFFFFF00000000;
		signature_length = 0;

	} else {
		std::cerr << "Version " << kmc_version << " not supported!"
				<< std::endl;
		fclose(file_pre);
		file_pre = NULL;
		return false;
	}

	//open suf-file
	if ((file_suf = fopen(file_name_suf.c_str(), "rb")) == NULL) {
		std::cerr << "Couldn't open '" << file_name_suf << "'" << std::endl;
		return false;
	}

	//check start marker
	result = fread(marker, 1, 4, file_suf);
	if (result != 4 || (strncmp(marker_suf, marker, 4) != 0)) {
		fclose(file_suf);
		file_suf = NULL;
		return false;
	}

	//move to end, register size
	fseek(file_suf, 0, SEEK_END);
	size_suf = ftell(file_suf);

	//check end marker
	fseek(file_suf, -4, SEEK_CUR);
	result = fread(marker, 1, 4, file_suf);
	if (result != 4 || (strncmp(marker_suf, marker, 4) != 0)) {
		fclose(file_suf);
		file_suf = NULL;
		return false;
	}

	//other variables
	signature_map_size = (1 << (2 * signature_length)) + 1;
	signature_map_position = size_pre - (signature_map_size * sizeof(uint32))
			- (header_offset + 8);
	prefixes_size = (size_pre - 4 - (header_offset + 8)
			- (signature_map_size * sizeof(uint32))) / sizeof(uint64);
	prefixes_list_size = 1 << (2 * prefix_length);
	if (kmc_version == 0x00) {
		prefixes_list_size -= 1;
	}
	suffix_size = (kmer_length - prefix_length) / 4;
	suffix_record_size = suffix_size + suffix_counter_size;
	suffixes_size = size_suf - 8;

	//object to create k-mers
	kmc_kmer = KMC_Kmer(kmer_length, signature_length, prefix_length,
			both_strands);

	//final check consistency for 0x200
	if (kmc_version == 0x200) {
		uint64 check_signature_map_size;
		check_signature_map_size = (prefixes_size - 1) * (prefixes_size - 1)
				/ (prefixes_list_size * prefixes_list_size) + 1;
		if (check_signature_map_size != signature_map_size) {
			return false;
		} else {
			initialised = true;
			return true;
		}
	} else {
		initialised = true;
		return true;
	}
}

void KMC_Db::info() {
	if(!initialised) {return;}
	std::cout << "Files" << std::endl;
	std::cout << "- prefix:" << file_name_pre << " (" << size_pre << " bytes)"
			<< std::endl;
	std::cout << "- suffix:" << file_name_suf << " (" << size_suf << " bytes)"
			<< std::endl;
	std::cout << "General" << std::endl;
	std::cout << "- size k-mer              : " << kmer_length << std::endl;
	std::cout << "- number of k-mers        : " << total_kmers << std::endl;
	std::cout << "- canonical form          : " << (both_strands ? "yes" : "no")
			<< std::endl;
	std::cout << "- minimum value counter   : " << min_count << std::endl;
	std::cout << "- maximum value counter   : " << max_count << std::endl;
	std::cout << "KMC" << std::endl;
	std::cout << "- version                 : 0x" << std::hex << kmc_version
			<< std::dec << std::endl;
	if (kmc_version == 0x200) {
		std::cout << "- signature length        : " << signature_length
				<< " symbols" << std::endl;
		std::cout << "- signature map position  : " << signature_map_position
				<< " in prefix file " << std::endl;
		std::cout << "- signature map size      : " << signature_map_size
				<< " signatures, "
				<< ((signature_map_size * sizeof(uint32)) + 1) << " bytes"
				<< std::endl;
	}
	std::cout << "- prefix length           : " << prefix_length << " symbols "
			<< std::endl;
	std::cout << "- prefixes lists position : " << prefixes_position
			<< " in prefix file" << std::endl;
	std::cout << "- prefixes lists size     : "
			<< (prefixes_size * sizeof(uint64)) << " bytes, "
			<< ((prefixes_size - 1) / prefixes_list_size) << " lists"
			<< std::endl;
	std::cout << "- prefixes list size      : " << prefixes_list_size
			<< " prefixes, " << (prefixes_list_size * sizeof(uint64))
			<< " bytes" << std::endl;
	if (kmc_version == 0x200) {
		std::cout << "- prefixes lists          : "
				<< ((prefixes_size - 1) / prefixes_list_size) << std::endl;
	}
	std::cout << "- suffix size             : " << suffix_size << " bytes, "
			<< (4 * suffix_size) << " symbols" << std::endl;
	std::cout << "- suffix counter_size     : " << suffix_counter_size
			<< " bytes" << std::endl;
	std::cout << "- suffix record size      : " << suffix_record_size
			<< " bytes" << std::endl;
	std::cout << "- suffix data position    : " << suffixes_position
			<< " in suffix file" << std::endl;
	std::cout << "- suffix data size        : " << suffixes_size << " bytes, "
			<< (suffixes_size / suffix_record_size) << " records" << std::endl;
}

void KMC_Db::dump(std::string output_file_name, uint32 min, uint32 max, bool rc) {
	if(!initialised) {return;}
	std::ofstream output_file(output_file_name);
	//create dump
	uint64 prefix_array_size = (1 << 2 * prefix_length);
	uint32 narrays;
	if (kmc_version == 0x200) {
		narrays = ((signature_map_position - prefixes_position - sizeof(uint64))
				/ (prefix_array_size * sizeof(uint64)));
	} else {
		narrays = 1;
	}
	//loop over prefixes
	uchar suffix_value[suffix_size];
	uint64 pos_suf, pos_suf_next, prefix_value;
	uint32 number, step_size = suffix_size + suffix_counter_size;
	fseek(file_pre, prefixes_position, SEEK_SET);
	if (kmc_version == 0x200) {
		fread(&pos_suf, sizeof(uint64), 1, file_pre);
		fseek(file_suf, 0, SEEK_SET);
	} else {
		pos_suf = 0;
	}
	for (uint32 carrays = 0; carrays < narrays; carrays++) {
		progress(1.0 * ftell(file_suf) / size_suf);
		for (prefix_value = 0; prefix_value < prefix_array_size;
				prefix_value++) {
			fread(&pos_suf_next, sizeof(uint64), 1, file_pre);
			if (pos_suf < pos_suf_next) {
				fseek(file_suf, suffixes_position + (pos_suf * (step_size)),
				SEEK_SET);
				uint64 position = pos_suf;
				while (position < pos_suf_next) {
					number = 0;
					fread(&suffix_value, sizeof(uchar), suffix_size, file_suf);
					fread(&number, suffix_counter_size, 1, file_suf);
					if ((!min || number >= min) && (!max || number <= max)) {
						if (kmc_version == 0x200) {
							output_file
									<< kmc_kmer.prefix2string(prefix_value,
											false);
						} else {
							output_file
									<< kmc_kmer.prefix2string(prefix_value - 1,
											false);
						}
						output_file
								<< kmc_kmer.suffix2string(suffix_value, false);
						output_file << "\t" << number << std::endl;
						if (rc) {
							output_file
										<< kmc_kmer.suffix2string(suffix_value, true);
							if (kmc_version == 0x200) {
								output_file
										<< kmc_kmer.prefix2string(prefix_value,
												true);
							} else {
								output_file
										<< kmc_kmer.prefix2string(prefix_value - 1,
												true);
							}
							output_file << "\t" << number << std::endl;
						}
					}
					position++;
				}
			}
			pos_suf = pos_suf_next;
		}
	}
	progress(std::string("Writing to file " + output_file_name));
}

void KMC_Db::search(std::vector<Kmer> kmers) {
	if(!initialised) {return;}
	std::map<Kmer, uint32> result;
	search(kmers, result);
	//print them
	std::map<Kmer, uint32>::iterator it;
	for (Kmer kmer : kmers) {
		if (result[kmer] > 0) {
			std::cout << kmc_kmer.kmer2string(kmer);
			std::cout << "\t" << result[kmer] << std::endl;
		}
	}

}

void KMC_Db::search(std::set<Kmer> kmers, std::map<Kmer, uint32> &result) {
	if(!initialised) {return;}
	std::vector<Kmer> vkmers(kmers.size());
	std::copy(kmers.begin(), kmers.end(), vkmers.begin());
	search(vkmers, result);
}

void KMC_Db::search(std::vector<Kmer> kmers, std::map<Kmer, uint32> &result) {
	if(!initialised) {return;}
	if (kmers.size() > 0) {
		//find signatures
		std::set<uint32> signatures;
		std::vector<Kmer>::const_iterator it;
		for (it = kmers.begin(); it != kmers.end(); it++) {
			signatures.insert((*it).signature);
		}
		//loop over signatures
		uint32 map_prefix_position;
		std::vector<Kmer> signature_kmers;
		if (kmc_version == 0x200) {
			for (uint32 signature : signatures) {
				//get prefixes location from signature map
				fseek(file_pre,
						signature_map_position + (signature * sizeof(uint32)),
						SEEK_SET);
				fread(&map_prefix_position, sizeof(uint32), 1, file_pre);
				//find and kmers for this signature
				signature_kmers.clear();
				for (it = kmers.begin(); it != kmers.end(); it++) {
					if ((*it).signature == signature) {
						signature_kmers.push_back((*it));
					}
				}
				//search them
				search_prefixes(signature_kmers, map_prefix_position, result);
			}
		} else {
			search_prefixes(kmers, 0, result);
		}
	}
}

void KMC_Db::search(Kmer kmer, std::map<Kmer, uint32> &result) {
	if(!initialised) {return;}
	uint32 map_prefix_position;
	std::vector<Kmer> signature_kmers;
	//get prefixes location from signature map
	fseek(file_pre, signature_map_position + (kmer.signature * sizeof(uint32)),
	SEEK_SET);
	fread(&map_prefix_position, sizeof(uint32), 1, file_pre);
	//search
	signature_kmers.push_back(kmer);
	search_prefixes(signature_kmers, map_prefix_position, result);
}

void KMC_Db::search_prefixes(std::vector<Kmer> kmers,
		const uint32 map_prefix_position, std::map<Kmer, uint32> &result) {
	//get prefixes
	if(!initialised) {return;}
	std::set<uint32> prefixes;
	std::vector<Kmer>::iterator it;
	for (it = kmers.begin(); it != kmers.end(); it++) {
		prefixes.insert((*it).prefix);
	}
	//loop over prefixes
	uint64 map_suffixes_position_start, map_suffixes_position_end;
	std::vector<Kmer> prefix_kmers;
	for (uint32 prefix : prefixes) {
		//get suffixes location from prefixes
		uint64 pp = prefixes_position;
		pp += map_prefix_position * (1 << (2 * prefix_length)) * sizeof(uint64);
		pp += prefix * (sizeof(uint64));
		fseek(file_pre, pp, SEEK_SET);
		fread(&map_suffixes_position_start, sizeof(uint64), 1, file_pre);
		fread(&map_suffixes_position_end, sizeof(uint64), 1, file_pre);
		if (map_suffixes_position_start < map_suffixes_position_end) {
			prefix_kmers.clear();
			for (it = kmers.begin(); it != kmers.end(); it++) {
				if ((*it).prefix == prefix) {
					prefix_kmers.push_back((*it));
				}
			}
			//search them
			search_suffixes(prefix, prefix_kmers, map_suffixes_position_start,
					map_suffixes_position_end, result);
		}
	}
}

void KMC_Db::search_suffixes(const uint32 prefix, std::vector<Kmer> kmers,
		const uint64 map_suffixes_position_start,
		const uint64 map_suffixes_position_end,
		std::map<Kmer, uint32> &result) {
	if(!initialised) {return;}
	//sort by suffixes
	sort(kmers.begin(), kmers.end());
	//search them
	uchar suffix_value[suffix_size];
	bool previous_reversed;
	uint64 position = map_suffixes_position_start;
	uint32 number;
	int comparison;
	bool match;
	fseek(file_suf,
			suffixes_position
					+ (map_suffixes_position_start * suffix_record_size),
			SEEK_SET);
	std::vector<Kmer>::iterator sit = kmers.begin();
	//loop over sorted suffixes
	while (position < map_suffixes_position_end) {
		fread(&suffix_value, sizeof(uchar), suffix_size, file_suf);
		fread(&number, suffix_counter_size, 1, file_suf);
		comparison = compare_suffix((*sit).suffix, suffix_value);
		match = false;
		//compare and forward sorted kmers
		while (comparison <= 0) {
			//match
			if (comparison == 0) {
				//check if previous same reversed status
				if (match && (*sit).reversed == previous_reversed) {
					match = false;
				}
				//no previous match for same (normalized) k-mer
				if (!match) {
					match = true;
					if (suffix_counter_size < 4) {
						number = number & ((1 << suffix_counter_size * 8) - 1);
					}
					previous_reversed = (*sit).reversed;
					//print original k-mer from query, not normalized version
					result[kmc_kmer.create((*sit))] = number;
				}
			}
			//update position in sorted k-mers
			sit++;
			if (sit == kmers.end())
				return;
			comparison = compare_suffix((*sit).suffix, suffix_value);
		}
		//update position in suffixes
		position++;
	}
}

int KMC_Db::compare_suffix(const uchar* suffix1, const uchar* suffix2) {
	uint32 i;
	for (i = 0; i < suffix_size; i++) {
		if (suffix1[i] < suffix2[i]) {
			return -1;
		} else if (suffix1[i] > suffix2[i]) {
			return 1;
		}
	}
	return 0;
}
