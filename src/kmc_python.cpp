#include "kmc_python.h"

extern "C" {

	bool status(const char *location, uint64_t *status) {
		KMC_Db *kmc_db = new KMC_Db(location);
		if (kmc_db->initialised) {
			status[0] = kmc_db->size_pre;
			status[1] = kmc_db->size_suf;
			status[2] = kmc_db->kmer_length;
			status[3] = kmc_db->mode;
			status[4] = kmc_db->suffix_counter_size;
			status[5] = kmc_db->prefix_length;
			status[6] = kmc_db->signature_length;
			status[7] = kmc_db->min_count;
			status[8] = kmc_db->max_count;
			status[9] = kmc_db->total_kmers;
			status[10] = kmc_db->both_strands;
			status[11] = kmc_db->kmc_version;
			status[12] = kmc_db->signature_map_size;
			status[13] = kmc_db->signature_map_position;
			status[14] = kmc_db->prefixes_list_size;
			status[15] = kmc_db->prefixes_size;
			status[16] = kmc_db->prefixes_position;
			status[17] = kmc_db->suffix_size;
			status[18] = kmc_db->suffix_record_size;
			status[19] = kmc_db->suffixes_position;
			status[20] = kmc_db->suffixes_size;
			delete kmc_db;
			return true;
		}
		delete kmc_db;
		return false;
	}

	bool kmer_frequencies(const char *location, const char **kmer_list,
			const size_t n, uint32_t *frequencies,
			uint32_t * stats, uint64_t *status) {
		std::vector<Kmer> kmers;
		std::map<Kmer, uint32> result;
		std::string kmer;
		KMC_Db *kmc_db = new KMC_Db(location);
		if (kmc_db->initialised) {
			status[0] = kmc_db->size_pre;
			status[1] = kmc_db->size_suf;
			status[2] = kmc_db->kmer_length;
			status[3] = kmc_db->mode;
			status[4] = kmc_db->suffix_counter_size;
			status[5] = kmc_db->prefix_length;
			status[6] = kmc_db->signature_length;
			status[7] = kmc_db->min_count;
			status[8] = kmc_db->max_count;
			status[9] = kmc_db->total_kmers;
			status[10] = kmc_db->both_strands;
			status[11] = kmc_db->kmc_version;
			status[12] = kmc_db->signature_map_size;
			status[13] = kmc_db->signature_map_position;
			status[14] = kmc_db->prefixes_list_size;
			status[15] = kmc_db->prefixes_size;
			status[16] = kmc_db->prefixes_position;
			status[17] = kmc_db->suffix_size;
			status[18] = kmc_db->suffix_record_size;
			status[19] = kmc_db->suffixes_position;
			status[20] = kmc_db->suffixes_size;
			stats[0]=0;//checked
			stats[1]=0;//positive
			stats[2]=0;//minimum
			stats[3]=0;//maximum
			for (uint i = 0; i < n; i++) {
				kmer = kmer_list[i];
				if (!std::regex_match(kmer, std::regex("^[ACGT]+$"))) {
					return false;
				} else if (kmer.size() != kmc_db->kmer_length) {
					return false;
				} else {
					kmers.push_back(kmc_db->kmc_kmer.create(kmer));
					stats[0]++;
				}
			}
			kmc_db->search(kmers, result);
			uint32_t f;
			stats[0]=kmers.size();
			for (uint i = 0; i < kmers.size(); i++) {
				f = result[kmers[i]];
				frequencies[i] = f;
				if (f>0) {
					if(stats[1]==0) {
						stats[1] = 1;
						stats[2] = f;
						stats[3] = f;
					} else {
						stats[1]++;
						stats[2] = std::min(stats[2],f);
						stats[3] = std::max(stats[3],f);
					}
				}
			}
			delete kmc_db;
			return true;
		}
		delete kmc_db;
		return false;
	}

	bool kmer_frequencies_mm(const char *location, const char **kmer_list,
			const size_t n, const size_t mm, const char *output,
			uint32_t * stats, uint64_t *status) {
		std::vector<Kmer> kmers;
		std::map<Kmer, uint32> result;
		std::string kmer;
		KMC_Db *kmc_db = new KMC_Db(location);
		if (kmc_db->initialised) {
			status[0] = kmc_db->size_pre;
			status[1] = kmc_db->size_suf;
			status[2] = kmc_db->kmer_length;
			status[3] = kmc_db->mode;
			status[4] = kmc_db->suffix_counter_size;
			status[5] = kmc_db->prefix_length;
			status[6] = kmc_db->signature_length;
			status[7] = kmc_db->min_count;
			status[8] = kmc_db->max_count;
			status[9] = kmc_db->total_kmers;
			status[10] = kmc_db->both_strands;
			status[11] = kmc_db->kmc_version;
			status[12] = kmc_db->signature_map_size;
			status[13] = kmc_db->signature_map_position;
			status[14] = kmc_db->prefixes_list_size;
			status[15] = kmc_db->prefixes_size;
			status[16] = kmc_db->prefixes_position;
			status[17] = kmc_db->suffix_size;
			status[18] = kmc_db->suffix_record_size;
			status[19] = kmc_db->suffixes_position;
			status[20] = kmc_db->suffixes_size;
			for (uint i = 0; i < n; i++) {
				kmer = kmer_list[i];
				if (!std::regex_match(kmer, std::regex("^[ACGT]+$"))) {
					return false;
				} else if (kmer.size() != kmc_db->kmer_length) {
					return false;
				} else {
					kmers.push_back(kmc_db->kmc_kmer.create(kmer));
				}
			}
			if (mm > 0) {
				std::set<Kmer> list;
				std::map<Kmer, std::string> types;
				for (Kmer kmer : kmers) {
					kmc_db -> kmc_kmer.expand_kmer(kmer, mm, list, types);
				}
				kmers.clear();
				for (Kmer item : list) {
					kmers.push_back(item);
				}
			}
			kmc_db->search(kmers, result);
			uint32_t f;
			stats[0]=kmers.size();
			std::ofstream myfile;
			myfile.open (output);
			for (uint i = 0; i < kmers.size(); i++) {
				f = result[kmers[i]];
				if (f>0) {
					myfile << kmc_db -> kmc_kmer.kmer2string(kmers[i]);
					myfile << "\t" << f;
					myfile << std::endl;
					if(stats[1]==0) {
						stats[1] = 1;
						stats[2] = f;
						stats[3] = f;
					} else {
						stats[1]++;
						stats[2] = std::min(stats[2],f);
						stats[3] = std::max(stats[3],f);
					}
				}
			}
			myfile.close();
			delete kmc_db;
			return true;
		}
		delete kmc_db;
		return false;
	}
}

