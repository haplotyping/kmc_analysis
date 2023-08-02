#ifndef CLASSES_KMC_DB_H_
#define CLASSES_KMC_DB_H_

#include <fstream>
#include "KMC_Defs.h"
#include "KMC_Kmer.h"

class KMC_Db {

private:

	FILE *file_pre;
	FILE *file_suf;
	FILE *file_suf_chunk;

	const char* marker_pre = "KMCP";
	const char* marker_suf = "KMCS";
	const char* file_name_pre_suffix = ".kmc_pre";
	const char* file_name_suf_suffix = ".kmc_suf";

	bool initialize();
	void search_prefixes(std::vector<Kmer> kmers, const uint32 map_prefix_position,
			std::map<Kmer, uint32> &result);
	void search_suffixes(const uint32 prefix, std::vector<Kmer> kmers,
			const uint64 map_suffixes_position_start,
			const uint64 map_suffixes_position_end, std::map<Kmer, uint32> &result);
	int compare_suffix(const uchar* suffix1, const uchar* suffix2);

public:

	bool initialised;

	uint64 size_pre;
	uint64 size_suf;

	uint32 kmer_length;
	uint32 mode;
	uint32 suffix_counter_size;
	uint32 prefix_length;
	uint32 signature_length;
	uint32 min_count;
	uint32 max_count;
	uint64 total_kmers;
	bool both_strands;
	uint32 kmc_version;

	uint64 signature_map_size;
	uint64 signature_map_position;
	uint32 prefixes_list_size;
	uint64 prefixes_size;
	uint64 prefixes_position = 4;

	uint32 suffix_size;
	uint32 suffix_record_size;
	uint64 suffixes_position = 4;
	uint64 suffixes_size;

	std::string file_base;
	std::string file_name_pre;
	std::string file_name_suf;

	KMC_Kmer kmc_kmer;

	KMC_Db(const std::string &file_name);
	~KMC_Db();

	void info();
	void dump(std::string output_file_name, uint32 min, uint32 max, bool rc);
	void search(std::vector<Kmer> kmers);
	void search(Kmer kmer, std::map<Kmer, uint32> &result);
	void search(std::set<Kmer> kmers, std::map<Kmer, uint32> &result);
	void search(std::vector<Kmer> kmers, std::map<Kmer, uint32> &result);
};

#endif /* CLASSES_KMC_DB_H_ */
