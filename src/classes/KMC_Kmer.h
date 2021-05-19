#ifndef CLASSES_KMC_KMER_H_
#define CLASSES_KMC_KMER_H_

#include <cstring>
#include <regex>
#include "KMC_Defs.h"

class KMC_Kmer {

private:

	uint32 kmer_length;
	uint32 signature_length;
	uint32 prefix_length;
	uint32 suffix_length;
	bool both_strands;

	uint32 values2signature(const uchar* values);
	bool canonical_kmer(char* symbols);
	void reverse_symbols(const char* symbols, char* rsymbols);

	bool signature_allowed(const uint32 signature, const uint32 length);
	void symbols2values(const char* symbols, const uint32 length,
			uchar* values);
	void values2symbols(const uchar* values, const uint32 length,
			char* symbols);

public:

	KMC_Kmer();
	KMC_Kmer(const uint32 kmer_length, const uint32 signature_length,
			const uint32 prefix_length, const bool both_strands);

	Kmer create(const uint32 prefix, const uchar* suffix);
	Kmer create(const char* symbols);
	Kmer create(const std::string symbols);
	Kmer create(const Kmer kmer1);

	std::string kmer2string(const Kmer kmer);
	void expand_kmer(const Kmer kmer, const uint32 mm,
			std::set<Kmer> &expansion, std::map<Kmer, std::string> &types);
	void substitutions(const Kmer kmer, std::set<Kmer> &list,
			std::map<Kmer, std::string> &types, std::map<Kmer, uchar> &numbers);
	void deletions(const Kmer kmer, std::set<Kmer> &list,
			std::map<Kmer, std::string> &types, std::map<Kmer, uchar> &numbers);
	void insertions(const Kmer kmer, std::set<Kmer> &list,
			std::map<Kmer, std::string> &types, std::map<Kmer, uchar> &numbers);

	std::string signature2string(const uint32 value);
	std::string prefix2string(const uint64 value, const bool reversed);
	std::string suffix2string(const uchar* value, const bool reversed);

};

#endif /* CLASSES_KMC_KMER_H_ */
