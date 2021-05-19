#include "KMC_Kmer.h"

KMC_Kmer::KMC_Kmer() {
	this->kmer_length = 0;
	this->signature_length = 0;
	this->prefix_length = 0;
	this->suffix_length = 0;
	this->both_strands = true;
}

KMC_Kmer::KMC_Kmer(const uint32 kmer_length, const uint32 signature_length,
		const uint32 prefix_length, const bool both_strands) {
	//checks
	uint32 suffix_length = kmer_length - prefix_length;
	if (suffix_length % 4 == 0) {
		this->kmer_length = kmer_length;
		this->signature_length = signature_length;
		this->prefix_length = prefix_length;
		this->suffix_length = suffix_length;
		this->both_strands = both_strands;
	} else {
		std::cerr << "Invalid suffix length" << std::endl;
		exit(1);
	}
}

uint32 KMC_Kmer::values2signature(const uchar* values) {
	uint32 value = 1 << 2 * signature_length;
	uint32 norm = 0, rnorm = 0, mask = (1 << 2 * signature_length) - 1;
	for (uint32 i = 0; i < kmer_length; i++) {
		norm = norm << 2;
		rnorm = rnorm << 2;
		norm += values[i];
		rnorm += 3 - values[kmer_length - 1 - i];
		norm &= mask;
		rnorm &= mask;
		if (i >= signature_length - 1) {
			if (norm < value && signature_allowed(norm, signature_length)) {
				value = norm;
			}
			if (rnorm < value && signature_allowed(rnorm, signature_length)) {
				value = rnorm;
			}
		}
	}
	return value;
}

Kmer KMC_Kmer::create(const uint32 prefix, const uchar* suffix) {
	return create(prefix2string(prefix, false) + suffix2string(suffix, false));
}

Kmer KMC_Kmer::create(const std::string symbols) {
	return create(symbols.c_str());
}

Kmer KMC_Kmer::create(const char* symbols) {
	Kmer kmer(kmer_length,prefix_length,signature_length,false);

	//normalize to canonical if necessary
	char nsymbols[kmer_length];
	strncpy(nsymbols, symbols, kmer_length);
	kmer.reversed = both_strands ? canonical_kmer(nsymbols) : false;

	//compute signature
	uchar values[kmer_length];
	symbols2values(nsymbols, kmer_length, values);
	kmer.signature = values2signature(values);

	//define prefix and suffix
	uint32 i, j;
	for (i = 0; i < prefix_length; i++) {
		kmer.prefix = kmer.prefix << 2;
		kmer.prefix += values[i];
	}
	for (i = 0; i < suffix_length; i += 4) {
		kmer.suffix[i / 4] = 0;
		for (j = 0; j < 4; j++) {
			kmer.suffix[i / 4] += values[prefix_length + i + j]
					<< (6 - (2 * j));
		}
	}
    return kmer;
}

Kmer KMC_Kmer::create(const Kmer kmer1) {
	Kmer kmer(kmer1.klength, kmer1.plength,kmer1.slength, kmer1.reversed);
    kmer.prefix = kmer1.prefix;
    kmer.signature = kmer1.signature;
    for(uchar i=0; i<(kmer1.klength - kmer1.plength) / 4; i++) {
		kmer.suffix[i]=kmer1.suffix[i];
	}
    return kmer;
}

bool KMC_Kmer::canonical_kmer(char* symbols) {
	uint32 i;
	char rsymbols[kmer_length];
	reverse_symbols(symbols, rsymbols);
	for (i = 0; i < kmer_length; i++) {
		if (symbols[i] < rsymbols[i]) {
			return false; //not reversed
		} else if (symbols[i] > rsymbols[i]) {
			strncpy(symbols, rsymbols, kmer_length);
			return true; //reversed
		}
	}
	return false; //not reversed
}

void KMC_Kmer::reverse_symbols(const char* symbols, char* rsymbols) {
	unsigned int i;
	for (i = 0; i < kmer_length; i++) {
		switch (symbols[i]) {
		case 'A':
			rsymbols[kmer_length - 1 - i] = 'T';
			break;
		case 'C':
			rsymbols[kmer_length - 1 - i] = 'G';
			break;
		case 'G':
			rsymbols[kmer_length - 1 - i] = 'C';
			break;
		case 'T':
			rsymbols[kmer_length - 1 - i] = 'A';
			break;
		default:
			rsymbols[kmer_length - 1 - i] = 'N';
			break;
		}
	}
}

void KMC_Kmer::expand_kmer(const Kmer kmer, const uint32 mm,
		std::set<Kmer> &expansion, std::map<Kmer, std::string> &types) {
	static std::string letters[4] = { "A", "C", "G", "T" };
	if (mm > 0) {
		std::string symbols = kmer2string(kmer);
		std::string nsymbols;
		for (uint32 i = 0; i < kmer.klength; i++) {
			for (std::string letter : letters) {
				nsymbols = std::string(symbols).replace(i, 1, letter);
				if (symbols[i] != letter[0]) {
					types[create(nsymbols)] = symbols[i] + letter;
				}
				expand_kmer(create(nsymbols), mm - 1, expansion, types);
			}
		}
	} else {
		expansion.insert(kmer);
	}
}

void KMC_Kmer::substitutions(const Kmer kmer, std::set<Kmer> &list,
		std::map<Kmer, std::string> &types, std::map<Kmer, uchar> &numbers) {
	static std::string letters[4] = { "A", "C", "G", "T" };
	std::string symbols = kmer2string(kmer);
	std::string nsymbols;
	for (uint32 i = 0; i < kmer.klength; i++) {
		for (std::string letter : letters) {
			nsymbols = std::string(symbols).replace(i, 1, letter);
			Kmer nkmer = create(nsymbols);
			types[nkmer] = symbols[i] + letter;
			numbers[nkmer] += 1;
			list.insert(nkmer);
		}
	}
}

void KMC_Kmer::deletions(const Kmer kmer, std::set<Kmer> &list,
		std::map<Kmer, std::string> &types, std::map<Kmer, uchar> &numbers) {
	static std::string letters[4] = { "A", "C", "G", "T" };
	std::string symbols = kmer2string(kmer);
	std::string nsymbols;
	//delete in the middle
	for (uint32 i = 1; i < kmer.klength - 1; i++) {
		for (std::string letter : letters) {
			nsymbols = std::string(symbols).erase(i, 1) + letter;
			Kmer nkmer = create(nsymbols);
			types[nkmer] = std::string(symbols).substr(i - 1, 3);
			numbers[nkmer] += 1;
			list.insert(nkmer);
		}
	}
	//delete last position, type depends on appended letter
	for (std::string letter : letters) {
		nsymbols = std::string(symbols).substr(0, kmer.klength - 1) + letter;
		Kmer nkmer = create(nsymbols);
		types[nkmer] = std::string(symbols).substr(kmer.klength - 2, 2)
				+ letter;
		numbers[nkmer] += 1;
		list.insert(nkmer);
	}
}

void KMC_Kmer::insertions(const Kmer kmer, std::set<Kmer> &list,
		std::map<Kmer, std::string> &types, std::map<Kmer, uchar> &numbers) {
	static std::string letters[4] = { "A", "C", "G", "T" };
	std::string symbols = kmer2string(kmer);
	std::string nsymbols;
	//add at begin
	for (std::string letter : letters) {
		nsymbols = letter + std::string(symbols).substr(0, kmer.klength - 1);
		Kmer nkmer = create(nsymbols);
		types[nkmer] = "*" + letter + std::string(symbols).substr(0, 1);
		numbers[nkmer] += 1;
		list.insert(nkmer);
	}
	//insert in the middle
	for (uint32 i = 1; i < kmer.klength; i++) {
		for (std::string letter : letters) {
			nsymbols = std::string(symbols).insert(i, letter).substr(0,
					kmer.klength);
			Kmer nkmer = create(nsymbols);
			types[nkmer] = std::string(nsymbols).substr(i - 1, 3);
			numbers[nkmer] += 1;
			list.insert(nkmer);
		}
	}
	//add at end
	for (std::string letter : letters) {
		nsymbols = std::string(symbols).substr(0, kmer.klength - 1) + letter;
		Kmer nkmer = create(nsymbols);
		types[nkmer] = std::string(symbols).substr(kmer.klength - 1, 1) + letter
				+ "*";
		numbers[nkmer] += 1;
		list.insert(nkmer);
	}
}

void KMC_Kmer::values2symbols(const uchar* values, const uint32 length,
		char* symbols) {
	unsigned int i;
	for (i = 0; i < length; i++) {
		switch (values[i]) {
		case 0:
			symbols[i] = 'A';
			break;
		case 1:
			symbols[i] = 'C';
			break;
		case 2:
			symbols[i] = 'G';
			break;
		case 3:
			symbols[i] = 'T';
			break;
		default:
			symbols[i] = 'N';
			break;
		}
	}
}

void KMC_Kmer::symbols2values(const char* symbols, const uint32 length,
		uchar* values) {
	unsigned int i;
	for (i = 0; i < length; i++) {
		switch (symbols[i]) {
		case 'A':
			values[i] = 0;
			break;
		case 'C':
			values[i] = 1;
			break;
		case 'G':
			values[i] = 2;
			break;
		case 'T':
			values[i] = 3;
			break;
		default:
			break;
		}
	}
}

bool KMC_Kmer::signature_allowed(const uint32 signature, const uint32 length) {
	uint32 value = signature;
	if ((value & 0x3f) == 0x3f) // TTT suffix
		return false;
	if ((value & 0x3f) == 0x3b) // TGT suffix
		return false;
	if ((value & 0x3c) == 0x3c) // TG* suffix
		return false;

	for (uint32 j = 0; j < length - 3; ++j)
		if ((value & 0xf) == 0) // AA inside
			return false;
		else
			value >>= 2;

	if (value == 0) // AAA prefix
		return false;
	if (value == 0x04) // ACA prefix
		return false;
	if ((value & 0xf) == 0) // *AA prefix
		return false;
	return true;
}

std::string KMC_Kmer::signature2string(const uint32 signature) {
	uchar* values = new uchar[signature_length];
	char* symbols = new char[signature_length];
	uint32 v = signature;
	unsigned int i;
	for (i = 0; i < signature_length; i++) {
		values[signature_length - 1 - i] = v & 3;
		v = v >> 2;
	}
	values2symbols(values, signature_length, symbols);
	std::string result = std::string(symbols, signature_length);
    delete[] values;
    delete[] symbols;
	return result;
}

std::string KMC_Kmer::prefix2string(const uint64 prefix, const bool reversed) {
	uchar* values = new uchar[prefix_length];
	char* symbols = new char[prefix_length];
	uint64 v = prefix;
	unsigned int i;
	if (reversed) {
		for (i = 0; i < prefix_length; i++) {
			values[i] = 3 - (v & 3);
			v = v >> 2;
		}
	} else {
		for (i = 0; i < prefix_length; i++) {
			values[prefix_length - i - 1] = v & 3;
			v = v >> 2;
		}
	}
	values2symbols(values, prefix_length, symbols);
	std::string result = std::string(symbols, prefix_length);
	delete[] values;
	delete[] symbols;
	return result;
}

std::string KMC_Kmer::suffix2string(const uchar* suffix, const bool reversed) {
	uint32 length = suffix_length / 4;
	uchar* values = new uchar[4 * length];
	char* symbols = new char[4 * length];
	uchar v;
	unsigned int i, j;
	if (reversed) {
		for (i = 0; i < length; i++) {
			v = suffix[length - i - 1];
			for (j = 0; j < 4; j++) {
				values[j + (4 * i)] = 3 - ((v >> (6 - (2 * (4 - j - 1)))) & 3);
			}
		}
	} else {
		for (i = 0; i < length; i++) {
			v = suffix[i];
			for (j = 0; j < 4; j++) {
				values[j + (4 * i)] = (v >> (6 - (2 * j))) & 3;
			}
		}
	}
	values2symbols(values, 4 * length, symbols);
	std::string result = std::string(symbols, 4 * length);
	delete[] values;
	delete[] symbols;
	return result;
}

std::string KMC_Kmer::kmer2string(const Kmer kmer) {
	std::stringstream ss;
	if (kmer.reversed) {
		ss << suffix2string(kmer.suffix, true);
		ss << prefix2string(kmer.prefix, true);
	} else {
		ss << prefix2string(kmer.prefix, false);
		ss << suffix2string(kmer.suffix, false);
	}
	return ss.str();
}
