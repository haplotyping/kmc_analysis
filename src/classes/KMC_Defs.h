#ifndef CLASSES_KMC_DEFS_H_
#define CLASSES_KMC_DEFS_H_

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <cstring>

typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;
typedef unsigned char uchar;

struct Kmer {
	uint32 signature;
	uint64 prefix;
	uint32 klength; //length k-mer
	uint32 slength; //length signature
	uint32 plength; //length prefix
	uchar suffix[255];
	bool reversed;  //original reversed to make canonical

	Kmer(const Kmer &kmer1) : signature(kmer1.signature),
			prefix(kmer1.prefix), klength(kmer1.klength),
			slength(kmer1.slength), plength(kmer1.plength),
			reversed(kmer1.reversed)
	{
		if (klength != plength) {
			for(uchar i=0; i<(klength - plength) / 4; i++) {
				suffix[i]=kmer1.suffix[i];
			}
		}
	}

	Kmer(uint32 klength_ = 0, uint32 plength_ = 0, uint32 slength_ = 0,
			bool reversed_ = false) : signature(0),
					prefix(0), klength(klength_),
					slength(slength_), plength(plength_),
					reversed(false) {
	}

	~Kmer() {
	}
};

bool operator !=(const Kmer kmer1, const Kmer kmer2);
bool operator ==(const Kmer kmer1, const Kmer kmer2);
bool operator <(const Kmer kmer1, const Kmer kmer2);

void progress(float progress);
void progress(std::string text);

#endif /* CLASSES_KMC_DEFS_H_ */
