#include "KMC_Defs.h"

bool operator !=(const Kmer kmer1, const Kmer kmer2) {
	return (kmer1 == kmer2) ? false : true;
}

bool operator ==(const Kmer kmer1, const Kmer kmer2) {
	if (kmer1.klength != kmer2.klength) {
		return false;
	} else if (kmer1.slength != kmer2.slength) {
		return false;
	} else if (kmer1.plength != kmer2.plength) {
		return false;
	} else if (kmer1.prefix != kmer2.prefix) {
		return false;
	} else if (kmer1.signature != kmer2.signature) {
		return false;
	} else {
		uchar i;
		uchar c = (kmer1.klength - kmer1.plength) / 4;
		for (i = 0; i < c; i++) {
			if (kmer1.suffix[i] != kmer2.suffix[i]) {
				return false;
			}
		}
	}
	return true;
}

bool operator <(const Kmer kmer1, const Kmer kmer2) {
	if (kmer1.klength != kmer2.klength) {
		return kmer1.klength < kmer2.klength;
	} else if (kmer1.slength != kmer2.slength) {
		return kmer1.slength < kmer2.slength;
	} else if (kmer1.plength != kmer2.plength) {
		return kmer1.plength < kmer2.plength;
	} else if (kmer1.prefix != kmer2.prefix) {
		return kmer1.prefix < kmer2.prefix;
	} else if (kmer1.signature != kmer2.signature) {
		return kmer1.signature < kmer2.signature;
	} else {
		uchar i;
		uchar c = (kmer1.klength - kmer1.plength) / 4;
		for (i = 0; i < c; i++) {
			if (kmer1.suffix[i] != kmer2.suffix[i]) {
				return kmer1.suffix[i] < kmer2.suffix[i];
			}
		}
	}
	return false;
}
void progress(float progress) {
	int barWidth = 70;
	std::cout << "[";
	int pos = barWidth * progress;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos)
			std::cout << "=";
		else if (i == pos)
			std::cout << ">";
		else
			std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();
}

void progress(std::string text) {
	std::cout << text;
	if (text.size() < 80) {
		std::cout << std::string(80 - text.size(), ' ');
	}
	std::cout << std::endl;
}
