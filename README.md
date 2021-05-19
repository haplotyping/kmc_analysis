# KMC analysis

Additional tools to query a [KMC](https://github.com/refresh-bio/KMC) database.

## Installation

```
make
```

## Query

Query for a specific k-mer with `kmc_query`

```
KMC query
---------
Usage:
 bin/kmc_query <input> [parameters]
With:
 input: path to KMC database
And parameters:
 -q, --query <kmer>: single k-mer query
 -f, --queryfile <file>: multiple k-mer query from file
 -mm, --mismatch <number>: number of mismatches (default 0)
```

Query for a single k-mer

```
matthijs$ ./bin/kmc_query database.kmc -q "AACAATTTGAAATGGGGATCAGCTTGTGTAT"
AACAATTTGAAATGGGGATCAGCTTGTGTAT 20
```

Query for multiple k-mers

```
matthijs$ ./bin/kmc_query database.kmc -q "AAAAACAATTTGAAATGGGGATCAGCTTGTG" -q "CAATTTGAAATGGGGATCAGCTTGTGTATAT"
AAAAACAATTTGAAATGGGGATCAGCTTGTG 20
CAATTTGAAATGGGGATCAGCTTGTGTATAT 21
```

Query for k-mers from a file

```
matthijs$ ./bin/kmc_query database.kmc -f "kmer_list.txt"
ATTTGAAATGGGGATCAGCTTGTGTATATCT 21
TTTGAAATGGGGATCAGCTTGTGTATATCTC 21
TTGAAATGGGGATCAGCTTGTGTATATCTCG 21
```

Query for k-mers with a maximum number of mismatches

```
matthijs$ ./bin/kmc_query database.kmc -q "AACAATTTGAAATGGGGATCAGCTTGTGTAT" -mm 2
AACAATGTGAAATGGGGATCAGCTTGTGTAT 1
AACAATTTGAACTGGGGATCATCTTGTGTAT 1
AACAATTTGAAATGGGGATCAGCTTGTGTAT 20
```

## Analysis

Analysis with `kmc_analysis`

```
KMC analysis
------------
Usage:
 kmc_analysis <operation> [operation parameters]
Available operations:
 info           - information about KMC database
 dump           - dump k-mers
```

Information about KMC database

```
matthijs$ ./bin/kmc_analysis info database.kmc
KMC analysis - info
-------------------
Files
- prefix:database.kmc.kmc_pre (1310812 bytes)
- suffix:database.kmc.kmc_suf (4846685 bytes)
General
- size k-mer              : 31
- number of k-mers        : 440607
- canonical form          : yes
- minimum value counter   : 1
- maximum value counter   : 100000000
KMC
- version                 : 0x200
- signature length        : 9 symbols
- signature map position  : 262156 in prefix file 
- signature map size      : 262145 signatures, 1048581 bytes
- prefix length           : 3 symbols 
- prefixes lists position : 4 in prefix file
- prefixes lists size     : 262152 bytes, 512 lists
- prefixes list size      : 64 prefixes, 512 bytes
- prefixes lists          : 512
- suffix size             : 7 bytes, 28 symbols
- suffix counter_size     : 4 bytes
- suffix record size      : 11 bytes
- suffix data position    : 4 in suffix file
- suffix data size        : 4846677 bytes, 440607 records
```

Create list of all k-mers

```
KMC analysis - dump
-------------------
Usage:
 kmc_analysis dump <input> <output> [parameters]
With:
 input         - path to KMC database
 output        - path to output file
And optional parameters:
 -min <value>: exclude k-mers occurring less than <value> times (default 2)
 -max <value>: exclude k-mers occurring more than <value> times (default 255)
 -rc: include reverse complements
```
