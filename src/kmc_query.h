#ifndef KMC_QUERY_H_
#define KMC_QUERY_H_

#include "classes/KMC_Defs.h"
#include "classes/KMC_Db.h"
#include "classes/KMC_Kmer.h"

void operation_info(int argc, char** argv);
void operation_dump(int argc, char** argv);
void operation_query(int argc, char** argv);
void operation_histogram(int argc, char** argv);
void operation_error(int argc, char** argv);

#endif /* KMC_QUERY_H_ */
