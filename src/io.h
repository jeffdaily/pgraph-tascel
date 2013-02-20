#ifndef FASTA_IO_H_
#define FASTA_IO_H_

#include <string>
#include <vector>

using std::string;
using std::vector;

/**
 * Reads the given file such that each process gets nl/p sequences.
 *
 * Uses MPIIO.
 */
void read_fasta(const char *filename, size_t budget, vector<string> &sequences);

#endif /* FASTA_IO_H_ */
