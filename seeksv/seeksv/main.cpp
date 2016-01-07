//#include "sam_info.h"
#include "clip_reads.h"
#include "getsv.h"
#include <fstream>
#include <iostream>
#include "somatic.h"
#include "bam2depth.h"
#include "gzstream.h"
#include "cluster.h"
#include "process_bwasw.h"



int main(int argc, char *argv[])
{
	ifstream fin(argv[1]);
	if (!fin)
	{
		cerr << "Cannot open file " << argv[1] << endl;
		exit(1);
	}
	string chr, cigar, seq_left, qual_left, seq_right, qual_right;
	int pos, read_number;
	char side;
	multimap<pair<string, int>, ReadsInfo> breakpoint2read_l, breakpoint2read_r;

	while (fin >> chr)
	{
		fin >> pos >> side >> cigar >> seq_left >> qual_left >> seq_right >> qual_right >> read_number;
		ReadsInfo reads_info(seq_left, qual_left, seq_right, qual_right, cigar, 1, 0);
		//	reads_info.add_materead_id(read_id);
		breakpoint2read.insert(make_pair(make_pair(chr, pos), reads_info));
	}

	
	return 0;
}
