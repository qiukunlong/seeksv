/*
 ***************************************
 * getsv.h                             *
 *                                     *
 *  Created on:   2011-5-24            *
 *  Version:      0.1.0                *
 *  seeksv version:      0.6.0         *
 *  Author:       qiukl                *
 *                                     *
 ***************************************
 */
#pragma once
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include "getsv.h"
#include "clip_reads.h"


using namespace std;



/*@function                   
 @abstract              read  breakpoint information of tumor sample and output somatic breakpoints
 @param  file           filename of normal breakpoint file
 @param  output_file    output filename of tumor breakpoint file
*/
void ReadTumorFileAndOutputSomaticInfo(string normal_bam_file, string file, string output_file, multimap<pair<string, int>, ReadsInfo> &aligned_pos2reads_3clip, multimap<pair<string, int>, ReadsInfo> &aligned_pos2reads_5clip, int min_mapQ, int offset, double min_map_rate, int mean_insert_size, int deviation, int times);

/*@function                      
 @abstract                            Read soft-clipped reads file and store the reads.
 @param  file                         soft-clipped reads file.
 @param  aligned_pos2reads_3clip      the container used to store soft-clipped reads of which 3' end is clipped.
 @param  aligned_pos2reads_5clip      the container used to store soft-clipped reads of which 5' end is clipped.
 @param  min_len_of_clipped_seq       some clipped sequences are too short to confirm a breakpoint , this parameter is used to filter the soft-clipped re                                      reads of which clipped sequence is too short
 DATE: 2012-08-06
*/
template <typename T> void ReadsClipReads(string file, multimap<pair<string, int>, ReadsInfo> &aligned_pos2reads_3clip, multimap<pair<string, int>, ReadsInfo> &aligned_pos2reads_5clip,int min_len_of_clipped_seq)
{
	string chr, read_id, cigar, aligned_seq, aligned_qual, clipped_seq, clipped_qual, temp;
	int support_count, pos, mismatch_no;
	char orientation;
	T fin(file.c_str());
	if (!fin) { cerr << "Cannot open file " << file << endl; exit(1); }

	while (fin >> chr)
	{
		fin >> pos >> orientation >> cigar >> aligned_seq >> aligned_qual >> clipped_seq >> clipped_qual >> support_count;
		getline(fin, temp);
		vector<pair<int, char> > cigar_vec = ChangeCigarType(cigar);

		if (clipped_seq.length() < min_len_of_clipped_seq) continue;
		if (orientation == '3')
		{
			ReadsInfo reads_info(aligned_seq, aligned_qual, clipped_seq, clipped_qual, cigar_vec, support_count, 0);
			aligned_pos2reads_3clip.insert(make_pair(make_pair(chr, pos), reads_info));
		}	
		else if (orientation == '5')
		{
			ReadsInfo reads_info(clipped_seq, clipped_qual, aligned_seq, aligned_qual, cigar_vec, support_count, 0);
			aligned_pos2reads_5clip.insert(make_pair(make_pair(chr, pos), reads_info));
		}
		else
		{
			cerr << "Error:The orientation of soft-clipped reads must be 3 or 5 in position " << chr << ":" << pos << endl;
		}
	}
}
