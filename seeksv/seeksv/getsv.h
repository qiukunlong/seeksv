/*
 *********************************************
 * getsv.h                                   *
 *                                           *
 *  Created on:   2011-5-24    v0.1.0        * 
 *  Modified on:  2013-1-10    v0.2.0        * 
 *  seeksv version:      0.8.0               * 
 *  Author:       qiukl                      * 
 *                                           * 
 *********************************************
 */

#pragma once
#include <ctype.h>
#include <map>
#include <fstream>
#include <vector>
#include "clip_reads.h"
#include "gzstream.h"
#include "cluster.h"
#include "getsv.h"


class AlignInfo
{
public:
	string chr;
	int pos;
	int len;
	char strand;
	vector<pair<int, char> > cigar_vec;

	string seq;
	int left_clipped_seq_length;
	int right_clipped_seq_length;
	//uniq = 'u', repeat = 'r', none = 'n'
	char type;


	AlignInfo() { }
	AlignInfo(const string c, const int p, const int l, const char s, const vector<pair<int, char> > &ci, const string &seq1, const int& lcsl, const int& rcsl, const char ty) { chr = c; pos = p; len = l; strand = s; cigar_vec = ci; seq = seq1; left_clipped_seq_length = lcsl; right_clipped_seq_length = rcsl; type = ty; }
	void set_value(const string c, const int p, const int l, const char s, const vector<pair<int, char> > &ci, const string &seq1, const int& lcsl, const int& rcsl, const char ty) { chr = c; pos = p; len = l; strand = s; cigar_vec = ci; seq = seq1; left_clipped_seq_length = lcsl; right_clipped_seq_length = rcsl; type = ty; }

//	~AlignInfo();
};


class SeqInfo
{
public:
	string seq;
	vector<pair<int, char> > cigar_vec;
	int left_clipped_seq_length;
	int right_clipped_seq_length;
	int support_read_no;
//  Not clipped seq                       is_clipped_seq_and_uniq_mapped = 0
//  Clipped seq and Multimap              is_clipped_seq_and_uniq_mapped = 1
//  Clipped seq and uniq mapped           is_clipped_seq_and_uniq_mapped = 2
	int is_clipped_seq_and_uniq_mapped;

	SeqInfo() { }
	SeqInfo(const string s, const vector<pair<int, char> > &ci, const int lcsl, const int rcsl, const int read_no, const int clip_uniq) { seq = s; cigar_vec = ci; left_clipped_seq_length = lcsl; right_clipped_seq_length = rcsl; support_read_no = read_no; is_clipped_seq_and_uniq_mapped = clip_uniq; }
	SeqInfo(const SeqInfo &seq_info) { seq = seq_info.seq; cigar_vec = seq_info.cigar_vec; left_clipped_seq_length = seq_info.left_clipped_seq_length; right_clipped_seq_length = seq_info.right_clipped_seq_length;  support_read_no = seq_info.support_read_no; is_clipped_seq_and_uniq_mapped = seq_info.is_clipped_seq_and_uniq_mapped; }
	SeqInfo& operator= (const SeqInfo &s) { seq = s.seq; cigar_vec = s.cigar_vec; left_clipped_seq_length = s.left_clipped_seq_length; right_clipped_seq_length = s.right_clipped_seq_length; support_read_no = s.support_read_no; is_clipped_seq_and_uniq_mapped = s.is_clipped_seq_and_uniq_mapped; return *this; }
	void set_value(const string s, const vector<pair<int, char> > &ci, const int lcsl, const int rcsl, const int read_no, const int clip_uniq) { seq = s; cigar_vec = ci; left_clipped_seq_length = lcsl; right_clipped_seq_length = rcsl; support_read_no = read_no; is_clipped_seq_and_uniq_mapped = clip_uniq; }

	
};

class Breakpoint
{
public:
	string down_chr;
	int down_pos;
	int microhomology_length;
	//added in version v0.2.0
	int abnormal_read_pair_no;
	SeqInfo up_seq_info;
	SeqInfo down_seq_info;
	
	Breakpoint(const string dc, const int dp, const int trl, const int arpn, const SeqInfo &up_seq, const SeqInfo &down_seq)
	{
		down_chr = dc; down_pos = dp; microhomology_length = trl; abnormal_read_pair_no = arpn; up_seq_info = up_seq; down_seq_info = down_seq;
	}
	Breakpoint& operator= (const Breakpoint &b) { down_chr = b.down_chr; down_pos = b.down_pos; microhomology_length = b.microhomology_length; abnormal_read_pair_no = b.abnormal_read_pair_no; up_seq_info = b.up_seq_info; down_seq_info = b.down_seq_info; return *this; }
	
};

class OtherInfo
{
public:
	SeqInfo up_seq_info;
	SeqInfo down_seq_info;
	int microhomology_length;
	int abnormal_read_pair_no;
	
	OtherInfo() { }	
	OtherInfo(const SeqInfo &up_seq, const SeqInfo &down_seq, const int micro_ho_len, const int arpn)
	{
		up_seq_info = up_seq; down_seq_info = down_seq; microhomology_length = micro_ho_len; abnormal_read_pair_no = arpn;
	}
	void set_value(const SeqInfo &up_seq, const SeqInfo &down_seq, const int micro_ho_len, const int arpn)
	{
		up_seq_info = up_seq; down_seq_info = down_seq; microhomology_length = micro_ho_len; abnormal_read_pair_no = arpn;
	}
	OtherInfo& operator= (const OtherInfo &oi) {up_seq_info = oi.up_seq_info; down_seq_info = oi.down_seq_info; microhomology_length = oi.microhomology_length; abnormal_read_pair_no = oi.abnormal_read_pair_no; return *this; }
};

class ClipReads
{
public:
	SeqInfo aligned_seq_info;
	char clipped_side;
	string clipped_seq;
	string clipped_qual;
	//repeat = 'r', none = 'n'
	char type;

	ClipReads(const SeqInfo &aligned, const char c_side, const string &c_seq, const string &c_qual, char t) { aligned_seq_info = aligned; clipped_side = c_side; clipped_seq = c_seq; clipped_qual = c_qual; type = t; }
};

class WholeSeqInfo
{
public:
	SeqInfo up_seq_info;
	SeqInfo down_seq_info;
	int up_pos;
	int down_pos;
	int microhomology_length;
	WholeSeqInfo() { }	
	//SeqInfo& operator= (const SeqInfo &s) { seq = s.seq; cigar = s.cigar; mismatch_no = s.mismatch_no; support_read_no = s.support_read_no; return *this; }
	WholeSeqInfo(const SeqInfo &up_seq_info1, const SeqInfo &down_seq_info1, int up_pos1, int down_pos1, int microhomology_length1){ up_seq_info = up_seq_info1; down_seq_info = down_seq_info1; up_pos = up_pos1; down_pos = down_pos1; microhomology_length = microhomology_length1; }
	void set_value(const SeqInfo &up_seq_info1, const SeqInfo &down_seq_info1, int up_pos1, int down_pos1, int microhomology_length1){ up_seq_info = up_seq_info1; down_seq_info = down_seq_info1; up_pos = up_pos1; down_pos = down_pos1; microhomology_length = microhomology_length1; }

	WholeSeqInfo& operator= (const WholeSeqInfo &wsi) { up_seq_info = wsi.up_seq_info; down_seq_info = wsi.down_seq_info; up_pos = wsi.up_pos; down_pos = wsi.down_pos; microhomology_length = wsi.microhomology_length; return *this; }
};

//class EndDepth
//{
//public:
//	int end_pos;
//	unsigned long sum_depth;
//	EndDepth(int ep, unsigned long sd) { end_pos = ep; sum_depth = sd; }
//};




class Junction
{
public:
	string up_chr;
	int up_pos;
	char up_strand;
	string down_chr;
	int down_pos;
	char down_strand;
	Junction() { }
	Junction(const string u_c, const int u_p, const char u_s, const string d_c, const int d_p, const char d_s) 
	{
		up_chr = u_c;
		up_pos = u_p;
		up_strand = u_s;
		down_chr = d_c;
		down_pos = d_p;
		down_strand = d_s;
	}
	Junction (const Junction &j)
	{
		up_chr = j.up_chr;
		up_pos = j.up_pos;
		up_strand = j.up_strand;
		down_chr = j.down_chr;
		down_pos = j.down_pos;
		down_strand = j.down_strand;
	}

	void set_value(const string &u_c, const int &u_p, const char &u_s, const string &d_c, const int& d_p, const char& d_s) 
	{
		up_chr = u_c;
		up_pos = u_p;
		up_strand = u_s;
		down_chr = d_c;
		down_pos = d_p;
		down_strand = d_s;
	}
	bool operator< (const Junction& junction) const
	{
		if (up_chr < junction.up_chr)
			return 1;
		else if (up_chr > junction.up_chr)
			return 0;
		else
		{
			if (down_chr < junction.down_chr)
				return 1;
			else if (down_chr > junction.down_chr)
				return 0;
			else
			{
				if (up_strand < junction.up_strand)
					return 1;
				else if (up_strand > junction.up_strand)
					return 0;
				else
				{
					if (down_strand < junction.down_strand)
						return 1;
					else if (down_strand > junction.down_strand)
						return 0;
					else
					{
						// same chr , same strand , both upstream breakend and downstream breakend
						if (up_pos < junction.up_pos)
							return 1;
						else if (up_pos > junction.up_pos)
							return 0;
						else if (down_pos < junction.down_pos)
							return 1;
						else return 0;
					}
				}
			}
		}
	}
	friend ostream& operator<< (ostream &os, Junction &nj);
};



class ChrRange
{
public:
	string chr;
	unsigned int begin; //1-based
	unsigned int end; //1-based
	ChrRange(const string &c, const unsigned int b, const unsigned int e) { chr = c; begin = b; end = e; }
	bool operator< (const ChrRange& chr_range) const// it is error without const at the end
	{
		if (chr < chr_range.chr)
			return 1;
		else if (chr == chr_range.chr)
		{
			if (begin < chr_range.begin)
				return 1;
			else if (begin == chr_range.begin)
			{
				if (end < chr_range.end) return 1;
				else return 0;
			}
			else
				return 0;
		}
		else
			return 0;
	}
}
;


/*!@function
  @abstract        Read the alignment results of clipped sequence(bam/sam format) and store them in a map<string, AlingInfo> type object
  @param  file     the path/name of alignment result file
  @param read_id2align_info     the object use to store useful information of alignment results
  @para min_mapQ     minimum map quality of clipped sequences
*/
void StoreSamOfClippedSeq(string file, multimap<string, AlignInfo> &read_id2align_info, int min_mapQ);
void GetAlignInfo(samfile_t *samfin, bam1_t *b, AlignInfo &align_info);

/*@function
 @abstract             1.Read soft-clipped reads file of which postfix is clip. 2.Store the breakpoints 
 @param clipfile            the name of soft-clipped reads file
 @param read_id2align_info     the object use to store useful information of alignment results
 @param aligned2clipped    store the clipped reads which is unmapped or aligned to more than one place
 @param match_rate       match rate of two sequence, match_rate = match_base_number/shorter_length 
*/
template <typename Tinstream> void InputSoftInfoStoreBreakpoint(string clipfile, multimap<string, AlignInfo> &read_id2align_info, multimap<Junction, OtherInfo> &junction2other, multimap<pair<string, int>, ClipReads> &aligned2clipped, double match_rate);
//void InputSoftInfoStoreBreakpoint(string clipfile, map<string, AlignInfo> &read_id2align_info, multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_same_strand, multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_diff_strand_5clip, multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_diff_strand_3clip);


/*@function
 @abstract                          
*/
bool ProcessSeqTandemRepeatUnequal(WholeSeqInfo &small_wseq_info, WholeSeqInfo &big_wseq_info, WholeSeqInfo &result, double match_rate);

/*@function
 @abstract             change an integer to string
 @param       l        the source integer
 @param       s        the target string
*/
void itoa(int l, string &s);

/*@function
  @abstract           change cigar string to cigar vector
  @param     cigar    cigar string
*/
vector<pair<int, char> > ChangeCigarType(const string &cigar);

/*@function
  @abstract          reverse the vector
*/
void ReverseCigar(vector<pair<int, char> > &cigar_vec); 

/*@function
  @abstract          calculate exceed alignment length 
*/
pair<int, int> UpSeqExceedLength(const vector<pair<int, char> > &cigar_vec);
pair<int, int> DownSeqExceedLength(const vector<pair<int, char> > &cigar_vec);


/*@function     
 @abstract                    change the upstream seq,downstream seq, upstream cigar, downstream cigar of a junction
 @param  small_seq_up         upstream seq to be changed 
 @param  small_seq_down       downstream seq to be changed  
 @param  big_seq_up           we will change upstream seq refer big_seq_up
 @param  big_seq_down         we will change downstream seq refer big_seq_down
 @param  small_cigar_up       upstream cigar to be changed 
 @param  small_cigar_down       downstream cigar to be changed 
 @param  microhomology_length
*/
void ChangeString(string &small_seq_up, string &small_seq_down, string big_seq_up, string big_seq_down, string &small_cigar_up, string &small_cigar_down, int microhomology_len);

/*@function
 @abstract                  the function of this function resembles ChangeString, this function is used to the 5 end clipped different strand junction
*/
void ChangeString1(string &small_seq_up, string &small_seq_down, string big_seq_up, string big_seq_down, string &small_cigar_up, string &small_cigar_down, int microhomology_len);

/*@function
 @abstract                  the function of this function resembles ChangeString, this function is used to the 3 end clipped different strand junction
*/
void ChangeString2(string &small_seq_up, string &small_seq_down, string big_seq_up, string big_seq_down, string &small_cigar_up, string &small_cigar_down, int microhomology_len);


/*@function
 @abstract                  the function of this function resembles ChangeString, this function is used to the same strand junction
*/
void ChangeString3(string &small_seq_up, string &big_seq_down, string big_seq_up, string small_seq_down, string &small_cigar_up, string &big_cigar_down, int tandep_repeat_len);

/*function
 @abstract               change cigar in the start end
 @param   cigar          the cigar to be changed
 @param   add_len        the length added to the start end of cigar
*/
void ChangeCigarStartEnd(string &cigar, int add_len);

/*function
 @abstract               change cigar in the back end
 @param   cigar          the cigar to be changed
 @param   add_len        the length added to the back end of cigar
*/
void ChangeCigarBackEnd(string &cigar, int add_len);


/*function
 @abstract               calculate the number of cigar, if cigar == 90M , ncigar = 1, if cigar == 10M10D70M, ncigar = 3;
 @param   cigar          the cigar to be calculate
*/
int NumberCigar(const string &cigar);


/*function
 @abstract               There are some clipped sequence aligned to more than one place, some of these sequences can be used to support the junction,This function is used to recover these sequences (the junctions are of same strand)
*/
void ChangeSameOrientationJuntion(multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_same_strand, multimap<pair<string, int>, ClipReads> &aligned2clipped, int flank, double match_rate);


/*function
 @abstract               There are some clipped sequence aligned to more than one place, some of these sequences can be used to support the junction,This function is used to recover these sequences (the junctions are of different strand)
*/
void ChangeDiffOrientationJuntion(multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_diff_strand_5clip, multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_diff_strand_3clip, multimap<pair<string, int>, ClipReads> &aligned2clipped, int flank, double match_rate);

/*function
 @abstract                     Display bed file of breakpoints. The bed file is use to calculate depth of the breakpoints.
 @param     fout               Ofstream of bed file
*/
void GetBreak(multimap<Junction, OtherInfo> &junction2other, map<pair<string, int>, int> &pos2depth, map<ChrRange , unsigned long> &range2depth, map<Junction, pair<pair<ChrRange, ChrRange>, pair<ChrRange, ChrRange> > > &junction2range_pair, int length);
void GetBreak(multimap<pair<string, int>, ClipReads> &aligned2clipped, map<pair<string, int>, int> &pos2depth);

void MergeOverlap(map<ChrRange, unsigned long> &range2depth, map<pair<string, int>, int> &begin2end);

/*function
 @abstract                                  Display the breakpoint results
 @param   rescure_mode                      Turn off rescue mode. When rescure mode is on, the a SV with only 1 side with enough soft-clipped reads is considered as a valid one instead of rejecting it.  Default on.
 @param   min_no_one_side_clipped_reads     When rescure mode is on, the minimum number of soft-clipped reads on one side [5] 
 @param   sum_min_no_both_clipped_reads     Minimum number of soft clipping read,it's the sum of left clipped reads and right clipped reads [3]
 @param   min_abnormal_read_pair_no         Minimum number of read pairs which support the junction [1]
 @param	  max_microhomology                 Maximum length of tandem repeat,  sv have tandem repeat length longer than [15] will be filtered
 @param   min_seq_len                       Minimum length of up_seq or down_seq near the breakpoint , if abnormal_read_pair_no is larger than 0, this parameter is invalid
 @param   max_seq_indel_no                  Maximum indel number of up_seq or down_seq near the breakpoint, if abnormal_read_pair_no is larger than 0, this parameter is invalid
*/
void OutputBreakpoint(ofstream &fout, multimap<Junction, OtherInfo> &junction2other, map<pair<string, int>, int> &pos2depth, map<ChrRange , unsigned long> &range2depth, map<Junction, pair<pair<ChrRange, ChrRange>, pair<ChrRange, ChrRange> > > &junction2range_pair, bool rescure_mode, int min_no_one_side_clipped_reads, int sum_min_no_both_clipped_reads, int min_abnormal_read_pair_no, int repeat_coverage, double frequency, int min_distance, int max_microhomology, int min_seq_len, int max_seq_indel_no);
void OutputOneendUnmapBreakpoint(ofstream &fout, ofstream &fout1, multimap<pair<string, int>, ClipReads> &aligned2clipped, map<pair<string, int>, int> &pos2depth);


/*function
 @abstract          Find number of discordant read pairs which support the structaral variation breakpoint                     
 @param
*/
void FindDiscordantReadPairs(samfile_t *samfin, bam_index_t *idx, multimap<Junction, OtherInfo> &junction2other, int min_mapQ, int mean_insert_size, int deviation, int times);
int FindDiscordantReadPairs(samfile_t *samfin, bam_index_t *idx, Junction &junction, int min_mapQ, int mean_insert_size, int deviation, int times);



/*function
 @abstract         Read breakpoint file and store the breakpoints in multimap<Junction, OtherInfo> junction2other
*/
void ReadBreakpoint(string file, multimap<Junction, OtherInfo> &junction2other);

double CountLargestBaseFrequency(string &seq);

/*function
 @abstract         Get junction from by synthesize information of clip.gz and clip.bam
*/
void GetJunction(AlignReadsInfo &align_reads_info, char orientation, AlignInfo &clipped_align_info, multimap<Junction, OtherInfo> &junction2other, multimap<pair<string, int>, ClipReads> &aligned2clipped);

void MergeJunction(multimap<Junction, OtherInfo> &junction2other, int search_length);

//New Version
template <typename Tinstream> void InputSoftInfoStoreBreakpoint(string clipfile, string file, multimap<string, AlignInfo> &read_id2align_info, multimap<Junction, OtherInfo> &junction2other, multimap<pair<string, int>, ClipReads> &aligned2clipped, double match_rate)
{
	string chr, read_id, cigar, aligned_seq, aligned_qual, clipped_seq, clipped_qual, temp, last_clipped_seq;
	int support_count, pos;
	char orientation;
	Tinstream fin(clipfile.c_str());
	multimap<string, pair<AlignReadsInfo, char> > clipped2align_reads_info;
	map<pair<string, pair<string, int> >, AlignInfo> read_id_pos2align_info;


//For reading clip.bam
	int len;
	samfile_t *samfin;
	bam1_t *b;
	char in_mode[5], *fn_list = 0;
	in_mode[0] = 'r';
	if (file.rfind(".bam") == file.size() - 4)
	{
		//if in_mode[1] == 'b', it will read a bam file
		in_mode[1] = 'b';
	}

    if ((samfin = samopen(file.c_str(), in_mode, fn_list)) == 0) { cerr << "[main_samview] fail to open file for reading." << endl; delete samfin; exit(1); }
	if (samfin->header == 0) { cerr << "[main_samview] fail to read the header." << endl; exit(1); }

	map<string, pair<string, string> > read_id2seq_pair;

	b = bam_init1();
	int ret = 0;


	while (fin >> chr)
	{
		fin >> pos >> orientation >> cigar >> aligned_seq >> aligned_qual >> clipped_seq >> clipped_qual >> support_count;
		getline(fin, temp);


		vector<pair<int, char> > cigar_vec = ChangeCigarType(cigar);
		ReadsInfo reads_info(aligned_seq, aligned_qual, clipped_seq, clipped_qual, cigar_vec, support_count, 0);
		AlignReadsInfo align_reads_info(chr, pos, reads_info);

		if (last_clipped_seq.empty() || last_clipped_seq == clipped_seq)
		{
			clipped2align_reads_info.insert(make_pair(clipped_seq, make_pair(align_reads_info, orientation)));
			last_clipped_seq = clipped_seq;
			continue;
		}
		else 
		{

			//read a line of the bam file, if the line is without soft-clipping read, ignore it
			while ((ret = samread(samfin, b) >= 0))
			{
				// when read failed continue
				if (__g_skip_aln(samfin->header, b)) continue;
				if (IsHardClip(b)) continue;

				AlignInfo clipped_align_info;
				GetAlignInfo(samfin, b, clipped_align_info);
			

				if (last_clipped_seq == bam1_qname(b))
				{
					read_id_pos2align_info.insert(make_pair(make_pair(last_clipped_seq, make_pair(clipped_align_info.chr, clipped_align_info.pos)), clipped_align_info));
				}
				else
				{
					pair<multimap<string, pair<AlignReadsInfo, char> >::iterator, multimap<string, pair<AlignReadsInfo, char> >::iterator> clipped2align_reads_info_it_pair = clipped2align_reads_info.equal_range(last_clipped_seq);
					map<pair<string, pair<string, int> >, AlignInfo>::iterator read_id_pos2align_info_it = read_id_pos2align_info.begin();
					while (clipped2align_reads_info_it_pair.first != clipped2align_reads_info_it_pair.second)
					{
						while (read_id_pos2align_info_it != read_id_pos2align_info.end())
						{
							GetJunction(clipped2align_reads_info_it_pair.first->second.first, clipped2align_reads_info_it_pair.first->second.second, read_id_pos2align_info_it->second, junction2other, aligned2clipped);
							++ read_id_pos2align_info_it;
						}
						++clipped2align_reads_info_it_pair.first;
					}
					clipped2align_reads_info.clear();
					read_id_pos2align_info.clear();
					clipped2align_reads_info.insert(make_pair(clipped_seq, make_pair(align_reads_info, orientation)));
					read_id_pos2align_info.insert(make_pair(make_pair(last_clipped_seq, make_pair(clipped_align_info.chr, clipped_align_info.pos)), clipped_align_info));
					last_clipped_seq = clipped_seq;
					break;
				}

			}

		}
	}
		//read a line of the bam file, if the line is without soft-clipping read, ignore it
	while ((ret = samread(samfin, b) >= 0))
	{
		// when read failed continue
		if (__g_skip_aln(samfin->header, b)) continue;
		AlignInfo clipped_align_info;
		GetAlignInfo(samfin, b, clipped_align_info);

		if (last_clipped_seq == bam1_qname(b))
		{
			read_id_pos2align_info.insert(make_pair(make_pair(last_clipped_seq, make_pair(clipped_align_info.chr, clipped_align_info.pos)), clipped_align_info));
		}
		else
		{
			break;
		}
	}
	pair<multimap<string, pair<AlignReadsInfo, char> >::iterator, multimap<string, pair<AlignReadsInfo, char> >::iterator> clipped2align_reads_info_it_pair = clipped2align_reads_info.equal_range(last_clipped_seq);
	map<pair<string, pair<string, int> >, AlignInfo>::iterator read_id_pos2align_info_it = read_id_pos2align_info.begin();
	while (clipped2align_reads_info_it_pair.first != clipped2align_reads_info_it_pair.second)
	{
		while (read_id_pos2align_info_it != read_id_pos2align_info.end())
		{
			GetJunction(clipped2align_reads_info_it_pair.first->second.first, clipped2align_reads_info_it_pair.first->second.second, read_id_pos2align_info_it->second, junction2other, aligned2clipped);
			++ read_id_pos2align_info_it;
		}
		++clipped2align_reads_info_it_pair.first;
	}
	clipped2align_reads_info.clear();
	read_id_pos2align_info.clear();
}









// Add in version 1.2.0 , for reading and storing soft-clipped reads from clip.gz file
/*@Function 
  @Abstract  Read *.clip.gz file and store it input breakpos2reads_info
  @param     clipfile          path of *.clip.gz file
  @param     breakpos2reads_info     multimap to store *.clip.gz file
 */
template <typename Tinstream> void InputAndStoreSoftInfo(string clipfile, multimap<pair<string, int>, pair<ReadsInfo, char> > &breakpos2reads_info)
{
	Tinstream fin(clipfile.c_str());
	string chr, cigar, aligned_seq, aligned_qual, clipped_seq, clipped_qual, temp;
	int support_count, pos;
	char orientation;

	while (fin >> chr)
	{
		fin >> pos >> orientation >> cigar >> aligned_seq >> aligned_qual >> clipped_seq >> clipped_qual >> support_count;
		getline(fin, temp);
		vector<pair<int, char> > cigar_vec = ChangeCigarType(cigar);
		if (orientation = '5') {
			ReadsInfo reads_info(clipped_seq, clipped_qual, aligned_seq, aligned_qual, cigar_vec, support_count, 0);
			breakpos2reads_info.insert(make_pair(make_pair(chr, pos) , make_pair(reads_info, orientation)));
		}
		else if (orientation = '3') {
			ReadsInfo reads_info(aligned_seq, aligned_qual, clipped_seq, clipped_qual, cigar_vec, support_count, 0);
			breakpos2reads_info.insert(make_pair(make_pair(chr, pos) , make_pair(reads_info, orientation)));
		}
		else {
			cerr << "Warning : " << chr << '\t' << pos << '\t' << orientation << '\t' << " orientation error" << endl;
		}

	}
}


