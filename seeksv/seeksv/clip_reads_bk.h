/*
 ***************************************
 * soft_reads.h                        *
 *                                     *
 *  Created on:   2011-9-9             *
 *  Modified on:  2012-5-2             *
 *  Modified on:  2012-5-9             *
 *  Modified on:  2013-5-5             *
 *  Version:      0.3.0                *
 *  Author:       qiukl                *
 *                                     *
 ***************************************
 */

#pragma once
#include <ctype.h>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <stdexcept>
#include "head.h"
#include "gzstream.h"
//#include "somatic.h"
//#include "process_bwasw.h"
//#include "sam_info.h"

#define RIGHT_CLIPPED 0
#define LEFT_CLIPPED 1
#define SAME_STRAND 0
#define DIFF_STRAND_RIGHT_CLIP 1
#define DIFF_STRAND_LEFT_CLIP 2
#define KMER_LEN 25

using namespace std;


//typedef  multimap<pair<int, char>, ReadsInfo>::iterator breakpoint2read_iterator;
//change in getclip version 0.1, date: 2011/10/27

//This class is used to store information of the clipped reads
class ReadsInfo
{
public:
	string seq_left;
	string qual_left;
	string seq_right;
	string qual_right;
	vector<pair<int, char> > cigar_vec;
	int support_read_no;
	int used;

	ReadsInfo() { }
	ReadsInfo(string s_l, string q_l, string s_r, string q_r, const vector<pair<int, char> > &cigar, int quantity, int u);
//	~ReadsInfo();
	string get_seq_left() const { return seq_left; }
	string get_qual_left() const { return qual_left; }
	string get_seq_right() const { return seq_right; }
	string get_qual_right() const { return qual_right; }
	vector<pair<int, char> > get_aligned_cigar() { return cigar_vec; }
	//added int version 0.7.0
	int QuantityOfMateread() const { return support_read_no; }
	int get_used() const { return used; }

	ReadsInfo& operator= (const ReadsInfo &reads_info)
	{
		seq_left = reads_info.seq_left;
		qual_left = reads_info.qual_left;
		seq_right = reads_info.seq_right;
		qual_right = reads_info.qual_right;
		cigar_vec = reads_info.cigar_vec;
		support_read_no = reads_info.support_read_no;
		used = reads_info.used;
		return *this;
	}
	
	bool set_used(int u) { used = u; }
	bool set_reads_info(string s_l, string q_l, string s_r, string q_r, const vector<pair<int, char> > &cigar, int quantity);
//	bool add_materead_id(string id);
	bool ChangeSeqAndQual(string s_l, string q_l, string s_r, string q_r, const vector<pair<int, char> > &cigar, bool aa);
	bool support_read_no_increase() { ++support_read_no; }
};

class AlignReadsInfo
{
public:
	string chr;
	int pos;
	ReadsInfo reads_info;
	AlignReadsInfo() { }
	AlignReadsInfo(string c, int p, const ReadsInfo &ri) {chr = c; pos = p; reads_info = ri; }

};



class SeqPair
{
	string seq1;
	string qual1;
	string seq2;
	string qual2;
public:
	bool set_seq1(const string &s1) { seq1 = s1; return 0; }
	bool set_qual1(const string &q1) { qual1 = q1; return 0; }
	bool set_seq2(const string &s2) { seq2 = s2; return 0; }
	bool set_qual2(const string &q2) { qual2 = q2; return 0; }
	string get_seq1() { return seq1; }
	string get_qual1() { return qual1; }
	string get_seq2() { return seq2; }
	string get_qual2() { return qual2; }
	
};


typedef multimap<pair<string, int>, ReadsInfo>::iterator breakpoint2read_iterator;

void GetSClipReads(const bam_header_t *header, const bam1_t *b, multimap<pair<string, int>, ReadsInfo> &breakpoint2read_l, multimap<pair<string, int>, ReadsInfo> &breakpoint2read_r, double limit, int min_mapQ, bool save_low_quality);

double CompareStringEndFirst(string str1, string str2);

double CompareStringBeginFirst(string str1, string str2);


bool IsUsefulSoftClip(const bam1_t *b);
bool InsertSeq(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, string chr, int pos, string seq_left, string qual_left, string seq_right, string qual_right, vector<pair<int, char> > &cigar, string read_id, double limit, bool aa);
//break a sequence into aligned part and clipped part
void GetSeq(const bam1_t *b, int begin_pos, int seq_left_len ,int seq_right_len, string &seq_left, string &qual_left, string &seq_right, string &qual_right, string &read_id);
vector<pair<int, char> > GenerateCigar(const bam1_t *b, int &l);
int Compare(const string seq1, const string seq2, const string seq3, const string seq4, const double match_rate);
bool GetSeqAndQual(const bam1_t *b, string &seq, string &qual);
void GetReverseComplementSeq(string &seq);

bool MinusCigarLeft(vector<pair<int, char> > &cigar_vec, int length);
bool MinusCigarRight(vector<pair<int, char> > &cigar_vec, int length);

bool AddCigarLeft(vector<pair<int, char> > &cigar_vec, int length);
bool AddCigarRight(vector<pair<int, char> > &cigar_vec, int length);

/*function
 @abstract         Get structural variation type
*/
string GetSVType(string up_chr, int up_pos, char up_strand, string down_chr, int down_pos, char down_strand);

/*@funcion
  @abstract    Merge left clipped soft-clipped reads with some alignment deviation
*/
void MergeLeftClippedSoftClippedReads(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, int search_length, double min_match_rate);

/*@funcion
  @abstract    Merge right clipped soft-clipped reads with some alignment deviation
*/
void MergeRightClippedSoftClippedReads(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, int search_length, double min_match_rate);

void GenerateKmer(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator, int> > &kmer2interator, int kmer_len, int left_or_right_clip);
//void GenerateKmer(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, multimap<string, multimap<pair<string, int>, ReadsInfo>::iterator> &kmer2interator, int kmer_len);

void CompareKmer(ofstream &fout, multimap<pair<string, int>, ReadsInfo> &breakpoint2read, multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator, int> > &kmer2interator, int kmer_len, int type);

int FindDistance(string &seq1, string &seq2, int kmer_len);


/*@function
 @abstract         Output cigar vector
*/
template <typename Tout> void DisplayCigarVector(Tout &fout, vector<pair<int, char> > &cigar_vec, int left_length, int right_length);

template <typename Tout> bool StoreUnmapSeqAndQual(const bam1_t *b, const string &seq, const string &qual, map<string, pair<pair<string, string>, char> > &id2seq_qual, Tout &fout1, Tout &fout2)
{
	map<string, pair<pair<string, string>, char> >::iterator map_it = id2seq_qual.find(bam1_qname(b));
	char end;
	if (map_it != id2seq_qual.end())
	{
		if (b->core.flag & BAM_FREAD1)
		{
			if (map_it->second.second == '2')
			{
				fout1 << "@" << map_it->first << "/1" << '\n'
					  << seq << '\n'
					  << "+\n"
					  << qual << endl;

				fout2 << "@" << map_it->first << "/2\n"
					  << map_it->second.first.first << '\n'
					  << "+\n"
					  << map_it->second.first.second << endl; 
				id2seq_qual.erase(bam1_qname(b));
			}
		}
		else
		{
			if (map_it->second.second == '1')
			{
				fout1 << "@" << map_it->first << "/1" << '\n'
					  << map_it->second.first.first << '\n'
					  << "+\n"
					  << map_it->second.first.second << endl; 

				fout2 << "@" << map_it->first << "/2\n"
					  << seq << '\n'
					  << "+\n"
					  << qual << endl;
				id2seq_qual.erase(bam1_qname(b));
			}
		}
	}
	else
	{
		if (b->core.flag & BAM_FREAD1)
			end = '1';
		else
			end = '2';
		id2seq_qual.insert(make_pair(bam1_qname(b), make_pair(make_pair(seq, qual), end)));
	}
}


//output the information of soft-clipping reads and the fq file of soft-clipping reads
template <typename Tout> void DisplaySClipReads(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, map<string, string> &seq2qual, char c,  Tout &softfout)
{
	breakpoint2read_iterator breakpoint2read_it = breakpoint2read.begin();
	map<string, string>::iterator seq2qual_it;
	//string qual;
	if (c == '5')
		while (breakpoint2read_it != breakpoint2read.end())
		{
			//out put the reads of which quality is high enough
			softfout << breakpoint2read_it->first.first << '\t' << breakpoint2read_it->first.second << '\t' << c << '\t';
			vector<pair<int, char> > cigar_vec = breakpoint2read_it->second.get_aligned_cigar();
			DisplayCigarVector<Tout> (softfout, cigar_vec, 0, 0);
			softfout << '\t' << breakpoint2read_it->second.get_seq_right() << '\t' << breakpoint2read_it->second.get_qual_right() << '\t' << breakpoint2read_it->second.get_seq_left() << '\t' << breakpoint2read_it->second.get_qual_left() << '\t' << breakpoint2read_it->second.QuantityOfMateread() << endl;


			if (breakpoint2read_it->second.used == 1)
			{
				breakpoint2read_it++;
				continue;
			}

			seq2qual_it = seq2qual.find(breakpoint2read_it->second.get_seq_left());
			if (seq2qual_it == seq2qual.end())
			{
				seq2qual.insert(make_pair(breakpoint2read_it->second.get_seq_left(), breakpoint2read_it->second.get_qual_left()));
			}
			else
			{
				string qual1 = seq2qual_it->second, qual2 = breakpoint2read_it->second.get_qual_left();

				for (int i = 0; i < qual1.length(); i++)
				{
					if (qual1[i] < qual2[i]) qual1[i] = qual2[i];
				}
				seq2qual_it->second = qual1;
			}
	
			++breakpoint2read_it;
		}
	else if (c == '3')
		while (breakpoint2read_it != breakpoint2read.end())
		{
			//out put the reads of which quality is high enough
			softfout << breakpoint2read_it->first.first << '\t' << breakpoint2read_it->first.second << '\t' << c << '\t';
			vector<pair<int, char> > cigar_vec = breakpoint2read_it->second.get_aligned_cigar();
			DisplayCigarVector<Tout> (softfout, cigar_vec, 0, 0);
			softfout << '\t' << breakpoint2read_it->second.get_seq_left() << '\t' << breakpoint2read_it->second.get_qual_left() << '\t' << breakpoint2read_it->second.get_seq_right() << '\t' << breakpoint2read_it->second.get_qual_right() << '\t' << breakpoint2read_it->second.QuantityOfMateread() << endl;



			if (breakpoint2read_it->second.used == 1)
			{
				breakpoint2read_it++;
				continue;
			}
	
			seq2qual_it = seq2qual.find(breakpoint2read_it->second.get_seq_right());
			if (seq2qual_it == seq2qual.end())
			{
				seq2qual.insert(make_pair(breakpoint2read_it->second.get_seq_right(), breakpoint2read_it->second.get_qual_right()));
			}
			else
			{
				string qual1 = seq2qual_it->second, qual2 = breakpoint2read_it->second.get_qual_right();

				for (int i = 0; i < qual1.length(); i++)
				{
					if (qual1[i] < qual2[i]) qual1[i] = qual2[i];
				}
				seq2qual_it->second = qual1;
			}

			++breakpoint2read_it;
		}
	
}


template <typename Tout> void DisplaySClipFq(map<string, string> &seq2qual, Tout &fout)
{
	map<string, string>::iterator seq2qual_it = seq2qual.begin();
	while (seq2qual_it != seq2qual.end())
	{
		fout << '@' << seq2qual_it->first << '\n'
			 << seq2qual_it->first << '\n'
			 << '+' << '\n'
			 << seq2qual_it->second << endl;
		++seq2qual_it;
	}
}

template <typename Tout> void InputBamOutputReads (string file, multimap<pair<string, int>, ReadsInfo> &breakpoint2read_l, multimap<pair<string, int>, ReadsInfo> &breakpoint2read_r, map<string, pair<pair<string, string>, char> > &id2seq_qual, double match_rate_of_2seq, int min_mapQ, bool save_low_quality, string prefix)
{
/*	tcb needed
	samfile_t *samfin;
	bam1_t *b;
	char in_mode[5], *fn_list = 0;
	in_mode[0] = 'r';
	if (file.rfind(".bam") == file.size() - 4)
	{
		//if in_mode[1] == 'b', it will read a bam file
		in_mode[1] = 'b';
	}

    if ((samfin = samopen(file.c_str(), in_mode, fn_list)) == 0)
	{
		cerr << "[main_samview] fail to open file for reading." << endl;
   		exit(1);
	}
	if (samfin->header == 0) 
	{
		cerr << "[main_samview] fail to read the header." << endl;
		exit(1);
	}

	b = bam_init1();
	int ret = 0;
	//map<string, pair<pair<string, string>, char> >::iterator id2seq_qual_it;
*/
	string soft_clip_reads, clip_seq, pe_mapped_end, pe_unmapped_end, pe_unmapped_cut20_end, pe_unmapped1, pe_unmapped2;
	soft_clip_reads = clip_seq = pe_mapped_end = pe_unmapped_end = pe_unmapped_cut20_end = pe_unmapped1 = pe_unmapped2 = prefix;
	soft_clip_reads.append(".clip.gz");
	clip_seq.append(".clip.fq.gz");
	pe_unmapped1.append(".unmapped_1.fq.gz");
	pe_unmapped2.append(".unmapped_2.fq.gz");
	string temp_breakpoint = prefix;
	temp_breakpoint.append(".temp.breakpoint");
	
	
	Tout softfout(soft_clip_reads.c_str()), fqfout(clip_seq.c_str()), fuout1(pe_unmapped1.c_str()), fuout2(pe_unmapped2.c_str());
	
	if (!softfout) { cerr << "Cannot open file " << soft_clip_reads << endl; exit(1);}
	if (!fqfout) { cerr << "Cannot open file " << clip_seq << endl; exit(1); }
	if (!fuout1) { cerr << "Cannot open file " << pe_unmapped1 << endl; exit(1); }
	if (!fuout2) { cerr << "Cannot open file " << pe_unmapped2 << endl; exit(1); }
/* tcb needed

	//read a line of the bam file, if the line is without soft-clipping read, ignore it
	while ((ret = samread(samfin, b) >= 0))
	{
		// when read failed continue
		if (__g_skip_aln(samfin->header, b)) continue;
		//if the reads is unmap
		if (b->core.flag & BAM_FUNMAP || b->core.flag & BAM_FMUNMAP)
		{
			string seq, qual;
			GetSeqAndQual(b, seq, qual);
			StoreUnmapSeqAndQual<Tout> (b, seq, qual, id2seq_qual, fuout1, fuout2);
		}
		else
			GetSClipReads(samfin->header, b, breakpoint2read_l, breakpoint2read_r, match_rate_of_2seq, min_mapQ, save_low_quality);
	}
*/
//tcb
	string fileclip = prefix;
	fileclip.append(".clip");
	ifstream fin(fileclip.c_str());
	if (!fin)
	{
		cerr << "Cannot open file " << fileclip << endl;
		exit(1);
	}
	string chr, cigar, aligned_seq, aligned_qual, clipped_seq, clipped_qual, temp;
	int pos, read_number;
	char side;
	vector<pair<int, char> > cigar_vec;

	while (fin >> chr)
	{
		fin >> pos >> side >> cigar >> aligned_seq >> aligned_qual >> clipped_seq >> clipped_qual >> read_number;
		getline(fin, temp);
		//	reads_info.add_materead_id(read_id);
		if (side == '5')
		{
			ReadsInfo reads_info(clipped_seq, clipped_qual, aligned_seq, aligned_qual, cigar_vec, read_number, 0);
			breakpoint2read_l.insert(make_pair(make_pair(chr, pos), reads_info));
		}
		else
		{
			ReadsInfo reads_info(aligned_seq, aligned_qual, clipped_seq, clipped_qual, cigar_vec, read_number, 0);
			breakpoint2read_r.insert(make_pair(make_pair(chr, pos), reads_info));
		}

	}

//tce
	cerr << "[GetSClipReads] finished" << endl;
	MergeRightClippedSoftClippedReads(breakpoint2read_r, 50, match_rate_of_2seq);
	cerr << "[MergeRightClippedSoftClippedReads] finished" << endl;
	MergeLeftClippedSoftClippedReads(breakpoint2read_l, 50, match_rate_of_2seq);
	cerr << "[MergeLeftClippedSoftClippedReads] finished" << endl;

	multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator, int> > kmer2interator;

	ofstream svfout(temp_breakpoint.c_str());
	if (!svfout)
	{
		cerr << "Cannot open file " << temp_breakpoint << endl;
		exit(1);
	}
	GenerateKmer(breakpoint2read_r, kmer2interator, KMER_LEN, RIGHT_CLIPPED);
	cerr << "[GenerateKmer] right clipped finished" << endl;

	CompareKmer(svfout, breakpoint2read_l, kmer2interator, KMER_LEN, SAME_STRAND);
	cerr << "[CompareKmer] SAME_STRAND finished" << endl;
	CompareKmer(svfout, breakpoint2read_r, kmer2interator, KMER_LEN, DIFF_STRAND_RIGHT_CLIP);
	cerr << "[CompareKmer] DIFF_STRAND_RIGHT_CLIP finished" << endl;

	



	kmer2interator.clear();
	GenerateKmer(breakpoint2read_l, kmer2interator, KMER_LEN, LEFT_CLIPPED);
	cerr << "[GenerateKmer] left clipped finished" << endl;
	CompareKmer(svfout, breakpoint2read_l, kmer2interator, KMER_LEN, DIFF_STRAND_LEFT_CLIP);
	cerr << "[CompareKmer] DIFF_STRAND_LEFT_CLIP finished" << endl;
	

	map<string, string> seq2qual;
	DisplaySClipReads<Tout> (breakpoint2read_l, seq2qual, '5', softfout);
	DisplaySClipReads<Tout> (breakpoint2read_r, seq2qual, '3', softfout);
	DisplaySClipFq(seq2qual, fqfout);

}

/*@function
 @abstract         Output cigar vector
*/
template <typename Tout> void DisplayCigarVector(Tout &fout, vector<pair<int, char> > &cigar_vec, int left_length, int right_length)
{
	vector<pair<int, char> >::iterator vec_it = cigar_vec.begin();
	if (left_length > 0)
	{
		fout << left_length << 'S';
	}
	while (vec_it != cigar_vec.end())
	{
		fout << (*vec_it).first << (*vec_it).second;
		++vec_it;
	}
	if (right_length > 0)
	{
		fout << right_length << 'S';
	}
}


