/*
 ***************************************
 * getsv.h                             *
 *                                     *
 *  Created on:   2011-5-24            *
 *  Version:      0.2.0                *
 *  seeksv version:      0.8.2         *
 *  Author:       qiukl                *
 *                                     *
 ***************************************
 */
#include "getsv.h"
#include <algorithm>

const int kCrossLength = 5;

ostream& operator<< (ostream &os, Junction &nj)
{
	os << nj.up_chr << '\t' << nj.up_pos << '\t' << nj.up_strand << '\t'
	   << nj.down_chr << '\t' << nj.down_pos << '\t' << nj.down_strand;
	return os;
}


void GetAlignInfo(samfile_t *samfin, bam1_t *b, AlignInfo &align_info)
{
	char type;
	int len;
	string seq;

	//if the reads is unmap
	if (b->core.flag & BAM_FUNMAP)
	{
		int left_clipped_seq_length = 0, right_clipped_seq_length = 0;
		string aligned_seq = bam1_qname(b), seq, qual;
		GetSeqAndQual(b, seq, qual);
		vector<pair<int, char> > cigar_vec;

		align_info.set_value("Exogenous", -1, -1, '*', cigar_vec, seq, left_clipped_seq_length, right_clipped_seq_length, 'n');
	}
	else 
	{
		if (b->core.flag & BAM_FSECONDARY || b->core.qual == 0) 
		{
			type = 'r';
		}
		else type = 'u';
	
		int left_clipped_seq_length = 0, right_clipped_seq_length = 0;
	
		int op1 = bam1_cigar(b)[0] & BAM_CIGAR_MASK, op2 = bam1_cigar(b)[b->core.n_cigar - 1] & BAM_CIGAR_MASK;
		if (op1 == BAM_CSOFT_CLIP || op1 == BAM_CHARD_CLIP) 
		{
			left_clipped_seq_length = (bam1_cigar(b)[0] >> BAM_CIGAR_SHIFT);
		}
		if (op2 == BAM_CSOFT_CLIP || op2 == BAM_CHARD_CLIP) 
		{
			right_clipped_seq_length = (bam1_cigar(b)[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);
		}
	
		seq = bam1_qname(b);
		//Because the sequence may include some deletions, so do not use b->core.l_qseq to count "len"
		vector<pair<int, char> > cigar_vec = GenerateCigar(b, len);
		char strand;
		if (bam1_strand(b))	{ strand = '-';	GetReverseComplementSeq(seq); }
		else { strand = '+'; }
		align_info.set_value(samfin->header->target_name[b->core.tid], b->core.pos + 1, len, strand, cigar_vec, seq, left_clipped_seq_length, right_clipped_seq_length, type);
		//bam1_strand(b) == 1, means its strand is '-'
	}

}

void StoreSamOfClippedSeq(string file, multimap<string, AlignInfo> &read_id2align_info, int min_mapQ)
{
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

	//read a line of the bam file, if the line is without soft-clipping read, ignore it
	while ((ret = samread(samfin, b) >= 0))
	{
		// when read failed continue
		if (__g_skip_aln(samfin->header, b)) continue;
		//if the reads is unmap
		if (b->core.flag & BAM_FUNMAP)
		{
			int left_clipped_seq_length = 0, right_clipped_seq_length = 0;
			string aligned_seq = bam1_qname(b), seq, qual;
			GetSeqAndQual(b, seq, qual);
			vector<pair<int, char> > cigar_vec;

			AlignInfo align_info("Exogenous", -1, -1, '*', cigar_vec, seq, left_clipped_seq_length, right_clipped_seq_length, 'n');
			read_id2align_info.insert(make_pair(aligned_seq, align_info));
		}
		//else if(b->core.qual >= min_mapQ || (b->core.qual & BAM_FSECONDARY))
		else //if(!(b->core.qual & BAM_FSECONDARY))
		{
			char type;
			string seq, qual, seq1;
			if (b->core.flag & BAM_FSECONDARY) 
			{
				if (read_id2seq_pair.count(bam1_qname(b)) == 0) continue;
				else type = 'r';
			}
			else if (b->core.qual < min_mapQ)
			{
				GetSeqAndQual(b, seq, qual);
				seq1 = seq;
				GetReverseComplementSeq(seq1);
				// strand == '-'
				if (bam1_strand(b)) read_id2seq_pair.insert(make_pair(bam1_qname(b), make_pair(seq1, seq)));
				else read_id2seq_pair.insert(make_pair(bam1_qname(b), make_pair(seq, seq1)));
				type = 'r';
			}
			else type = 'u';


			//uint8_t *s = bam_aux_get(b, "XM");
		//	int mismatch_no = bam_aux2i(s);

			GetSeqAndQual(b, seq, qual);
			int left_clipped_seq_length = 0, right_clipped_seq_length = 0;

			int op1 = bam1_cigar(b)[0] & BAM_CIGAR_MASK, op2 = bam1_cigar(b)[b->core.n_cigar - 1] & BAM_CIGAR_MASK;
			if (op1 == BAM_CSOFT_CLIP) 
			{
				left_clipped_seq_length = (bam1_cigar(b)[0] >> BAM_CIGAR_SHIFT);
			}
			if (op2 == BAM_CSOFT_CLIP) 
			{
				right_clipped_seq_length = (bam1_cigar(b)[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);
			}

			//Because the sequence may include some deletions, so do not use b->core.l_qseq to count "len"
			vector<pair<int, char> > cigar_vec = GenerateCigar(b, len);
			char strand;
			if (bam1_strand(b))
			{
				strand = '-';
				if (b->core.flag & BAM_FSECONDARY && read_id2seq_pair.count(bam1_qname(b)) == 1) seq = read_id2seq_pair[bam1_qname(b)].second;
			}
			else
			{
				strand = '+';
				if (b->core.flag & BAM_FSECONDARY && read_id2seq_pair.count(bam1_qname(b)) == 1) seq = read_id2seq_pair[bam1_qname(b)].first;
			}


			AlignInfo align_info(samfin->header->target_name[b->core.tid], b->core.pos + 1, len, strand, cigar_vec, seq, left_clipped_seq_length, right_clipped_seq_length, type);
			//bam1_strand(b) == 1, means its strand is '-'
			read_id2align_info.insert(make_pair(bam1_qname(b), align_info));


		}
	}
}


/*
bool ProcessSeqTandemRepeatUnequal(WholeSeqInfo &small_wseq_info, WholeSeqInfo &big_wseq_info, WholeSeqInfo &result, double match_rate)
{
	
	if (NumberCigar(small_wseq_info.up_seq_info.cigar) == 1 && NumberCigar(small_wseq_info.down_seq_info.cigar) == 1)
	{
		int microhomology_len = -1;
		if (NumberCigar(big_wseq_info.up_seq_info.cigar) == 1 && NumberCigar(big_wseq_info.down_seq_info.cigar) > 1)
			microhomology_len = big_wseq_info.up_pos - small_wseq_info.up_pos;
		else if (NumberCigar(big_wseq_info.up_seq_info.cigar) > 1 && NumberCigar(big_wseq_info.down_seq_info.cigar) == 1)
			microhomology_len = big_wseq_info.down_pos - small_wseq_info.down_pos;
		if (microhomology_len == -1) { return 1; }
		if (big_wseq_info.up_seq_info.seq.length() <= microhomology_len) { return 1; }
		string seq1 = big_wseq_info.up_seq_info.seq;
		string seq2 = seq1.substr(seq1.length() - microhomology_len);
		seq2.append(big_wseq_info.down_seq_info.seq);
		seq1 = seq1.substr(0, seq1.length() - microhomology_len);
		
		if (CompareStringEndFirst(small_wseq_info.up_seq_info.seq, seq1) >= match_rate && CompareStringBeginFirst(small_wseq_info.down_seq_info.seq, seq2) >= match_rate)
		{
			result = small_wseq_info;
//			ChangeString(result.up_seq, result.down_seq, big_wseq_info.up_seq, big_wseq_info.down_seq, result.up_cigar, result.down_cigar, microhomology_len);
			result.up_seq_info.support_read_no = big_wseq_info.up_seq_info.support_read_no;
			result.microhomology_length = 0;
			return 0;
		}
		else return 1;
	}
	else if (NumberCigar(big_wseq_info.up_seq_info.cigar) == 1 && NumberCigar(big_wseq_info.down_seq_info.cigar) == 1)
	{
		int microhomology_len = -1;
		if (NumberCigar(small_wseq_info.up_seq_info.cigar) == 1 && NumberCigar(small_wseq_info.down_seq_info.cigar) > 1)
			microhomology_len = big_wseq_info.up_pos - small_wseq_info.up_pos;
		else if (NumberCigar(small_wseq_info.up_seq_info.cigar) > 1 && NumberCigar(small_wseq_info.down_seq_info.cigar) == 1)
			microhomology_len = big_wseq_info.down_pos - small_wseq_info.down_pos;
		if (microhomology_len == -1) { return 1;}
		if (big_wseq_info.up_seq_info.seq.length() <= microhomology_len) { return 1;}

		string seq1 = big_wseq_info.up_seq_info.seq;
		string seq2 = seq1.substr(seq1.length() - microhomology_len);
		seq2.append(big_wseq_info.down_seq_info.seq);
		seq1 = seq1.substr(0, seq1.length() - microhomology_len);
		if (CompareStringEndFirst(small_wseq_info.up_seq_info.seq, seq1) >= match_rate && CompareStringBeginFirst(small_wseq_info.down_seq_info.seq, seq2) >= match_rate)
		{

			result = big_wseq_info;
	//		result.up_seq = small_wseq_info.up_seq;
	//		result.up_cigar = small_wseq_info.up_cigar;
	//		result.up_pos = small_wseq_info.up_pos;

			//ChangeString3(result.up_seq, result.down_seq, big_wseq_info.up_seq, small_wseq_info.down_seq, result.up_cigar, result.down_cigar, microhomology_len);
			result.down_seq_info.support_read_no = small_wseq_info.down_seq_info.support_read_no;
			//result.down_pos = result.down_pos - microhomology_len;
			result.microhomology_length = 0;
			return 0;
		}
		else return 1;
	}
	else return 1;
}


*/



void ChangeString(string &small_seq_up, string &small_seq_down, string big_seq_up, string big_seq_down, string &small_cigar_up, string &small_cigar_down, int microhomology_len)
{
	if (big_seq_up.length() > microhomology_len + small_seq_up.length())
	{
		int add_len = big_seq_up.length() - microhomology_len - small_seq_up.length();
		string temp = small_seq_up;
		small_seq_up = big_seq_up.substr(0, add_len);
		small_seq_up.append(temp);
		ChangeCigarStartEnd(small_cigar_up, add_len);
	}
	if (big_seq_down.length() + microhomology_len > small_seq_down.length())
	{
		int add_len = big_seq_down.length() + microhomology_len - small_seq_down.length();
		if (big_seq_down.length() >= add_len)
		{
			small_seq_down.append(big_seq_down, big_seq_down.length() - add_len, add_len);
			ChangeCigarBackEnd(small_cigar_down, add_len);
		}
		else if (big_seq_up.length() >= microhomology_len)
		{
			string temp = big_seq_up.substr(big_seq_up.length() - microhomology_len);
			temp = temp.append(big_seq_down);
			small_seq_down.append(temp, temp.length() - add_len, add_len);
			ChangeCigarBackEnd(small_cigar_down, add_len);
		}
	}
}


void ChangeString1(string &small_seq_up, string &small_seq_down, string big_seq_up, string big_seq_down, string &small_cigar_up, string &small_cigar_down, int microhomology_len)
{
	if (big_seq_up.length() + microhomology_len > small_seq_up.length())
	{
		int add_len = big_seq_up.length() + microhomology_len - small_seq_up.length();
		string temp = small_seq_up;
		small_seq_up = big_seq_up.substr(0, add_len);
		small_seq_up.append(temp);
		ChangeCigarBackEnd(small_cigar_up, add_len);
	}
	if (big_seq_down.length() > microhomology_len + small_seq_down.length())
	{
		int add_len = big_seq_down.length() - microhomology_len - small_seq_down.length();
		if (big_seq_down.length() >= add_len)
		{
			small_seq_down.append(big_seq_down, big_seq_down.length() - add_len, add_len);
			ChangeCigarStartEnd(small_cigar_down, add_len);
		}
		else if (big_seq_up.length() >= microhomology_len)
		{
			string temp = big_seq_up.substr(big_seq_up.length() - microhomology_len);
			temp = temp.append(big_seq_down);
			small_seq_down.append(temp, temp.length() - add_len, add_len);
			ChangeCigarBackEnd(small_cigar_down, add_len);
		}
	}
}


void ChangeString2(string &small_seq_up, string &small_seq_down, string big_seq_up, string big_seq_down, string &small_cigar_up, string &small_cigar_down, int microhomology_len)
{
	if (big_seq_up.length() > microhomology_len + small_seq_up.length())
	{
		int add_len = big_seq_up.length() - microhomology_len - small_seq_up.length();
		string temp = small_seq_up;
		small_seq_up = big_seq_up.substr(0, add_len);
		small_seq_up.append(temp);
		ChangeCigarStartEnd(small_cigar_up, add_len);
	}
	if (big_seq_down.length() + microhomology_len > small_seq_down.length())
	{
		int add_len = big_seq_down.length() + microhomology_len - small_seq_down.length();
		if (big_seq_down.length() >= add_len)
		{
			small_seq_down.append(big_seq_down, big_seq_down.length() - add_len, add_len);
			ChangeCigarStartEnd(small_cigar_down, add_len);
		}
		else if (big_seq_up.length() >= microhomology_len)
		{
			string temp = big_seq_up.substr(big_seq_up.length() - microhomology_len);
			temp = temp.append(big_seq_down);
			small_seq_down.append(temp, temp.length() - add_len, add_len);
			ChangeCigarBackEnd(small_cigar_down, add_len);
		}
	}
}

void ChangeString3(string &small_seq_up, string &big_seq_down, string big_seq_up, string small_seq_down, string &small_cigar_up, string &big_cigar_down, int microhomology_len)
{
	if (big_seq_up.length() > microhomology_len + small_seq_up.length())
	{
		int add_len = big_seq_up.length() - microhomology_len - small_seq_up.length();
		string temp = small_seq_up;
		small_seq_up = big_seq_up.substr(0, add_len);
		small_seq_up.append(temp);
		ChangeCigarStartEnd(small_cigar_up, add_len);
	}
	if (small_seq_down.length() > microhomology_len + big_seq_down.length())
	{
		int add_len = small_seq_down.length() - microhomology_len - big_seq_down.length();
		big_seq_down.append(small_seq_down, small_seq_down.length() - add_len, add_len);
		ChangeCigarBackEnd(big_cigar_down, add_len);
	}
	string temp = big_seq_up.substr(big_seq_up.length() - microhomology_len);
	big_seq_down = temp.append(big_seq_down);
	ChangeCigarStartEnd(big_cigar_down, microhomology_len);
}
/*
void ChangeString3(string &big_seq_up, string &big_seq_down, string small_seq_up, string small_seq_down, string &big_cigar_up, string &big_cigar_down, int microhomology_len)
{
	if (small_seq_up.length() + microhomology_len > big_seq_up.length())
	{
		int add_len = small_seq_up.length() + microhomology_len - big_seq_up.length();
		string temp = big_seq_up;
		big_seq_up = small_seq_up.substr(0, add_len);
		big_seq_up.append(temp);
		ChangeCigarStartEnd(big_cigar_up, add_len);
	}
	if (small_seq_down.length() > microhomology_len + big_seq_down.length())
	{
		int add_len = small_seq_down.length() - microhomology_len - big_seq_down.length();
		big_seq_up.append(small_seq_down, small_seq_down.length() - add_len, add_len);
		ChangeCigarStartEnd(big_cigar_down, add_len);
	}

}
*/

void ChangeCigarStartEnd(string &cigar, int add_len)
{
	int l = 0, i;
	string temp;
	for (i = 0; i < cigar.length(); i++)
	{
		if (isdigit(cigar[i]))
		{
			l = l * 10 + (cigar[i] - 48);
		}
		else
			break;
	}
	l += add_len;
	string str_new_len;
	itoa(l, str_new_len);
	temp = cigar;
	cigar = str_new_len;
	cigar.append(temp, i, temp.length() - i);
}


void ChangeCigarBackEnd(string &cigar, int add_len)
{
	int l, i;
	char m = cigar[cigar.length() - 1];
	string temp, str_new_len;
	for (i = 2; i <= cigar.length(); ++i)
	{
		if (!isdigit(cigar[cigar.length() - i]))
			break;
	}
	temp = cigar.substr(cigar.length() - i + 1);
	l = atoi(temp.c_str());
	l += add_len;
	itoa(l, str_new_len); 
	cigar = cigar.substr(0, cigar.length() - i + 1);
	cigar.append(str_new_len);
	cigar.append(1, m);
}


int NumberCigar(const string &cigar)
{
	int ncigar = 0;
	for (int i = 0; i < cigar.length(); ++i)
		if (!isdigit(cigar[i])) ++ncigar;
	return ncigar;
}

void itoa(int l, string &s)
{
	int i = 0;
	char a[15];
	while (l > 0)
	{
		a[i] = l % 10 + 48;
		l = l / 10;
		++i;
	}
	a[i] = 0;
	s = a;
	std::reverse(s.begin(), s.end());
}


vector<pair<int, char> > ChangeCigarType(const string &cigar)
{
	vector<pair<int, char> > cigar_vec;
	int cigar_len = cigar.length();
	int length = 0;
	for (int i = 0; i < cigar_len; ++i)
	{
		if (isdigit(cigar[i]))
		{
			length = length * 10 + (cigar[i] - 48);
		}
		else
		{
			cigar_vec.push_back(make_pair(length, cigar[i]));
			length = 0;
		}
	}
	return cigar_vec;
}

void ReverseCigar(vector<pair<int, char> > &cigar_vec)
{
	pair<int, char> len_type_pair;
	int size = cigar_vec.size();
	for (int i = 0; i < size/2; i++)
	{
		len_type_pair = cigar_vec[i];
		cigar_vec[i] = cigar_vec[size-i-1];
		cigar_vec[size-i-1] = len_type_pair;
	}
}

pair<int, int> UpSeqExceedLength(const vector<pair<int, char> > &cigar_vec)
{
	int size = cigar_vec.size();
	int len = cigar_vec[size-1].first;
	pair<int, int> len_pair = make_pair(len, len);
	if (cigar_vec[size-2].second == 'I')
	{
		len_pair.second += cigar_vec[size-2].first; 
	}
	else if (cigar_vec[size-2].second == 'D')
	{
		len_pair.first += cigar_vec[size-2].first;
	}
	return len_pair;
}


pair<int, int> DownSeqExceedLength(const vector<pair<int, char> > &cigar_vec)
{
	int len = cigar_vec[0].first;
	pair<int, int> len_pair = make_pair(len, len);
	if (cigar_vec[1].second == 'I')
	{
		len_pair.first += cigar_vec[1].first; 
	}
	else if (cigar_vec[1].second == 'D')
	{
		len_pair.second += cigar_vec[1].first;
	}
	return len_pair;
}

/*
void ChangeSameOrientationJuntion(multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_same_strand, multimap<pair<string, int>, ClipReads> &aligned2clipped, int flank, double match_rate)
{
	multimap<pair<string, int>, pair<string, int> > down_pos2up_pos;

	multimap<pair<string, int>, Breakpoint>::iterator mmap_it = uppos2breakpoint_same_strand.begin();
	while (mmap_it != uppos2breakpoint_same_strand.end())
	{
		if (mmap_it->second.microhomology_length == -1)
			down_pos2up_pos.insert(make_pair(make_pair(mmap_it->second.down_chr, mmap_it->second.down_pos), mmap_it->first));
		++mmap_it;
	}

	multimap<pair<string, int>, pair<string, int> >::iterator mmap_it2;
	multimap<pair<string, int>, ClipReads>::iterator mmap_it1  = aligned2clipped.begin();
	pair<multimap<pair<string, int>, Breakpoint>::iterator, multimap<pair<string, int>, Breakpoint>::iterator> iter_pair;
	while (mmap_it1 != aligned2clipped.end())
	{
		switch (mmap_it1->second.clipped_side)
		{
		case '5':
		{
			mmap_it2 = down_pos2up_pos.lower_bound(mmap_it1->first);
			while (mmap_it2 != down_pos2up_pos.end() && mmap_it2->first.first == mmap_it1->first.first && mmap_it2->first.second <= mmap_it1->first.second + flank)
			{
				//int microhomology_len = mmap_it2->first.second - mmap_it1->first.second;
				//int microhomology_len = Compare(
				iter_pair = uppos2breakpoint_same_strand.equal_range(mmap_it2->second);
				while (iter_pair.first != iter_pair.second)
				{
					if (make_pair(iter_pair.first->second.down_chr, iter_pair.first->second.down_pos)== mmap_it2->first)
					{
						int microhomology_len = Compare(iter_pair.first->second.up_seq_info.seq, iter_pair.first->second.down_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it1->second.aligned_seq_info.seq, match_rate);
						if (microhomology_len != -1)
						{
							int i, l;
							string temp, str_new_len, down_cigar = mmap_it1->second.aligned_seq_info.cigar, up_cigar = iter_pair.first->second.up_seq_info.cigar;
							if (up_cigar[up_cigar.length() - 1] != 'M') break;
							for (i = 2; i <= up_cigar.length(); ++i)
							{
							    if (!isdigit(up_cigar[up_cigar.length() - i]))
							         break;
							}
							temp = up_cigar.substr(up_cigar.length() - i + 1);
							l = atoi(temp.c_str());
							if (l <= microhomology_len) break;
							itoa(l, str_new_len);
							up_cigar = up_cigar.substr(0, up_cigar.length() - i + 1);
							up_cigar.append(str_new_len);
							up_cigar.append(1, 'M');

							string up_seq = iter_pair.first->second.up_seq_info.seq, down_seq = iter_pair.first->second.down_seq_info.seq;
							if (mmap_it1->second.clipped_seq.length() + microhomology_len >= up_seq.length())
								up_seq = mmap_it1->second.clipped_seq;
							else
								up_seq = up_seq.substr(0, up_seq.length() - microhomology_len);
							if (down_seq.length() + microhomology_len <= mmap_it1->second.aligned_seq_info.seq.length())
								down_seq = mmap_it1->second.aligned_seq_info.seq;
							else
							{
								ChangeCigarBackEnd(down_cigar, down_seq.length() + microhomology_len - mmap_it1->second.aligned_seq_info.seq.length());
								temp = down_seq;
								down_seq = mmap_it1->second.aligned_seq_info.seq.substr(0, microhomology_len);
								down_seq.append(temp);
							}
							
							SeqInfo up_seq_info(up_seq, up_cigar, iter_pair.first->second.up_seq_info.support_read_no), down_seq_info(down_seq, down_cigar, mmap_it1->second.aligned_seq_info.support_read_no);
							Breakpoint breakpoint(iter_pair.first->second.down_chr, iter_pair.first->second.down_pos - microhomology_len, microhomology_len, up_seq_info, down_seq_info);
							string up_chr = iter_pair.first->first.first;
							int up_pos = iter_pair.first->first.second - microhomology_len;
							uppos2breakpoint_same_strand.erase(iter_pair.first);
							uppos2breakpoint_same_strand.insert(make_pair(make_pair(up_chr, up_pos), breakpoint));
						//	else
						//	{
						//	iter_pair.first->second.down_seq_info.support_read_no = mmap_it1->second.aligned_seq_info.support_read_no;
						//	iter_pair.first->second.microhomology_length = microhomology_len;
						//	mmap_it1->second.type = mmap_it1->second.type - 32;	
						//	}
						//	
							break;
						}
					}
					iter_pair.first++;
				}
				++mmap_it2;
			}
			break;
		}
		case '3':
		{
			mmap_it = uppos2breakpoint_same_strand.lower_bound(make_pair(mmap_it1->first.first, mmap_it1->first.second - flank));
			while (mmap_it != uppos2breakpoint_same_strand.end() && mmap_it->first.first == mmap_it1->first.first && mmap_it->first.second <= mmap_it1->first.second)
			{
				if (mmap_it->second.microhomology_length == -1 && mmap_it->second.up_seq_info.support_read_no == 0)
				{
					//judge whether there are come from the same junction(breakpoint)
					int microhomology_len = Compare(mmap_it1->second.aligned_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it->second.up_seq_info.seq, mmap_it->second.down_seq_info.seq, match_rate); 
					if (microhomology_len != -1)
					{
						ChangeString(mmap_it->second.up_seq_info.seq, mmap_it->second.down_seq_info.seq, mmap_it1->second.aligned_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it->second.up_seq_info.cigar, mmap_it->second.down_seq_info.cigar, microhomology_len);
						mmap_it->second.up_seq_info.support_read_no = mmap_it1->second.aligned_seq_info.support_read_no;
						mmap_it->second.microhomology_length = microhomology_len;
						mmap_it1->second.type = mmap_it1->second.type - 32;	
						break;
					}

				}
				++mmap_it;
			}
			break;
		}
		default:
			break;
		}
		++mmap_it1;
	}
}


void ChangeDiffOrientationJuntion(multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_diff_strand_5clip, multimap<pair<string, int>, Breakpoint> &uppos2breakpoint_diff_strand_3clip, multimap<pair<string, int>, ClipReads> &aligned2clipped, int flank, double match_rate)
{
	multimap<pair<string, int>, pair<string, int> > down_pos2up_pos, down_pos2up_pos1;
	multimap<pair<string, int>, pair<string, int> >::iterator mmap_it2;
	multimap<pair<string, int>, Breakpoint>::iterator mmap_it;
	multimap<pair<string, int>, ClipReads>::iterator mmap_it1  = aligned2clipped.begin();
	pair<multimap<pair<string, int>, Breakpoint>::iterator, multimap<pair<string, int>, Breakpoint>::iterator> iter_pair;

	mmap_it = uppos2breakpoint_diff_strand_5clip.begin();
	while (mmap_it != uppos2breakpoint_diff_strand_5clip.end())
	{
		if (mmap_it->second.microhomology_length == -1)
			down_pos2up_pos.insert(make_pair(make_pair(mmap_it->second.down_chr, mmap_it->second.down_pos), mmap_it->first));
		++mmap_it;
	}

	mmap_it = uppos2breakpoint_diff_strand_3clip.begin();
	while (mmap_it != uppos2breakpoint_diff_strand_3clip.end())
	{
		if (mmap_it->second.microhomology_length == -1)
			down_pos2up_pos1.insert(make_pair(make_pair(mmap_it->second.down_chr, mmap_it->second.down_pos), mmap_it->first));
		++mmap_it;
	}


	while (mmap_it1 != aligned2clipped.end())
	{
		switch (mmap_it1->second.clipped_side)
		{
		case '5':
		{
			mmap_it = uppos2breakpoint_diff_strand_5clip.lower_bound(mmap_it1->first);
			while (mmap_it != uppos2breakpoint_diff_strand_5clip.end() && mmap_it->first.first == mmap_it1->first.first && mmap_it->first.second <= mmap_it1->first.second + flank)
			{
				if (mmap_it->second.microhomology_length == -1 && mmap_it->second.up_seq_info.support_read_no == 0)
				{
					//judge whether there are come from the same junction(breakpoint)
					int microhomology_len = Compare(mmap_it1->second.aligned_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it->second.up_seq_info.seq, mmap_it->second.down_seq_info.seq, match_rate); 
					if (microhomology_len != -1)
					{
						mmap_it->second.up_seq_info.support_read_no = mmap_it1->second.aligned_seq_info.support_read_no;
						mmap_it->second.microhomology_length = microhomology_len;
						mmap_it1->second.type = mmap_it1->second.type - 32;	
						break;
					}

				}
				++mmap_it;
			}
			if (mmap_it == uppos2breakpoint_diff_strand_5clip.end() || mmap_it->first.first != mmap_it1->first.first || mmap_it->first.second > mmap_it1->first.second + flank)
			{
				string aligned_seq = mmap_it1->second.aligned_seq_info.seq, clipped_seq = mmap_it1->second.clipped_seq;
				GetReverseComplementSeq(aligned_seq);
				GetReverseComplementSeq(clipped_seq);
				mmap_it2 = down_pos2up_pos.lower_bound(mmap_it1->first);
				while (mmap_it2 != down_pos2up_pos.end() && mmap_it2->first.first == mmap_it1->first.first && mmap_it2->first.second <= mmap_it1->first.second + flank)
				{
					iter_pair = uppos2breakpoint_diff_strand_5clip.equal_range(mmap_it2->second);
					while (iter_pair.first != iter_pair.second)
					{
						if (make_pair(iter_pair.first->second.down_chr, iter_pair.first->second.down_pos)== mmap_it2->first)
						{
							int microhomology_len = Compare(iter_pair.first->second.up_seq_info.seq, iter_pair.first->second.down_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it1->second.aligned_seq_info.seq, match_rate);
							if (microhomology_len != -1)
							{
								ChangeString1(iter_pair.first->second.up_seq_info.seq, iter_pair.first->second.down_seq_info.seq, clipped_seq, aligned_seq, iter_pair.first->second.up_seq_info.cigar, iter_pair.first->second.down_seq_info.cigar, microhomology_len);
								iter_pair.first->second.down_seq_info.support_read_no = mmap_it1->second.aligned_seq_info.support_read_no;
								iter_pair.first->second.microhomology_length = microhomology_len;
								mmap_it1->second.type = mmap_it1->second.type - 32;	
								break;
							}
						}
						iter_pair.first++;
					}
					++mmap_it2;
				}
			}
			break;
		}
		case '3':
		{
			mmap_it = uppos2breakpoint_diff_strand_3clip.lower_bound(make_pair(mmap_it1->first.first, mmap_it1->first.second - flank));
			while (mmap_it != uppos2breakpoint_diff_strand_3clip.end() && mmap_it->first.first == mmap_it1->first.first && mmap_it->first.second <= mmap_it1->first.second)
			{
				if (mmap_it->second.microhomology_length == -1 && mmap_it->second.up_seq_info.support_read_no == 0)
				{
					//judge whether there are come from the same junction(breakpoint)
					int microhomology_len = Compare(mmap_it1->second.aligned_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it->second.up_seq_info.seq, mmap_it->second.down_seq_info.seq, match_rate); 
					if (microhomology_len != -1)
					{
						ChangeString2(mmap_it->second.up_seq_info.seq, mmap_it->second.down_seq_info.seq, mmap_it1->second.aligned_seq_info.seq, mmap_it1->second.clipped_seq, mmap_it->second.up_seq_info.cigar, mmap_it->second.down_seq_info.cigar, microhomology_len);
						mmap_it->second.up_seq_info.support_read_no = mmap_it1->second.aligned_seq_info.support_read_no;
						mmap_it->second.microhomology_length = microhomology_len;
						mmap_it1->second.type = mmap_it1->second.type - 32;	
						break;
					}

				}
				++mmap_it;
			}
			if (mmap_it == uppos2breakpoint_diff_strand_3clip.end() || mmap_it->first.first != mmap_it1->first.first || mmap_it->first.second > mmap_it1->first.second)
			{
				string aligned_seq = mmap_it1->second.aligned_seq_info.seq, clipped_seq = mmap_it1->second.clipped_seq;
				GetReverseComplementSeq(aligned_seq);
				GetReverseComplementSeq(clipped_seq);
				mmap_it2 = down_pos2up_pos1.lower_bound(make_pair(mmap_it1->first.first, mmap_it1->first.second - flank));
				while (mmap_it2 != down_pos2up_pos1.end() && mmap_it2->first.first == mmap_it1->first.first && mmap_it2->first.second <= mmap_it1->first.second)
				{
					iter_pair = uppos2breakpoint_diff_strand_3clip.equal_range(mmap_it2->second);
					while (iter_pair.first != iter_pair.second)
					{
						if (make_pair(iter_pair.first->second.down_chr, iter_pair.first->second.down_pos)== mmap_it2->first)
						{
							int microhomology_len = Compare(iter_pair.first->second.up_seq_info.seq, iter_pair.first->second.down_seq_info.seq, clipped_seq, aligned_seq, match_rate);
							if (microhomology_len != -1)
							{
								iter_pair.first->second.down_seq_info.support_read_no = mmap_it1->second.aligned_seq_info.support_read_no;
								iter_pair.first->second.microhomology_length = microhomology_len;
								mmap_it1->second.type = mmap_it1->second.type - 32;	
								break;
							}
						}
						iter_pair.first++;
					}
					++mmap_it2;
				}
			}
			break;
		}
		default:
			break;
		}
		++mmap_it1;
	}
}
*/
void GetBreak(multimap<Junction, OtherInfo> &junction2other, map<pair<string, int>, int> &pos2depth, map<ChrRange , unsigned long> &range2depth, map<Junction, pair<pair<ChrRange, ChrRange>, pair<ChrRange, ChrRange> > > &junction2range_pair, int length)
{
	multimap<Junction, OtherInfo>::iterator mmap_it = junction2other.begin();
	while (mmap_it != junction2other.end())
	{
		pos2depth.insert(make_pair(make_pair(mmap_it->first.up_chr, mmap_it->first.up_pos), 0));
		pos2depth.insert(make_pair(make_pair(mmap_it->first.down_chr, mmap_it->first.down_pos), 0));

		unsigned int up_up_begin, up_up_end, up_down_begin, up_down_end, down_up_begin, down_up_end, down_down_begin, down_down_end;
		int l = 0;
		if (mmap_it->first.up_chr == mmap_it->first.down_chr && mmap_it->first.up_strand == mmap_it->first.down_strand)
		{
			abs(mmap_it->first.down_pos - 1 - mmap_it->first.up_pos) < length ? l = abs(mmap_it->first.down_pos - 1 - mmap_it->first.up_pos) : l = length;
		}
		else
		{
			l = length;
		}
		up_up_begin = mmap_it->first.up_pos - l + 1;
		up_up_end = mmap_it->first.up_pos;
		up_down_begin = mmap_it->first.up_pos + 1;
		up_down_end = mmap_it->first.up_pos + l;
		down_up_begin = mmap_it->first.down_pos - l;
		down_up_end = mmap_it->first.down_pos - 1;
		down_down_begin = mmap_it->first.down_pos;
		down_down_end = mmap_it->first.down_pos + l - 1;

		ChrRange chr_range1(mmap_it->first.up_chr, up_up_begin, up_up_end), chr_range2(mmap_it->first.up_chr, up_down_begin, up_down_end), chr_range3(mmap_it->first.down_chr, down_up_begin, down_up_end), chr_range4(mmap_it->first.down_chr, down_down_begin, down_down_end);	
		range2depth.insert(make_pair(chr_range1, 0));
		range2depth.insert(make_pair(chr_range2, 0));
		range2depth.insert(make_pair(chr_range3, 0));
		range2depth.insert(make_pair(chr_range4, 0));

		junction2range_pair.insert(make_pair(mmap_it->first, make_pair(make_pair(chr_range1, chr_range2), make_pair(chr_range3, chr_range4))));

		++mmap_it;
	}
}

void GetBreak(multimap<pair<string, int>, ClipReads> &aligned2clipped, map<pair<string, int>, int> &pos2depth)
{
	multimap<pair<string, int>, ClipReads>::iterator mmap_it1 = aligned2clipped.begin();
	while (mmap_it1 != aligned2clipped.end())
	{
		if (mmap_it1->second.type == 'n')
		{
			pos2depth.insert(make_pair(make_pair(mmap_it1->first.first, mmap_it1->first.second), 0));
		}
		++mmap_it1;
	}
}

void MergeOverlap(map<ChrRange, unsigned long> &range2depth, map<pair<string, int>, int> &begin2end)
{
	string chr;
	unsigned int begin, end;
	map<ChrRange, unsigned long>::iterator map_it = range2depth.begin();
	while (map_it != range2depth.end())
	{
		if (map_it == range2depth.begin())
		{
			chr = map_it->first.chr;
			begin = map_it->first.begin;
			end = map_it->first.end;
		}
		else
		{
			if (chr == map_it->first.chr && begin <= map_it->first.begin && end + 1 >= map_it->first.begin)
			{
				map_it->first.end > end && (end = map_it->first.end);
			}
			else
			{
				begin2end.insert(make_pair(make_pair(chr, begin), end));
				chr = map_it->first.chr;
				begin = map_it->first.begin;
				end = map_it->first.end;
			}
		}
		++map_it;
	}
	begin2end.insert(make_pair(make_pair(chr, begin), end));

}

//int CountDepth(char *file[], const char *chr, int p, int baseQ, int mapQ, int n);
void OutputBreakpoint(ofstream &fout, multimap<Junction, OtherInfo> &junction2other, map<pair<string, int>, int> &pos2depth, map<ChrRange , unsigned long> &range2depth, map<Junction, pair<pair<ChrRange, ChrRange>, pair<ChrRange, ChrRange> > > &junction2range_pair, bool rescure_mode, int min_no_one_side_clipped_reads, int sum_min_no_both_clipped_reads, int min_abnormal_read_pair_no, double frequency, int min_distance, int max_microhomology, int min_seq_len, int max_seq_indel_no)
{
	multimap<Junction, OtherInfo>::iterator junction2other_it = junction2other.begin();
	map<Junction, pair<pair<ChrRange, ChrRange>, pair<ChrRange, ChrRange> > >::iterator junction2range_pair_it;
	map<ChrRange , unsigned long>::iterator range2depth_it;
	map<pair<string, int>, int>::iterator map_it;
	int updepth = 0, downdepth = 0;
	double rate1, rate2;

	//output head of the output file

	while (junction2other_it != junction2other.end())
	{
		map_it = pos2depth.find(make_pair(junction2other_it->first.up_chr, junction2other_it->first.up_pos));
		if (map_it == pos2depth.end())
		{
			cerr << "Error: There is something wrong in upstream position " << junction2other_it->first.up_chr << ":" << junction2other_it->first.up_pos << endl;
			updepth = 0;
		}
		else
			updepth = map_it->second + junction2other_it->second.down_seq_info.support_read_no;
		map_it = pos2depth.find(make_pair(junction2other_it->first.down_chr, junction2other_it->first.down_pos));
		if (map_it == pos2depth.end())
		{
			cerr << "Error: There is something wrong in downstream position " << junction2other_it->first.down_chr << ":" << junction2other_it->first.down_pos << endl;
			downdepth = 0;
		}
		else
			downdepth = map_it->second + junction2other_it->second.up_seq_info.support_read_no;
		int junction_reads_no = junction2other_it->second.up_seq_info.support_read_no + junction2other_it->second.down_seq_info.support_read_no;
		if (updepth == 0) rate1 = 0;
		else rate1 = (double)junction_reads_no / updepth;
		if (downdepth == 0) rate2 = 0;
		else rate2 = (double)junction_reads_no / downdepth;
		
		bool align_qual_pass = 0;
		//if (junction2other_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped + junction2other_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped >= 2 || junction2other_it->second.abnormal_read_pair_no >= 2) align_qual_pass = 1;
		if (junction2other_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped + junction2other_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped >= 2) align_qual_pass = 1;


		if ((junction2other_it->first.up_chr != junction2other_it->first.down_chr || (junction2other_it->first.up_chr == junction2other_it->first.down_chr && abs(junction2other_it->first.up_pos - junction2other_it->first.down_pos) >= min_distance)) && junction2other_it->second.microhomology_length <= max_microhomology && junction2other_it->second.abnormal_read_pair_no >= min_abnormal_read_pair_no && (rate1 >= frequency || rate2 >= frequency) && (junction2other_it->second.microhomology_length == -1 && rescure_mode && (junction2other_it->second.up_seq_info.support_read_no >= min_no_one_side_clipped_reads || junction2other_it->second.down_seq_info.support_read_no >= min_no_one_side_clipped_reads) || (junction2other_it->second.microhomology_length != -1 && junction2other_it->second.up_seq_info.support_read_no + junction2other_it->second.down_seq_info.support_read_no >= sum_min_no_both_clipped_reads)) && ((junction2other_it->second.abnormal_read_pair_no == 0 && junction2other_it->second.up_seq_info.seq.length() >= junction2other_it->second.up_seq_info.left_clipped_seq_length + junction2other_it->second.up_seq_info.right_clipped_seq_length + min_seq_len && junction2other_it->second.down_seq_info.seq.length() >= junction2other_it->second.down_seq_info.left_clipped_seq_length + junction2other_it->second.down_seq_info.right_clipped_seq_length + min_seq_len && junction2other_it->second.up_seq_info.cigar_vec.size() <= 2 * max_seq_indel_no + 1 && junction2other_it->second.down_seq_info.cigar_vec.size() <= 2 * max_seq_indel_no + 1 && CountLargestBaseFrequency(junction2other_it->second.up_seq_info.seq) < 0.8 && CountLargestBaseFrequency(junction2other_it->second.down_seq_info.seq) < 0.8) || junction2other_it->second.abnormal_read_pair_no > 0) && align_qual_pass)
		{
			unsigned int up_up_depth = 0, up_down_depth = 0, down_up_depth = 0, down_down_depth = 0;
			junction2range_pair_it = junction2range_pair.find(Junction(junction2other_it->first.up_chr, junction2other_it->first.up_pos, junction2other_it->first.up_strand, junction2other_it->first.down_chr, junction2other_it->first.down_pos, junction2other_it->first.down_strand));
			if (junction2range_pair_it != junction2range_pair.end())
			{
				range2depth_it = range2depth.find(junction2range_pair_it->second.first.first);	
				if (range2depth_it == range2depth.end())
				{
					cerr << "Error depth in the vicinity of junction " << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t'  << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << endl;
				}
				else
				{
					up_up_depth = range2depth_it->second / (range2depth_it->first.end - range2depth_it->first.begin + 1);
				}

				range2depth_it = range2depth.find(junction2range_pair_it->second.first.second);	
				if (range2depth_it == range2depth.end())
				{
					cerr << "Error depth in the vicinity of junction " << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t'  << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << endl;
				}
				else
				{
					up_down_depth = range2depth_it->second / (range2depth_it->first.end - range2depth_it->first.begin + 1);
				}

				range2depth_it = range2depth.find(junction2range_pair_it->second.second.first);	
				if (range2depth_it == range2depth.end())
				{
					cerr << "Error depth in the vicinity of junction " << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t'  << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << endl;
				}
				else
				{
					down_up_depth = range2depth_it->second / (range2depth_it->first.end - range2depth_it->first.begin + 1);
				}

				range2depth_it = range2depth.find(junction2range_pair_it->second.second.second);	
				if (range2depth_it == range2depth.end())
				{
					cerr << "Error depth in the vicinity of junction " << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t'  << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << endl;
				}
				else
				{
					down_down_depth = range2depth_it->second / (range2depth_it->first.end - range2depth_it->first.begin + 1);
				}
			}
			else
			{
				cerr << "Error depth in the vicinity of junction " << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t'  << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << endl;
			}
			string sv_type = GetSVType(junction2other_it->first.up_chr, junction2other_it->first.up_pos, junction2other_it->first.up_strand, junction2other_it->first.down_chr, junction2other_it->first.down_pos, junction2other_it->first.down_strand);
			fout << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t' << junction2other_it->second.up_seq_info.support_read_no << '\t' << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << '\t' << junction2other_it->second.down_seq_info.support_read_no << '\t' << junction2other_it->second.microhomology_length << '\t' << junction2other_it->second.abnormal_read_pair_no << '\t' << sv_type << '\t' << updepth << '\t' << downdepth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << rate1 << '\t' << rate2 << '\t';
			DisplayCigarVector<ostream> (fout, junction2other_it->second.up_seq_info.cigar_vec, junction2other_it->second.up_seq_info.left_clipped_seq_length, junction2other_it->second.up_seq_info.right_clipped_seq_length);
			fout << '\t';
			DisplayCigarVector<ostream> (fout, junction2other_it->second.down_seq_info.cigar_vec, junction2other_it->second.down_seq_info.left_clipped_seq_length, junction2other_it->second.down_seq_info.right_clipped_seq_length);
			fout << '\t' << junction2other_it->second.up_seq_info.seq << '\t' << junction2other_it->second.down_seq_info.seq << endl;

		}
		//test
		else {
			cout << "updepth:" << updepth << '\t' << "downdepth" << downdepth << endl;
			string sv_type = GetSVType(junction2other_it->first.up_chr, junction2other_it->first.up_pos, junction2other_it->first.up_strand, junction2other_it->first.down_chr, junction2other_it->first.down_pos, junction2other_it->first.down_strand);
			cout << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t' << junction2other_it->second.up_seq_info.support_read_no << '\t' << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << '\t' << junction2other_it->second.down_seq_info.support_read_no << '\t' << junction2other_it->second.microhomology_length << '\t' << junction2other_it->second.abnormal_read_pair_no << '\t' << sv_type << endl;
		
		}
		++junction2other_it;
	}
	
}


void FindDiscordantReadPairs(samfile_t *samfin, bam_index_t *idx, multimap<Junction, OtherInfo> &junction2other, int min_mapQ, int mean_insert_size, int deviation, int times)
{

//	samfile_t *samfin;
	bam1_t *b;
//	char in_mode[5], *fn_list = 0;
//	in_mode[0] = 'r';
//	int is_bamin = 0;
//	if (file.rfind(".bam") == file.size() - 4)
//	{
//		//if in_mode[1] == 'b', it will read a bam file
//		in_mode[1] = 'b';
//		is_bamin = 1;
//	}
//
//    if ((samfin = samopen(file.c_str(), in_mode, fn_list)) == 0)
//	{
//		cerr << "[main_samview] fail to open file for reading." << endl;
//   		exit(1);
//	}
//	if (samfin->header == 0) 
//	{
//		cerr << "[main_samview] fail to read the header." << endl;
//		exit(1);
//	}
//
	int ret = 0;
//
//	bam_index_t *idx = 0;
//	if (is_bamin) idx = bam_index_load(file.c_str()); // load BAM index
//	if (idx == 0)
//	{ // index is unavailable
//		cerr << "[main_samview] random alignment retrieval only works for indexed BAM files.\n" << endl;
//		cerr << file << endl;
//		exit(1);
//	}

	g_min_mapQ = min_mapQ;
	map<string, int> seq_name2tid;
	StoreSeqName2Tid(samfin->header, seq_name2tid);
	cerr << "'StoreSeqName2Tid' finished" << endl;

	int min_insert_size = mean_insert_size - deviation * times, max_insert_size = mean_insert_size + deviation * times;	
	if (min_insert_size < 0)
		min_insert_size = 0;

	//beg is 0-base, end is 1-base
	int tid = -1, beg = - 1, end = -1, mtid = -1, mbeg = -1, mend = -1;
	multimap<Junction, OtherInfo>::iterator junction2other_it = junction2other.begin();
	while (junction2other_it != junction2other.end())
	{
		tid = BamGetTid(seq_name2tid, junction2other_it->first.up_chr);
		unsigned int chr_length = samfin->header->target_len[tid];
		if (tid == -1) { ++ junction2other_it; continue; }
		if (junction2other_it->first.up_strand == '+')
		{
			end = junction2other_it->first.up_pos;
			beg = end - max_insert_size;
		}
		else if (junction2other_it->first.up_strand == '-')
		{
			beg = junction2other_it->first.up_pos - 1 - kCrossLength;
			end = junction2other_it->first.up_pos - 1 + max_insert_size;
		}
		else 
		{
			cerr << "Error strand of junction " << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << endl;
			continue;
		}
		if (beg <= 0) beg = 1; 
		if (end > chr_length) end = chr_length;


		bam_iter_t iter = bam_iter_query(idx, tid, beg, end);
		bamFile fp = samfin->x.bam;
		b = bam_init1();
		int read_pair_no = 0;
		while ((ret = bam_iter_read(fp, iter, b)) >= 0)
		{
			if (__g_skip_aln(samfin->header, b)) continue;
			// Add this to ignore hard clipped reads 2016-07-07
			if (IsHardClip(b)) continue;
			if ((b->core.flag&BAM_FDUP) || (b->core.flag&BAM_FUNMAP) || (b->core.flag&BAM_FMUNMAP)|| IsConcordant(samfin->header, b, mean_insert_size, deviation, times)) continue;
			mtid = BamGetTid(seq_name2tid, junction2other_it->first.down_chr);
			if (mtid != -1 && mtid == b->core.mtid) //&& b->core.pos + b->core.l_qseq <= junction2other_it->first.up_pos + kCrossLength && b->core.mpos + 1 >= junction2other_it->first.down_pos - kCrossLength)
			{
				if (junction2other_it->first.up_strand == '+' && junction2other_it->first.down_strand == '+' && b->core.pos + b->core.l_qseq <= junction2other_it->first.up_pos + kCrossLength && b->core.mpos + 1 >= junction2other_it->first.down_pos - kCrossLength)
				{
					if (!(b->core.flag&BAM_FREVERSE) && (b->core.flag&BAM_FMREVERSE))
					{
						//max insert size is shorter than tandem_duplication_length + 2 * insert_size
						if (tid == mtid && junction2other_it->first.up_pos > junction2other_it->first.down_pos && junction2other_it->first.up_pos - junction2other_it->first.down_pos + 1 + 2 * b->core.l_qseq <= max_insert_size)
						{
							int insert_size = junction2other_it->first.up_pos - b->core.pos + b->core.mpos + b->core.l_qseq - junction2other_it->first.down_pos + 1;
							int mark = 0;
							while (insert_size <= max_insert_size)
							{
								if (insert_size >= min_insert_size) { mark = 1;	break; }
								else { insert_size += junction2other_it->first.up_pos - junction2other_it->first.down_pos + 1;}
							}
							if (mark) { ++ read_pair_no; }
						}
						else 
						{
							int insert_size = junction2other_it->first.up_pos - b->core.pos + b->core.mpos + b->core.l_qseq - junction2other_it->first.down_pos + 1;
							if (min_insert_size <= insert_size && insert_size <= max_insert_size) { ++read_pair_no; }
						}
					}
				}
				// left sequence is inverted sequence, strand of discordant read pair is - -
				else if (junction2other_it->first.up_strand == '-' && junction2other_it->first.down_strand == '+' &&  (b->core.flag&BAM_FREVERSE) && (b->core.flag&BAM_FMREVERSE) && b->core.mpos + 1 >= junction2other_it->first.down_pos - kCrossLength)
				{

					int insert_size = b->core.pos + 1 - junction2other_it->first.up_pos + 1 + b->core.mpos + b->core.l_qseq - junction2other_it->first.down_pos + 1;
					if (min_insert_size <= insert_size && insert_size <= max_insert_size) { ++read_pair_no; }
				}
				// right sequence is inverted sequence, strand of discordant read pair is + + 
				else if (junction2other_it->first.up_strand == '+' && junction2other_it->first.down_strand == '-' && !(b->core.flag&BAM_FREVERSE) && !(b->core.flag&BAM_FMREVERSE) && b->core.pos + b->core.l_qseq <= junction2other_it->first.up_pos + kCrossLength && b->core.mpos + b->core.l_qseq <= junction2other_it->first.down_pos + kCrossLength)
				{
					int insert_size = junction2other_it->first.up_pos - b->core.pos + junction2other_it->first.down_pos - (b->core.mpos + b->core.l_qseq) + 1;
					if (min_insert_size <= insert_size && insert_size <= max_insert_size) { ++read_pair_no; }
				}
			}
			//else {cerr << bam1_qname(b) << endl;}
			
		}
		junction2other_it->second.abnormal_read_pair_no = read_pair_no;
		++ junction2other_it;
	}

}


int FindDiscordantReadPairs(samfile_t *samfin, bam_index_t *idx, Junction &junction, int min_mapQ, int mean_insert_size, int deviation, int times)
{
//	samfile_t *samfin;
	bam1_t *b;
//	char in_mode[5], *fn_list = 0;
//	in_mode[0] = 'r';
//	int is_bamin = 0;
//	if (file.rfind(".bam") == file.size() - 4)
//	{
//		//if in_mode[1] == 'b', it will read a bam file
//		in_mode[1] = 'b';
//		is_bamin = 1;
//	}
//
//    if ((samfin = samopen(file.c_str(), in_mode, fn_list)) == 0)
//	{
//		cerr << "[main_samview] fail to open file for reading." << endl;
//   		exit(1);
//	}
//	if (samfin->header == 0) 
//	{
//		cerr << "[main_samview] fail to read the header." << endl;
//		exit(1);
//	}
//
	int ret = 0;
//
//	bam_index_t *idx = 0;
//	if (is_bamin) idx = bam_index_load(file.c_str()); // load BAM index
//	if (idx == 0)
//	{ // index is unavailable
//		cerr << "[main_samview] random alignment retrieval only works for indexed BAM files.\n" << endl;
//		cerr << idx << '\t' <<  is_bamin << '\t' << file << endl;
//		exit(1);
//	}

	g_min_mapQ = min_mapQ;
	map<string, int> seq_name2tid;
	StoreSeqName2Tid(samfin->header, seq_name2tid);
	int min_insert_size = mean_insert_size - deviation * times, max_insert_size = mean_insert_size + deviation * times;	
	if (min_insert_size < 0)
		min_insert_size = 0;
	//beg is 0-base, end is 1-base
	int tid = -1, beg = - 1, end = -1, mtid = -1, mbeg = -1, mend = -1;

	tid = BamGetTid(seq_name2tid, junction.up_chr);
	unsigned int chr_length = samfin->header->target_len[tid];
	if (tid == -1) 
	{
		cerr << "Cannot find " << junction.up_chr << ", please check whether you input a wrong file" << endl;
	}
	if (junction.up_strand == '+')
	{
		end = junction.up_pos;
		beg = end - max_insert_size;
	}
	else if (junction.up_strand == '-')
	{
		beg = junction.up_pos - 1 - kCrossLength;
		end = junction.up_pos - 1 + max_insert_size;
	}
	else 
	{
		cerr << "Error strand of junction " << junction.up_chr << '\t' << junction.up_pos << endl;
	}
	if (beg <= 0) beg = 1; 
	if (end > chr_length) end = chr_length;


	bam_iter_t iter = bam_iter_query(idx, tid, beg, end);
	bamFile fp = samfin->x.bam;
	b = bam_init1();
	int read_pair_no = 0;
	while ((ret = bam_iter_read(fp, iter, b)) >= 0)
	{
		if (__g_skip_aln(samfin->header, b)) continue;
			// Add this to ignore hard clipped reads 2016-07-07
		if (IsHardClip(b)) continue;
		if ((b->core.flag&BAM_FDUP) || (b->core.flag&BAM_FUNMAP) || (b->core.flag&BAM_FMUNMAP)|| IsConcordant(samfin->header, b, mean_insert_size, deviation, times)) continue;
		mtid = BamGetTid(seq_name2tid, junction.down_chr);
		if (mtid != -1 && mtid == b->core.mtid) //&& b->core.pos + b->core.l_qseq <= junction.up_pos + kCrossLength && b->core.mpos + 1 >= junction.down_pos - kCrossLength)
		{
			if (junction.up_strand == '+' && junction.down_strand == '+' && b->core.pos + b->core.l_qseq <= junction.up_pos + kCrossLength && b->core.mpos + 1 >= junction.down_pos - kCrossLength)
			{
				if (!(b->core.flag&BAM_FREVERSE) && (b->core.flag&BAM_FMREVERSE))
				{
					//max insert size is shorter than tandem_duplication_length + 2 * insert_size
					if (tid == mtid && junction.up_pos > junction.down_pos && junction.up_pos - junction.down_pos + 1 + 2 * b->core.l_qseq <= max_insert_size)
					{
						int insert_size = junction.up_pos - b->core.pos + b->core.mpos + b->core.l_qseq - junction.down_pos + 1;
						int mark = 0;
						while (insert_size <= max_insert_size)
						{
							if (insert_size >= min_insert_size) { mark = 1;	break; }
							else { insert_size += junction.up_pos - junction.down_pos + 1;}
						}
						if (mark) { ++ read_pair_no; }
					}
					else 
					{
						int insert_size = junction.up_pos - b->core.pos + b->core.mpos + b->core.l_qseq - junction.down_pos + 1;
						if (min_insert_size <= insert_size && insert_size <= max_insert_size) { ++read_pair_no; }
					}
				}
			}
			// left sequence is inverted sequence, strand of discordant read pair is - -
			else if (junction.up_strand == '-' && junction.down_strand == '+' &&  (b->core.flag&BAM_FREVERSE) && (b->core.flag&BAM_FMREVERSE) && b->core.mpos + 1 >= junction.down_pos - kCrossLength)
			{

				int insert_size = b->core.pos + 1 - junction.up_pos + 1 + b->core.mpos + b->core.l_qseq - junction.down_pos + 1;
				if (min_insert_size <= insert_size && insert_size <= max_insert_size) { ++read_pair_no; }
			}
			// right sequence is inverted sequence, strand of discordant read pair is + + 
			else if (junction.up_strand == '+' && junction.down_strand == '-' && !(b->core.flag&BAM_FREVERSE) && !(b->core.flag&BAM_FMREVERSE) && b->core.pos + b->core.l_qseq <= junction.up_pos + kCrossLength && b->core.mpos + b->core.l_qseq <= junction.down_pos + kCrossLength)
			{
				int insert_size = junction.up_pos - b->core.pos + junction.down_pos - (b->core.mpos + b->core.l_qseq) + 1;
				if (min_insert_size <= insert_size && insert_size <= max_insert_size) { ++read_pair_no; }
			}
		}
		//else {cerr << bam1_qname(b) << endl;}
		
	}

	return read_pair_no;
}




void OutputOneendUnmapBreakpoint(ofstream &fout, ofstream &fout1, multimap<pair<string, int>, ClipReads> &aligned2clipped, map<pair<string, int>, int> &pos2depth)
{
	multimap<pair<string, int>, ClipReads>::iterator mmap_it1 = aligned2clipped.begin();
	map<pair<string, int>, int>::iterator map_it;
	double rate;
	while (mmap_it1 != aligned2clipped.end())
	{
		if (mmap_it1->second.type == 'n')
		{
			map_it = pos2depth.find(mmap_it1->first);
			if (map_it == pos2depth.end())
			{
				cerr << "Error: There is something wrong in position " << mmap_it1->first.first << ":" << mmap_it1->first.second << endl;
				++mmap_it1;
				continue;
			}
			int depth = map_it->second;;
			if (depth == 0) rate = 0;
			else rate = (double)mmap_it1->second.aligned_seq_info.support_read_no / depth;
			/*
			if (mmap_it1->second.clipped_side == '5')
			{
				fout << "*\t-1\t*\t0\t" << mmap_it1->first.first << '\t' << mmap_it1->first.second << "\t+\t" << mmap_it1->second.aligned_seq_info.support_read_no << "\t-1\t*\t-1\t" << depth << "\t0\t" << rate << '\t' << mmap_it1->second.aligned_seq_info.cigar << '\t' << mmap_it1->second.clipped_seq << '\t' << mmap_it1->second.aligned_seq_info.seq << endl;
			}
			else if (mmap_it1->second.clipped_side == '3')
			{
				fout << mmap_it1->first.first << '\t' << mmap_it1->first.second << "\t+\t" << mmap_it1->second.aligned_seq_info.support_read_no << "\t*\t-1\t*\t0\t-1\t" << depth << "\t-1\t" << rate << "\t0\t" << mmap_it1->second.aligned_seq_info.cigar << "\t*\t" << mmap_it1->second.aligned_seq_info.seq << '\t' << mmap_it1->second.clipped_seq << endl;
			}
			*/

			fout1 << "@" << mmap_it1->second.clipped_seq << '\n' << mmap_it1->second.clipped_seq << "\n+\n" << mmap_it1->second.clipped_qual << endl;
		}
		++mmap_it1;
	}
	fout.close();
	fout.clear();
}



void ReadBreakpoint(string file, multimap<Junction, OtherInfo> &junction2other)
{
	ifstream fin(file.c_str());
	if (!fin) {cerr << "Cannot open file " << file << endl; return; }
	
	Junction junction;
	SeqInfo up_seq_info;
	SeqInfo down_seq_info;
	string up_chr, down_chr, sv_type, up_cigar, down_cigar, up_seq, down_seq, temp;
	int up_pos, up_reads_no, down_pos, down_reads_no, microhomology_length, abnormal_read_pair_no, up_depth, down_depth, up_up_depth, up_down_depth, down_up_depth, down_down_depth;
	char up_strand, down_strand;
	double up_clip_rate, down_clip_rate;

	while (fin >> up_chr)
	{
		if (up_chr[0] == '@') 
		{
			getline(fin,temp);
			continue;
		}
		fin >> up_pos >> up_strand >> up_reads_no >> down_chr >> down_pos >> down_strand >> down_reads_no >> microhomology_length >> abnormal_read_pair_no >> sv_type >> up_depth >> down_depth >> up_up_depth >> up_down_depth >> down_up_depth >> down_down_depth >> up_clip_rate >> down_clip_rate >> up_cigar >> down_cigar >> up_seq >> down_seq;
		getline(fin, temp);
		junction.set_value(up_chr, up_pos, up_strand, down_chr, down_pos, down_strand);
		vector<pair<int, char> > up_cigar_vec = ChangeCigarType(up_cigar), down_cigar_vec = ChangeCigarType(down_cigar);
		up_seq_info.set_value(up_seq, up_cigar_vec, 0, 0, up_reads_no, 0);
		down_seq_info.set_value(down_seq, down_cigar_vec, 0, 0, down_reads_no, 0);
		OtherInfo other_info(up_seq_info, down_seq_info, microhomology_length, abnormal_read_pair_no);
		junction2other.insert(make_pair(junction, other_info));


	}
}

void MergeJunction(multimap<Junction, OtherInfo> &junction2other, int search_length)
{
	
	multimap<Junction, OtherInfo>::iterator junction2other_it = junction2other.begin(), junction2other_backup_it;

	//tcb
//	while (junction2other_it != junction2other.end())
//	{
//		Junction junction_tmp(junction2other_it->first);
//		cerr << junction_tmp << '\t' << junction2other_it->second.up_seq_info.support_read_no << '\t' << junction2other_it->second.down_seq_info.support_read_no << '\t' << junction2other_it->second.microhomology_length << '\t';
//		DisplayCigarVector<ostream>(cerr, junction2other_it->second.up_seq_info.cigar_vec, junction2other_it->second.up_seq_info.left_clipped_seq_length, junction2other_it->second.down_seq_info.right_clipped_seq_length);
//		cerr << '\t';
//		DisplayCigarVector<ostream>(cerr, junction2other_it->second.down_seq_info.cigar_vec, junction2other_it->second.down_seq_info.left_clipped_seq_length, junction2other_it->second.down_seq_info.right_clipped_seq_length);
//		cerr << '\t' << junction2other_it->second.up_seq_info.seq << '\t' << junction2other_it->second.down_seq_info.seq << endl;
//		++junction2other_it;
//	}
//	junction2other_it = junction2other.begin();
	//tce

	while (junction2other_it != junction2other.end())
	{
		//May be there is a bug here
		if (junction2other_it->second.up_seq_info.right_clipped_seq_length > 0 || junction2other_it->second.up_seq_info.left_clipped_seq_length > 0)
		{
			++junction2other_it;
			continue;
		}
		junction2other_backup_it = junction2other_it;
		junction2other_backup_it++;
		int mark = 0;
		while (junction2other_backup_it != junction2other.end() && junction2other_it->first.up_chr == junction2other_backup_it->first.up_chr && junction2other_it->first.down_chr == junction2other_backup_it->first.down_chr && junction2other_it->first.up_strand == junction2other_backup_it->first.up_strand && junction2other_it->first.down_strand == junction2other_backup_it->first.down_strand && junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos <= search_length)
		{
			if (abs(junction2other_backup_it->first.down_pos - junction2other_it->first.down_pos) <= search_length && junction2other_backup_it->second.down_seq_info.left_clipped_seq_length == 0)
			{
				string up_seq1, down_seq1, up_seq2, down_seq2, temp;
				if (junction2other_it->second.up_seq_info.cigar_vec.size() == 1 && junction2other_backup_it->second.up_seq_info.cigar_vec.size() == 1)
				{
					int microhomology_length = junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos;
					if (junction2other_it->first.up_strand == '+' && junction2other_backup_it->second.up_seq_info.seq.length() < microhomology_length + 5 || junction2other_it->first.up_strand == '-' && junction2other_it->second.up_seq_info.seq.length() < microhomology_length + 5)
					{
						++junction2other_backup_it;
						continue;
					}
					if (junction2other_it->first.up_strand == '+')
					{
						up_seq1 = junction2other_it->second.up_seq_info.seq;
						down_seq1 = junction2other_it->second.down_seq_info.seq;
						up_seq2 = junction2other_backup_it->second.up_seq_info.seq.substr(0, junction2other_backup_it->second.up_seq_info.seq.length() - microhomology_length);
						down_seq2 = junction2other_backup_it->second.up_seq_info.seq.substr(junction2other_backup_it->second.up_seq_info.seq.length() - microhomology_length);
						down_seq2.append(junction2other_backup_it->second.down_seq_info.seq);
					}
					else
					{
						up_seq1 = junction2other_it->second.up_seq_info.seq.substr(0, junction2other_it->second.up_seq_info.seq.length() - microhomology_length);
						down_seq1 = junction2other_it->second.up_seq_info.seq.substr(junction2other_it->second.up_seq_info.seq.length() - microhomology_length);
						down_seq1.append(junction2other_it->second.down_seq_info.seq);
						up_seq2 = junction2other_backup_it->second.up_seq_info.seq;
						down_seq2 = junction2other_backup_it->second.down_seq_info.seq;
					}
				}
				else if (junction2other_it->second.down_seq_info.cigar_vec.size() == 1 && junction2other_backup_it->second.down_seq_info.cigar_vec.size() == 1)
				{
					int microhomology_length = abs(junction2other_backup_it->first.down_pos - junction2other_it->first.down_pos);
					if (junction2other_it->first.up_strand == '+' && junction2other_it->second.down_seq_info.seq.length() < microhomology_length + 5 || junction2other_it->first.up_strand == '-' && junction2other_backup_it->second.down_seq_info.seq.length() < microhomology_length + 5)
					{
						++junction2other_backup_it;
						continue;
					}
					
					if (junction2other_it->first.up_strand == '+')
					{
						down_seq1 = junction2other_it->second.down_seq_info.seq.substr(microhomology_length);
						down_seq2 = junction2other_backup_it->second.down_seq_info.seq;
						up_seq1 = junction2other_it->second.up_seq_info.seq;
						up_seq1.append(junction2other_it->second.down_seq_info.seq, 0, microhomology_length);
						up_seq2 = junction2other_backup_it->second.up_seq_info.seq;
					}
					else
					{
						down_seq1 = junction2other_it->second.down_seq_info.seq;
						down_seq2 = junction2other_backup_it->second.down_seq_info.seq.substr(microhomology_length);
						up_seq1 = junction2other_it->second.up_seq_info.seq;
						up_seq2 = junction2other_backup_it->second.up_seq_info.seq;
						up_seq2.append(junction2other_backup_it->second.down_seq_info.seq, 0, microhomology_length);
					}
				}
				if (CompareStringEndFirst(up_seq1, up_seq2) >= 0.85 && CompareStringBeginFirst(down_seq1, down_seq2) >= 0.85)
				{
					junction2other_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped = junction2other_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped > junction2other_backup_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped ? junction2other_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped : junction2other_backup_it->second.up_seq_info.is_clipped_seq_and_uniq_mapped;
					junction2other_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped = junction2other_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped > junction2other_backup_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped ? junction2other_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped : junction2other_backup_it->second.down_seq_info.is_clipped_seq_and_uniq_mapped;

					if (junction2other_it->second.microhomology_length == -1 && junction2other_backup_it->second.microhomology_length == -1)
					{
						junction2other_it->second.up_seq_info.support_read_no += junction2other_backup_it->second.up_seq_info.support_read_no;
						junction2other_it->second.down_seq_info.support_read_no += junction2other_backup_it->second.down_seq_info.support_read_no;



						if (junction2other_it->second.up_seq_info.support_read_no != 0 && junction2other_backup_it->second.down_seq_info.support_read_no != 0 || junction2other_it->second.down_seq_info.support_read_no != 0 && junction2other_backup_it->second.up_seq_info.support_read_no != 0)
						{
							junction2other_it->second.microhomology_length = junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos;
						}
						junction2other.erase(junction2other_backup_it++);
					}
					else if (junction2other_it->second.microhomology_length != -1 && junction2other_backup_it->second.microhomology_length == -1)
					{
						junction2other_it->second.up_seq_info.support_read_no += junction2other_backup_it->second.up_seq_info.support_read_no;
						junction2other_it->second.down_seq_info.support_read_no += junction2other_backup_it->second.down_seq_info.support_read_no;
						junction2other.erase(junction2other_backup_it++);
					}
					else if (junction2other_it->second.microhomology_length == -1 && junction2other_backup_it->second.microhomology_length != -1)
					{
						junction2other_backup_it->second.up_seq_info.support_read_no += junction2other_it->second.up_seq_info.support_read_no;
						junction2other_backup_it->second.down_seq_info.support_read_no += junction2other_it->second.down_seq_info.support_read_no;
						mark = 1;	
					}
					else
					{
						if (junction2other_it->second.up_seq_info.support_read_no > junction2other_backup_it->second.up_seq_info.support_read_no || junction2other_it->second.down_seq_info.support_read_no == junction2other_backup_it->second.down_seq_info.support_read_no)
						{
							junction2other_it->second.up_seq_info.support_read_no += junction2other_backup_it->second.up_seq_info.support_read_no;
							junction2other.erase(junction2other_backup_it++);
						}
						else if (junction2other_it->second.up_seq_info.support_read_no == junction2other_backup_it->second.up_seq_info.support_read_no || junction2other_it->second.down_seq_info.support_read_no > junction2other_backup_it->second.down_seq_info.support_read_no)
						{
							junction2other_it->second.down_seq_info.support_read_no += junction2other_backup_it->second.down_seq_info.support_read_no;
							junction2other.erase(junction2other_backup_it++);
						}
						else if (junction2other_backup_it->second.up_seq_info.support_read_no > junction2other_it->second.up_seq_info.support_read_no && junction2other_it->second.down_seq_info.support_read_no == junction2other_backup_it->second.down_seq_info.support_read_no)
						{
							junction2other_backup_it->second.up_seq_info.support_read_no += junction2other_it->second.up_seq_info.support_read_no;
							mark = 1;
						}
						else if (junction2other_backup_it->second.down_seq_info.support_read_no > junction2other_it->second.down_seq_info.support_read_no && junction2other_backup_it->second.up_seq_info.support_read_no == junction2other_it->second.up_seq_info.support_read_no)
						{
							junction2other_backup_it->second.down_seq_info.support_read_no += junction2other_it->second.down_seq_info.support_read_no;
							mark = 1;
						}
						else 
							junction2other_backup_it++;
					}
					if (mark == 1)
						break;
				}
				else
				{
					junction2other_backup_it++;
				}
			}
			else
				junction2other_backup_it++;
		}
		if (mark == 1)
			junction2other.erase(junction2other_it++);
		else
			junction2other_it++;
	}
}


double CountLargestBaseFrequency(string &seq)
{
	int length = seq.length();
	int Anumber = 0, Tnumber = 0, Cnumber = 0, Gnumber = 0, Nnumber = 0, largest = 0;
	for (int i = 0; i < length; i++)
	{
		if (seq[i] == 'A' || seq[i] == 'a')
			Anumber ++;
		else if (seq[i] == 'T' || seq[i] == 't')
			Tnumber ++;
		else if (seq[i] == 'C' || seq[i] == 'c')
			Cnumber ++;
		else if (seq[i] == 'G' || seq[i] == 'g')
			Gnumber ++;
		else 
			Nnumber ++;
	}
	if (largest < Anumber) largest = Anumber;
	if (largest < Tnumber) largest = Tnumber;
	if (largest < Cnumber) largest = Cnumber;
	if (largest < Gnumber) largest = Gnumber;
	if (largest < Nnumber) largest = Nnumber;
	
	double rate = largest/(double)length;
	return rate;

}

/*
void MergeJunction(multimap<Junction, OtherInfo> &junction2other, int search_length)
{
	multimap<Junction, OtherInfo>::iterator junction2other_it = junction2other.begin(), junction2other_backup_it;
	while (junction2other_it != junction2other.end())
	{
		junction2other_backup_it = junction2other_it;
		junction2other_backup_it++;
		while (junction2other_backup_it != junction2other.end() && junction2other_it->first.up_chr == junction2other_backup_it->first.up_chr && junction2other_it->first.down_chr == junction2other_backup_it->first.down_chr && junction2other_it->first.up_strand == junction2other_backup_it->first.up_strand && junction2other_it->first.down_strand == junction2other_backup_it->first.down_strand && junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos <= search_length)
		{
			if (junction2other_it->first.up_strand == junction2other_it->first.down_strand && junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos == junction2other_backup_it->first.down_pos - junction2other_it->first.down_pos || junction2other_it->first.up_strand != junction2other_it->first.down_strand && junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos == junction2other_it->first.down_pos - junction2other_backup_it->first.down_pos)
			{
				junction2other_it->second.up_seq_info.support_read_no += junction2other_backup_it->second.up_seq_info.support_read_no;
				junction2other_it->second.down_seq_info.support_read_no += junction2other_backup_it->second.down_seq_info.support_read_no;
				junction2other_it->second.microhomology_length = junction2other_backup_it->first.up_pos - junction2other_it->first.up_pos;
				junction2other.erase(junction2other_backup_it);
				break;
			}
			else if (abs(junction2other_backup_it->first.down_pos - junction2other_it->first.down_pos) <= search_length)
			{
				int new_up_pos = 0, new_down_pos = 0;
				multimap<Junction, OtherInfo>::iterator junction2other_not_equal_it, junction2other_equal_it;
				pair<int, int> len_pair;
				int size1 = junction2other_it->second.up_seq_info.cigar_vec.size();
				int size2 = junction2other_it->second.down_seq_info.cigar_vec.size();
				int size3 = junction2other_backup_it->second.up_seq_info.cigar_vec.size(); 
				int size4 = junction2other_backup_it->second.down_seq_info.cigar_vec.size();
				if (size1 == 1 && size2 == 1 && (size3 == 1 && size4 != 1 || size3 != 1 && size4 == 1) ||
					size3 == 1 && size4 == 1 && (size1 == 1 && size2 != 1 || size1 != 1 && size2 == 1))

				{
					if (size3 != 1 || size4 != 1) 
					{
						junction2other_not_equal_it = junction2other_backup_it;
						junction2other_equal_it = junction2other_it;

					}
					else
					{ 
						junction2other_not_equal_it = junction2other_it;
						junction2other_equal_it = junction2other_backup_it;
					}

					if (junction2other_not_equal_it->second.down_seq_info.cigar_vec.size() != 1)
					{
						len_pair = DownSeqExceedLength(junction2other_not_equal_it->second.down_seq_info.cigar_vec);
					}
					else
					{
						len_pair = UpSeqExceedLength(junction2other_not_equal_it->second.down_seq_info.cigar_vec);
						len_pair.first = -len_pair.first;
						len_pair.second = -len_pair.second;
					}

					if (junction2other_not_equal_it->first.up_strand == junction2other_not_equal_it->first.down_strand)
					{
						new_up_pos = junction2other_not_equal_it->first.up_pos + len_pair.first;
						new_down_pos = junction2other_not_equal_it->first.down_pos + len_pair.second;
					}
					else if (junction2other_not_equal_it->first.up_strand == '-' && junction2other_not_equal_it->first.down_strand == '+')
					{
						new_up_pos = junction2other_not_equal_it->first.up_pos - len_pair.first;
						new_down_pos = junction2other_not_equal_it->first.down_pos + len_pair.second;
					}
					else if (junction2other_not_equal_it->first.up_strand == '+' && junction2other_not_equal_it->first.down_strand == '-')
					{
						new_up_pos = junction2other_not_equal_it->first.up_pos + len_pair.first;
						new_down_pos = junction2other_not_equal_it->first.down_pos - len_pair.second;
					}

					if (junction2other_equal_it->first.up_strand == junction2other_equal_it->first.down_strand && 
						new_up_pos - junction2other_equal_it->first.up_pos == new_down_pos - junction2other_equal_it->first.down_pos ||
						junction2other_equal_it->first.up_strand != junction2other_equal_it->first.down_strand &&
						new_up_pos + new_down_pos == junction2other_equal_it->first.up_pos + junction2other_equal_it->first.down_pos)
					{
						if (new_up_pos >= junction2other_equal_it->first.up_pos)
						{
							junction2other_equal_it->second.up_seq_info.support_read_no += junction2other_not_equal_it->second.up_seq_info.support_read_no;
							junction2other_equal_it->second.down_seq_info.support_read_no += junction2other_not_equal_it->second.down_seq_info.support_read_no;
							junction2other_equal_it->second.microhomology_length = new_up_pos - junction2other_equal_it->first.up_pos;


							if (junction2other_not_equal_it == junction2other_it)
							{
								if (junction2other_it == junction2other.begin()) 
								{
									junction2other.erase(junction2other_not_equal_it);
									junction2other_it = junction2other.begin();
								}
								else
								{
									junction2other_it--;
									junction2other.erase(junction2other_not_equal_it);
								}
							}
							else
								junction2other.erase(junction2other_not_equal_it);
						}
						else
						{
           					Junction new_junciton(junction2other_equal_it->first.up_chr, new_up_pos, junction2other_equal_it->first.up_strand, junction2other_equal_it->first.down_chr, new_down_pos, junction2other_equal_it->first.down_strand);
							int microhomology_length = junction2other_equal_it->first.up_pos - new_up_pos;
							string new_up_seq = junction2other_equal_it->second.up_seq_info.seq;
							string new_down_seq = junction2other_equal_it->second.down_seq_info.seq;
							string temp =  new_up_seq.substr(new_up_seq.length() - 2);
							temp.append(new_down_seq);
							new_down_seq = temp;
							new_up_seq = new_up_seq.substr(0, new_up_seq.length() - 2);
							vector<pair<int, char> > up_cigar_vec = junction2other_equal_it->second.up_seq_info.cigar_vec;
							vector<pair<int, char> > down_cigar_vec = junction2other_equal_it->second.up_seq_info.cigar_vec;
							up_cigar_vec[0].first -= microhomology_length;
							down_cigar_vec[0].first += microhomology_length;
							SeqInfo up_seq_info(new_up_seq, up_cigar_vec, junction2other_equal_it->second.up_seq_info.left_clipped_seq_length, junction2other_equal_it->second.up_seq_info.right_clipped_seq_length, junction2other_equal_it->second.up_seq_info.support_read_no + junction2other_not_equal_it->second.up_seq_info.support_read_no);
							SeqInfo down_seq_info(new_down_seq, down_cigar_vec, junction2other_equal_it->second.down_seq_info.left_clipped_seq_length, junction2other_equal_it->second.down_seq_info.right_clipped_seq_length, junction2other_equal_it->second.down_seq_info.support_read_no + junction2other_not_equal_it->second.down_seq_info.support_read_no);
							OtherInfo new_other_info(up_seq_info, down_seq_info, microhomology_length, 0);


							if (junction2other_it == junction2other.begin()) 
							{
								junction2other.erase(junction2other_it);
								junction2other.erase(junction2other_backup_it);
								junction2other_it = junction2other.begin();
							}
							else
							{
								junction2other_it--;
								junction2other.erase(junction2other_not_equal_it);
								junction2other.erase(junction2other_equal_it);
							}
							junction2other.insert(make_pair(new_junciton, new_other_info));
						}
						break;
					}
					else
					{
					}
				}

			}
			++junction2other_backup_it;
		}
		++ junction2other_it;
	}
}

*/

/*
bool GetSmallPosAndCigar(int &pos, string &cigar, int microhomology_len)
{
	int l, i;
	char m = cigar[cigar.length() - 1];
	int len = cigar.length();
	string temp, str_new_len;
	for (i = 2; i <= len; ++i)
	{
		if (!isdigit(cigar[len - i]))
		{

			break;
		}
	}
	temp = cigar.substr(len - i + 1);
	l = atoi(temp.c_str());
	if (l > minus_len)
	{
		l -= add_len;
		itoa(l, str_new_len); 
		cigar = cigar.substr(0, len - i + 1);
		cigar.append(str_new_len);
		cigar.append(1, m);
		pos -= microhomology_len;
		return 0;
	}
	else
	{
		if (cigar[len - i] == 'D' || cigar[len - i] == 'N')
		{
			for (i = i - 1; i <= len; ++i)
			{
				if (!isdigit(cigar[len - i]))
					break;
			}
		}
	}
}
*/
//double CompareStringEndFirst(string str1, string str2);
//double CompareStringBeginFirst(string str1, string str2);



void GetJunction(AlignReadsInfo &align_reads_info, char orientation, AlignInfo &clipped_align_info, multimap<Junction, OtherInfo> &junction2other, multimap<pair<string, int>, ClipReads> &aligned2clipped)
{
//if (read_id2align_info_it_pair.first != read_id2align_info_it_pair.second)

	string chr = align_reads_info.chr;
	int pos = align_reads_info.pos;
	vector<pair<int, char> > cigar_vec = align_reads_info.reads_info.cigar_vec;
	string aligned_seq = align_reads_info.reads_info.seq_left;
	string clipped_seq = align_reads_info.reads_info.seq_right;
	string clipped_qual = align_reads_info.reads_info.qual_right;
	int support_count = align_reads_info.reads_info.support_read_no;


	//uniq = 'u', repeat = 'r', none = 'n'
	char type;



	int is_clipped_seq_and_uniq_mapped;
	if ('u' == clipped_align_info.type) { is_clipped_seq_and_uniq_mapped = 2; }
	else if ('r' == clipped_align_info.type) { is_clipped_seq_and_uniq_mapped = 1; }
	else { is_clipped_seq_and_uniq_mapped = 0; return; }
	


	if ('u' == clipped_align_info.type || 'r' == clipped_align_info.type)
	{
		Junction junction;
		SeqInfo up_seq_info, down_seq_info;
		string up_chr, down_chr;
		int up_pos, down_pos;
		if ('+' == clipped_align_info.strand)
		{
			switch(orientation)
			{
			case '5':
				up_chr = clipped_align_info.chr;
				up_pos = clipped_align_info.pos + clipped_align_info.len - 1;
				junction.set_value(up_chr, up_pos, '+', chr, pos, '+');
				up_seq_info.set_value(clipped_seq, clipped_align_info.cigar_vec, clipped_align_info.left_clipped_seq_length, clipped_align_info.right_clipped_seq_length, 0, is_clipped_seq_and_uniq_mapped);
				down_seq_info.set_value(aligned_seq, cigar_vec, 0, 0, support_count, 0);
				break;
			case '3':

				down_chr = clipped_align_info.chr;
				down_pos = clipped_align_info.pos;
				junction.set_value(chr, pos, '+', down_chr, down_pos, '+');
				up_seq_info.set_value(aligned_seq, cigar_vec, 0, 0, support_count, 0);
				down_seq_info.set_value(clipped_seq, clipped_align_info.cigar_vec, clipped_align_info.left_clipped_seq_length, clipped_align_info.right_clipped_seq_length, 0, is_clipped_seq_and_uniq_mapped);
				break;
			}
		}
		else if ('-' == clipped_align_info.strand)
		{
			switch(orientation)
			{
			case '5':
				if (make_pair(clipped_align_info.chr, clipped_align_info.pos) <= make_pair(chr, pos))
				{
					junction.set_value(clipped_align_info.chr, clipped_align_info.pos, '-', chr, pos, '+');
					up_seq_info.set_value(clipped_seq, clipped_align_info.cigar_vec, clipped_align_info.left_clipped_seq_length, clipped_align_info.right_clipped_seq_length, 0, is_clipped_seq_and_uniq_mapped);
					down_seq_info.set_value(aligned_seq, cigar_vec, 0, 0, support_count, 0);
				}
				else
				{
					junction.set_value(chr, pos, '-', clipped_align_info.chr, clipped_align_info.pos, '+');
					GetReverseComplementSeq(aligned_seq);
					GetReverseComplementSeq(clipped_seq);
					ReverseCigar(cigar_vec);
					ReverseCigar(clipped_align_info.cigar_vec);
					up_seq_info.set_value(aligned_seq, cigar_vec, 0, 0, support_count, 0);
					down_seq_info.set_value(clipped_seq, clipped_align_info.cigar_vec, clipped_align_info.right_clipped_seq_length, clipped_align_info.left_clipped_seq_length, 0, is_clipped_seq_and_uniq_mapped);
				}
				break;
			case '3':
				if (make_pair(chr, pos) <= make_pair(clipped_align_info.chr, clipped_align_info.pos + clipped_align_info.len - 1))
				{
					junction.set_value(chr, pos, '+', clipped_align_info.chr, clipped_align_info.pos + clipped_align_info.len - 1, '-');
					up_seq_info.set_value(aligned_seq, cigar_vec, 0, 0, support_count, 0);
					down_seq_info.set_value(clipped_seq, clipped_align_info.cigar_vec, clipped_align_info.left_clipped_seq_length, clipped_align_info.right_clipped_seq_length, 0, is_clipped_seq_and_uniq_mapped);
				}
				else
				{
					junction.set_value(clipped_align_info.chr, clipped_align_info.pos + clipped_align_info.len - 1, '+', chr, pos, '-');
					GetReverseComplementSeq(aligned_seq);
					GetReverseComplementSeq(clipped_seq);
					ReverseCigar(clipped_align_info.cigar_vec);
					ReverseCigar(cigar_vec);
					up_seq_info.set_value(clipped_seq, clipped_align_info.cigar_vec, clipped_align_info.right_clipped_seq_length, clipped_align_info.left_clipped_seq_length, 0, is_clipped_seq_and_uniq_mapped);
					down_seq_info.set_value(aligned_seq, cigar_vec, 0, 0, support_count, 0);
				}
				break;
			}
		}
		else
		{
			cerr << "The alignment orientation maybe error when the sequence is " << clipped_seq
			     << ". The error orientation is " << clipped_align_info.strand << endl;
		}

		pair<multimap<Junction, OtherInfo>::iterator, multimap<Junction, OtherInfo>::iterator> junction2other_it_pair = junction2other.equal_range(junction);
		if (junction2other_it_pair.first == junction2other_it_pair.second)
		{
			OtherInfo other_info(up_seq_info, down_seq_info, -1, 0);
			junction2other.insert(make_pair(junction, other_info));
	
		}
		else
		{
			bool status = 1;
			while (junction2other_it_pair.first != junction2other_it_pair.second) 
			{
				if (junction2other_it_pair.first->second.up_seq_info.right_clipped_seq_length == down_seq_info.left_clipped_seq_length &&  junction2other_it_pair.first->second.down_seq_info.left_clipped_seq_length == up_seq_info.right_clipped_seq_length)
				{
					junction2other_it_pair.first->second.up_seq_info.is_clipped_seq_and_uniq_mapped = junction2other_it_pair.first->second.up_seq_info.is_clipped_seq_and_uniq_mapped > up_seq_info.is_clipped_seq_and_uniq_mapped ? junction2other_it_pair.first->second.up_seq_info.is_clipped_seq_and_uniq_mapped : up_seq_info.is_clipped_seq_and_uniq_mapped;
					junction2other_it_pair.first->second.down_seq_info.is_clipped_seq_and_uniq_mapped = junction2other_it_pair.first->second.down_seq_info.is_clipped_seq_and_uniq_mapped > down_seq_info.is_clipped_seq_and_uniq_mapped ? junction2other_it_pair.first->second.down_seq_info.is_clipped_seq_and_uniq_mapped : down_seq_info.is_clipped_seq_and_uniq_mapped;
	
					junction2other_it_pair.first->second.up_seq_info.support_read_no += up_seq_info.support_read_no;
					junction2other_it_pair.first->second.down_seq_info.support_read_no += down_seq_info.support_read_no;
					if (junction2other_it_pair.first->second.microhomology_length == -1)
						junction2other_it_pair.first->second.microhomology_length = junction2other_it_pair.first->first.up_pos - junction.up_pos;
					status = 0;
				}
			junction2other_it_pair.first++;
			}
			if (status)
			{
				OtherInfo other_info(up_seq_info, down_seq_info, -1, 0);
				junction2other.insert(make_pair(junction, other_info));
			}
		}

	}
	else // type == 'r' or 'n'
	{
		//Store the breakpoint in another map object
		SeqInfo aligned_seq_info(aligned_seq, cigar_vec, 0, 0, support_count, 0);
		ClipReads clip_reads(aligned_seq_info, orientation, clipped_seq, clipped_qual, clipped_align_info.type);
		aligned2clipped.insert(make_pair(make_pair(chr, pos), clip_reads));
	}
}
