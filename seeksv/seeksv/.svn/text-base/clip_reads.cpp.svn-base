/*
 ***************************************
 * soft_reads.h                        *
 *                                     *
 *  Created on:   2011-9-9             *
 *  Modified on:  2012-5-2             *
 *  Modified on:  2012-5-9             *
 *  Version:      0.3.0                *
 *  Author:       qiukl                *
 *                                     *
 ***************************************
 */

#include "clip_reads.h"



ReadsInfo::ReadsInfo(string s_l, string q_l, string s_r, string q_r, const vector<pair<int, char> > &cigar, int quantity, int u)
{
	seq_left = s_l;
	qual_left = q_l;
	seq_right = s_r;
	qual_right = q_r;
	cigar_vec = cigar;
	support_read_no = quantity;
	used = u;
}

/*
ReadsInfo::~ReadsInfo()
{
	vec_materead_id.clear();
}
*/

bool ReadsInfo::set_reads_info(string s_l, string q_l, string s_r, string q_r, const vector<pair<int, char> > &cigar, int quantity)
{
	seq_left = s_l;
	qual_left = q_l;
	seq_right = s_r;
	qual_right = q_r;
	cigar_vec = cigar;
	support_read_no = quantity;
	return 0;
}

/*
bool ReadsInfo::add_materead_id(string id)
{
	vec_materead_id.push_back(id);
	return 0;
}
*/

//aa true means left clipped
bool ReadsInfo::ChangeSeqAndQual(string s_l, string q_l, string s_r, string q_r, const vector<pair<int, char> > &cigar, bool aa)
{
	if (seq_left.length() < s_l.length())
	{
		seq_left = s_l;
		qual_left = q_l;
		if (aa == RIGHT_CLIPPED)
		{
			cigar_vec = cigar;
		}
	}
	if (seq_right.length() < s_r.length())
	{
		seq_right = s_r;
		qual_right = q_r;
		if (aa == LEFT_CLIPPED)
		{
			cigar_vec = cigar;
		}
	}
	return 0;
}


//get the candidate breakpoints and assemble the clipped sequence
void GetSClipReads(const bam_header_t *header, const bam1_t *b, multimap<pair<string, int>, ReadsInfo> &breakpoint2read_l, multimap<pair<string, int>, ReadsInfo> &breakpoint2read_r, double limit, int min_mapQ, bool save_low_quality)
{
	//if this breakpoint hasn't been store ,then store it
	// bam1_strand(b) == 0 means strand is '+'
	int op1 = bam1_cigar(b)[0] & BAM_CIGAR_MASK, op2 = bam1_cigar(b)[b->core.n_cigar - 1] & BAM_CIGAR_MASK;
	//ignore the query with hard clip
	if (op1 == BAM_CHARD_CLIP || op2 == BAM_CHARD_CLIP || b->core.qual < min_mapQ || b->core.flag & BAM_FDUP) return;
	string chr, seq_left, qual_left, seq_right, qual_right, read_id;
	int seq_left_len, seq_right_len, pos, map_len_in_ref;
	int XC_value;
	char XT_type;
	uint8_t *s;
	if (op1 == BAM_CSOFT_CLIP && op2 != BAM_CSOFT_CLIP || op1 != BAM_CSOFT_CLIP && op2 == BAM_CSOFT_CLIP)
	{
		s = bam_aux_get(b, "XC");
		XC_value = bam_aux2i(s);

		if (XC_value != 0 && !save_low_quality) return;

		chr = header->target_name[b->core.tid];
		vector<pair<int, char> > cigar_vec = GenerateCigar(b, map_len_in_ref);
		if (op1 == BAM_CSOFT_CLIP)
		{
			seq_left_len = (bam1_cigar(b)[0] >> BAM_CIGAR_SHIFT);
			seq_right_len = b->core.l_qseq - seq_left_len;
			pos = b->core.pos + 1;
			GetSeq(b, 0, seq_left_len , seq_right_len, seq_left, qual_left, seq_right, qual_right, read_id);
			InsertSeq(breakpoint2read_l, chr, pos, seq_left, qual_left, seq_right, qual_right, cigar_vec, read_id, limit, LEFT_CLIPPED);
		}
		else
		{
			seq_right_len = (bam1_cigar(b)[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);
			seq_left_len = b->core.l_qseq - seq_right_len;
			GetSeq(b, 0, seq_left_len , seq_right_len, seq_left, qual_left, seq_right, qual_right, read_id);
			pos = b->core.pos + map_len_in_ref;
			InsertSeq(breakpoint2read_r, chr, pos, seq_left, qual_left, seq_right, qual_right, cigar_vec, read_id, limit, RIGHT_CLIPPED);
		}
	}
	else if (op1 == BAM_CSOFT_CLIP && op2 == BAM_CSOFT_CLIP)
	{
		seq_left_len = (bam1_cigar(b)[0] >> BAM_CIGAR_SHIFT);
		int seq_right_clipped_len = (bam1_cigar(b)[b->core.n_cigar -1] >> BAM_CIGAR_SHIFT);
		seq_right_len = b->core.l_qseq - seq_left_len - seq_right_clipped_len;
		chr = header->target_name[b->core.tid];
		vector<pair<int, char> > cigar_vec = GenerateCigar(b, map_len_in_ref);

		s = bam_aux_get(b, "XC");
		XC_value = bam_aux2i(s);

		if (XC_value != 0 && !save_low_quality)
		{
			//the left clipped sequence is useful
			if (bam1_strand(b) == 0)
			{
				GetSeq(b, 0, seq_left_len , seq_right_len, seq_left, qual_left, seq_right, qual_right, read_id);
				pos = b->core.pos + 1;
				InsertSeq(breakpoint2read_l, chr, pos, seq_left, qual_left, seq_right, qual_right, cigar_vec, read_id, limit, LEFT_CLIPPED);
			}
			else //the right clipped sequence is useful
			{
				GetSeq(b, seq_left_len, seq_right_len, seq_right_clipped_len, seq_left, qual_left, seq_right, qual_right, read_id);
				pos = b->core.pos + map_len_in_ref;
				InsertSeq(breakpoint2read_r, chr, pos, seq_left, qual_left, seq_right, qual_right, cigar_vec, read_id, limit, RIGHT_CLIPPED);
			}
		}
		else //both clipped sequence is useful
		{
			GetSeq(b, 0, seq_left_len , seq_right_len, seq_left, qual_left, seq_right, qual_right, read_id);
			pos = b->core.pos + 1;
			InsertSeq(breakpoint2read_l, chr, pos, seq_left, qual_left, seq_right, qual_right, cigar_vec, read_id, limit, LEFT_CLIPPED);
			seq_left.clear();
			qual_left.clear();
			seq_right.clear();
			qual_right.clear();
			GetSeq(b, seq_left_len, seq_right_len, seq_right_clipped_len, seq_left, qual_left, seq_right, qual_right, read_id);
			pos = b->core.pos + map_len_in_ref;
			InsertSeq(breakpoint2read_r, chr, pos, seq_left, qual_left, seq_right, qual_right, cigar_vec, read_id, limit, RIGHT_CLIPPED);
		}
	}

}

double CompareStringEndFirst(string str1, string str2)
{
	int len1 = str1.size(), len2 = str2.size();
	int len = len1 < len2 ? len1 : len2;
	int total_match = 0;
	for (int i = 0; i < len; i++)
	{
		if (str1[len1 - 1 - i] == str2[len2 - 1 - i])
			total_match++;
	}
	return (double)total_match/len;
}

double CompareStringBeginFirst(string str1, string str2)
{
	int len = str1.size() < str2.size() ? str1.size() : str2.size();
	int total_match = 0;
	for (int i = 0; i < len; i++)
	{
		if (str1[i] == str2[i])
			total_match++;
	}
	return (double)total_match/len;
}





bool IsUsefulSoftClip(const bam1_t *b)
{
	uint32_t k;
	uint32_t *cigar = bam1_cigar(b);
	int n_cigar = b->core.n_cigar;
	if (n_cigar > 3 || n_cigar <= 1)
		return 0;
	int op1, op2;
	if (!bam1_strand(b))
	{
		op1 = (cigar[0] & BAM_CIGAR_MASK);
		op2 = (cigar[1] & BAM_CIGAR_MASK);
	}
	else
	{
		op1 = (cigar[n_cigar- 1] & BAM_CIGAR_MASK);
		op2 = (cigar[n_cigar - 2] & BAM_CIGAR_MASK);
	}
	if (op1 == BAM_CSOFT_CLIP && op2 == BAM_CMATCH)
		return 1;
	else
		return 0;
}


bool InsertSeq(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, string chr, int pos, string seq_left, string qual_left, string seq_right, string qual_right, vector<pair<int, char> > &cigar, string read_id, double limit, bool aa)
{
	pair<breakpoint2read_iterator, breakpoint2read_iterator> pos_pair = breakpoint2read.equal_range(make_pair(chr, pos));
	while (pos_pair.first != pos_pair.second)
	{
		//if this is true, it means the two seq can be combined to one seq
		if (CompareStringEndFirst(seq_left, pos_pair.first->second.get_seq_left()) >= limit && CompareStringBeginFirst(seq_right, pos_pair.first->second.get_seq_right()) >= limit)
		{
			pos_pair.first->second.ChangeSeqAndQual(seq_left, qual_left, seq_right, qual_right, cigar, aa);
			//pos_pair.first->second.add_materead_id(read_id);
			pos_pair.first->second.support_read_no_increase();
			break;
		}
		++pos_pair.first;
	}
	// if this seq do not have overlap with all seq in the map, insert the seq 
	if (pos_pair.first == pos_pair.second)
	{
		ReadsInfo reads_info(seq_left, qual_left, seq_right, qual_right, cigar, 1, 0);
//		reads_info.add_materead_id(read_id);
		breakpoint2read.insert(make_pair(make_pair(chr, pos), reads_info));
	}
	return 0;
}


void GetSeq(const bam1_t *b, int begin_pos, int seq_left_len ,int seq_right_len, string &seq_left, string &qual_left, string &seq_right, string &qual_right, string &read_id)
{
	seq_left.clear();
	qual_left.clear();
	seq_right.clear();
	qual_right.clear();

	uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
	for (int i = begin_pos; i < seq_left_len + begin_pos; ++i) seq_left.append(1, bam_nt16_rev_table[bam1_seqi(s, i)]);
	for (int i = seq_left_len + begin_pos; i < seq_left_len + seq_right_len + begin_pos; ++i) seq_right.append(1, bam_nt16_rev_table[bam1_seqi(s, i)]);
	if (t[0] == 0xff) { qual_left = "*"; qual_right = "*"; }
	else
	{
		for (int i = begin_pos; i < seq_left_len + begin_pos; ++i) qual_left.append(1, t[i] + 33);
		for (int i = seq_left_len + begin_pos; i < seq_left_len + seq_right_len + begin_pos; ++i) qual_right.append(1, t[i] + 33);
	}
	transform(seq_left.begin(), seq_left.end(), seq_left.begin(), ::toupper);	
	transform(seq_right.begin(), seq_right.end(), seq_right.begin(), ::toupper);	
			
	read_id = bam1_qname(b);
}

//l, such as 10M2I15M1D20M, then l = 46 = 10 + 15 + 1 + 20
vector<pair<int, char> > GenerateCigar(const bam1_t *b, int &l)
{
	l = 0;
	vector<pair<int, char> > cigar_vec;
	uint32_t *cigar = bam1_cigar(b);
	unsigned part_len;
	char buf[16];
	int j;

	for (int i = 0; i < b->core.n_cigar; ++i)
	{
		int op = cigar[i] & BAM_CIGAR_MASK;
		if (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) continue;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CREF_SKIP)
			l += cigar[i] >> BAM_CIGAR_SHIFT;
		part_len = cigar[i]>>BAM_CIGAR_SHIFT;
		
		cigar_vec.push_back(make_pair(part_len, "MIDNSHP=X"[cigar[i]&BAM_CIGAR_MASK]));
	}
	return cigar_vec;
}


//seq2 is 3'clipped  seq4 is 3' aligned
int Compare(const string seq1, const string seq2, const string seq3, const string seq4, const double match_rate)
{
	int len = seq2.length();
	if (len < 10)
		return -1;
	string sub_seq2;
	try {
		sub_seq2 = seq2.substr(0, 10);
	}
	catch (const out_of_range& oor) {
		cerr << "Out of Range error: " << oor.what() << " sub_seq2 = seq2.substr(0, 10); in Compare" << endl;
	}
	int pos = seq4.find(sub_seq2);
	
	if (pos == string::npos)
		return -1;
	else
	{
		string seq5 = seq3;
		try {
			seq5.append(seq4.substr(0, pos));
		}
		catch (const out_of_range& oor) {
			cerr << "Out of Range error: " << oor.what() << "seq5.append(seq4.substr(0, pos)) in Compare()" << endl;
		}
		string seq6;
		try {
			seq6 = seq4.substr(pos);
		}
		catch (const out_of_range& oor) {
			cerr << "Out of Range error: " << oor.what() << "seq6 = seq4.substr(pos); in Compare()" << endl;
		}
		if (CompareStringEndFirst(seq1, seq5) >= match_rate && CompareStringBeginFirst(seq2, seq6) >= match_rate)
		{
			return pos;
		}
		else
			return -1;
	}
}


bool GetSeqAndQual(const bam1_t *b, string &seq, string &qual)
{
	uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
	seq.clear();
	qual.clear();
	if (b->core.l_qseq)
	{
		for (int i = 0; i < b->core.l_qseq; ++i) seq.append(1, bam_nt16_rev_table[bam1_seqi(s, i)]);
		if (t[0] == 0xff) qual = "*";
		else for (int i = 0; i < b->core.l_qseq; ++i) qual.append(1, t[i] + 33);
	}

	return 0;
}

/*
bool DisplayUnmapFastq(map<string, SeqPair> &id2seq_qual, string fq1, string fq2)
{
	ofstream fout1(fq1.c_str()), fout2(fq2.c_str());
	map<string, SeqPair>::iterator map_it = id2seq_qual.begin();
	while (map_it != id2seq_qual.end())
	{
		if (!map_it->second.get_seq1().empty() && !map_it->second.get_seq2().empty())
		{
			fout1 << "@" << map_it->first << "/1" << '\n'
				  << map_it->second.get_seq1() << '\n'
				  << "+\n"
				  << map_it->second.get_qual1() << endl;

			fout2 << "@" << map_it->first << "/2\n"
				  << map_it->second.get_seq2() << '\n'
				  << "+\n"
				  << map_it->second.get_qual2() << endl; 
		}
		++map_it;
	}
}
*/

void GetReverseComplementSeq(string &seq)
{
	char c;
	int l = seq.size();
	for(int i = 0; i < l / 2; i++)
	{
		c = seq[i];
		switch (c)
		{
		case 'A': c = 'T'; break;
		case 'a': c = 'T'; break;
		case 'T': c = 'A'; break;
		case 't': c = 'A'; break;
		case 'C': c = 'G'; break;
		case 'c': c = 'G'; break;
		case 'G': c = 'C'; break;
		case 'g': c = 'C'; break;
		case 'N': c = 'N'; break;
		case 'n': c = 'N';
		}
		seq[i] = seq[l-i-1];
		switch (seq[i])
		{
		case 'A': seq[i] = 'T'; break;
		case 'a': seq[i] = 'T'; break;
		case 'T': seq[i] = 'A'; break;
		case 't': seq[i] = 'A'; break;
		case 'C': seq[i] = 'G'; break;
		case 'c': seq[i] = 'G'; break;
		case 'G': seq[i] = 'C'; break;
		case 'g': seq[i] = 'C'; break;
		case 'N': seq[i] = 'N'; break;
		case 'n': seq[i] = 'N';
		}
		seq[l-i-1] = c;
	}
	if (l % 2 == 1)
	{
		switch (seq[l/2])
		{
		case 'A': seq[l/2] = 'T'; break;
		case 'a': seq[l/2] = 'T'; break;
		case 'T': seq[l/2] = 'A'; break;
		case 't': seq[l/2] = 'A'; break;
		case 'C': seq[l/2] = 'G'; break;
		case 'c': seq[l/2] = 'G'; break;
		case 'G': seq[l/2] = 'C'; break;
		case 'g': seq[l/2] = 'C'; break;
		case 'N': seq[l/2] = 'N'; break;
		case 'n': seq[l/2] = 'N';
		}
	}
}
bool MinusCigarLeft(vector<pair<int, char> > &cigar_vec, int length)
{
	int len = 0;
	vector<pair<int, char> >::iterator vec_it = cigar_vec.begin();
	while (vec_it != cigar_vec.end())
	{
		if ((*vec_it).second == 'M' || (*vec_it).second == 'I')
			len += (*vec_it).first;
		++vec_it;
	}
	if (len <= length)
	{
		return 0;
	}
	else
	{
		len = length;
		vec_it = cigar_vec.begin();
		while (vec_it != cigar_vec.end())
		{
			if ((*vec_it).second == 'M' || (*vec_it).second == 'I')
			{
				if ((*vec_it).first > len)
				{
					(*vec_it).first -= len;
					break;
				}
				else
				{
					len -= (*vec_it).first;
					cigar_vec.erase(vec_it++);
				}
			}
			else
				cigar_vec.erase(vec_it++);
		}
		return 1;
	}
}

bool MinusCigarRight(vector<pair<int, char> > &cigar_vec, int length)
{
	int len = 0;
	vector<pair<int, char> >::iterator vec_it = cigar_vec.begin();
	while (vec_it != cigar_vec.end())
	{
		if ((*vec_it).second == 'M' || (*vec_it).second == 'I')
			len += (*vec_it).first;
		++vec_it;
	}
	if (len <= length)
	{
		return 0;
	}
	else
	{
		int l = len - length;
		vec_it = cigar_vec.begin();
		while (vec_it != cigar_vec.end())
		{
			if ((*vec_it).second == 'M' || (*vec_it).second == 'I')
			{
				if ((*vec_it).first >= l)
				{
					(*vec_it).first = l;
					cigar_vec.erase(++vec_it, cigar_vec.end());
					break;
				}
				else
				{
					l -= (*vec_it).first;
					++vec_it;
				}
			}
			else
				++vec_it;
		}
		return 1;
	}
}

bool AddCigarLeft(vector<pair<int, char> > &cigar_vec, int length)
{
	if (cigar_vec[0].second == 'M')
	{
		cigar_vec[0].first += length;
	}
	else
	{
		cigar_vec.insert(cigar_vec.begin(), make_pair(length, 'M'));
	}
}

bool AddCigarRight(vector<pair<int, char> > &cigar_vec, int length)
{
	if (cigar_vec[cigar_vec.size() - 1].second == 'M')
	{
		cigar_vec[cigar_vec.size() - 1].first += length;
	}
	else
	{
		cigar_vec.push_back(make_pair(length, 'M'));
	}
}

string GetSVType(string up_chr, int up_pos, char up_strand, string down_chr, int down_pos, char down_strand)
{
	string sv_type;
	if (up_chr != down_chr) { sv_type = "CTX"; }
	else if (up_strand != down_strand) { sv_type = "INV"; }
	else if (up_pos < down_pos) { sv_type = "DEL"; }
	else if (up_pos > down_pos) { sv_type = "INS"; }
	else { sv_type = "Unknown"; }
	return sv_type;
}



void MergeRightClippedSoftClippedReads(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, int search_length, double min_match_rate)
{
	multimap<pair<string, int>, ReadsInfo>::iterator breakpoint2read_it = breakpoint2read.end(), breakpoint2read_it1, breakpoint2read_it2;
	--breakpoint2read_it;
	while (breakpoint2read_it != breakpoint2read.begin())
	{
		breakpoint2read_it1 = breakpoint2read_it;
		while (breakpoint2read_it1-- != breakpoint2read.begin() && breakpoint2read_it1->first.first == breakpoint2read_it->first.first && breakpoint2read_it1->first.second > breakpoint2read_it->first.second - search_length)
		{
			
			pair<string, int> chr_pos = breakpoint2read_it1->first;
			ReadsInfo reads_info = breakpoint2read_it1->second;
			int shift_length = breakpoint2read_it->first.second - chr_pos.second;
			if (reads_info.seq_right.length() >= shift_length + 5)
			{
				reads_info.seq_left = reads_info.seq_left.append(reads_info.seq_right, 0, shift_length);
				try {
					reads_info.seq_right = reads_info.seq_right.substr(shift_length);
				}
				catch (const out_of_range& oor) {
					cerr << "Out of Range error: " << oor.what() << " reads_info.seq_right = reads_info.seq_right.substr(shift_length);"  << endl;
				}
				
				reads_info.qual_left = reads_info.qual_left.append(reads_info.qual_right, 0, shift_length);
				try {
					reads_info.qual_right = reads_info.qual_right.substr(shift_length);
				}
				catch (const out_of_range& oor) {
					cerr << "Out of Range error: " << oor.what() << " reads_info.qual_right = reads_info.qual_right.substr(shift_length);" << endl;
				}
				if (CompareStringEndFirst(reads_info.seq_left, breakpoint2read_it->second.seq_left) >= min_match_rate && 
					CompareStringBeginFirst(reads_info.seq_right, breakpoint2read_it->second.seq_right) >= min_match_rate)
				{
					string temp;
					if (reads_info.seq_left.length() > breakpoint2read_it->second.seq_left.length())
					{
						int add_length = reads_info.seq_left.length() - breakpoint2read_it->second.seq_left.length();
						try {
							temp = reads_info.seq_left.substr(0, add_length);
						}
						catch (const out_of_range& oor) {
							cerr << "Out of Range error: " << oor.what() << " temp = reads_info.seq_left.substr(0, add_length);" << endl;
						}
						temp.append(breakpoint2read_it->second.seq_left);
						breakpoint2read_it->second.seq_left = temp;
						try {
							temp = reads_info.qual_left.substr(0, add_length);
						}
						catch (const out_of_range& oor) {
							cerr << "Out of Range error: " << oor.what() << " temp = reads_info.qual_left.substr(0, add_length);" <<endl;
						}
						temp.append(breakpoint2read_it->second.qual_left);
						breakpoint2read_it->second.qual_left = temp;
						AddCigarLeft(breakpoint2read_it->second.cigar_vec, add_length);
					}
					if (reads_info.seq_right.length() > breakpoint2read_it->second.seq_right.length())
					{
						int add_length = reads_info.seq_right.length() - breakpoint2read_it->second.seq_right.length();
						breakpoint2read_it->second.seq_right.append(reads_info.seq_right, breakpoint2read_it->second.seq_right.length(), add_length);
						breakpoint2read_it->second.qual_right.append(reads_info.qual_right, breakpoint2read_it->second.qual_right.length(), add_length);
					}
					breakpoint2read_it->second.support_read_no += reads_info.support_read_no;

					breakpoint2read_it2 = breakpoint2read_it1;
					++breakpoint2read_it2;
					breakpoint2read.erase(breakpoint2read_it1);
					breakpoint2read_it1 = breakpoint2read_it2;
				}
			}
			else
			{
				int overlap_length = breakpoint2read_it->second.seq_left.length() - shift_length;
				if (overlap_length >= 0)
				{
					string seq_right;
					try {
						seq_right = breakpoint2read_it->second.seq_left.substr(overlap_length);
					}
					catch (const out_of_range& oor) {
						cerr << "Out of Range error: " << oor.what() << " seq_right = breakpoint2read_it->second.seq_left.substr(overlap_length);" << endl;
					}
					seq_right.append(breakpoint2read_it->second.seq_right);

					if (CompareStringBeginFirst(reads_info.seq_right, seq_right) >= min_match_rate)
					{
						breakpoint2read_it2 = breakpoint2read_it1;
						++breakpoint2read_it2;
						breakpoint2read.erase(breakpoint2read_it1);
						breakpoint2read_it1 = breakpoint2read_it2;
					}
				}
			}
			
		}
		//if erase the first element of breakpoint2read , breakpoint2read_it will be equal to breakpoint2read.begin()
		if (breakpoint2read_it == breakpoint2read.begin())
			break;
		
		--breakpoint2read_it;
	}
}


void MergeLeftClippedSoftClippedReads(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, int search_length, double min_match_rate)
{
	multimap<pair<string, int>, ReadsInfo>::iterator breakpoint2read_it = breakpoint2read.begin(), breakpoint2read_it1;
	while (breakpoint2read_it != breakpoint2read.end())
	{
		breakpoint2read_it1 = breakpoint2read_it;
		++breakpoint2read_it1;
		while (breakpoint2read_it1 != breakpoint2read.end() && breakpoint2read_it1->first.first == breakpoint2read_it->first.first && breakpoint2read_it1->first.second < breakpoint2read_it->first.second + search_length)
		{
			pair<string, int> chr_pos = breakpoint2read_it1->first;
			ReadsInfo reads_info = breakpoint2read_it1->second;
			int shift_length = chr_pos.second - breakpoint2read_it->first.second;
			//for situation that there are some small INDELs near the breakpoint,such as 
			//chr17   77656162        5       15M2D66M        AAAGAAAGAAGCTGTGCAGCACTGGACTGGCAGAGAGACCCATGGGACAGAACAGAGCCAGGGCAGCCCCATGCACTGATG	AAAAAAAGAAGAAGAAAGAGAAAGACAAAGAAAGAAGACAAAGAAAG
			//chr17   77656179        5       83M     GCAGCACTGGACTGGCAGAGAGACCCATGGGACAGAACAGAGCCAGGGCAGCCCCATGCACTGATGTGAGCTTGGTCATGACC	AAGCTGT
			if (breakpoint2read_it->second.cigar_vec.size() > 1)
			{
				int len = 0, coordinate_diff = shift_length;
				for (int i = 0; i < breakpoint2read_it->second.cigar_vec.size(); i++)
				{
					if (breakpoint2read_it->second.cigar_vec[i].second == 'M')
					{
						len += breakpoint2read_it->second.cigar_vec[i].first;
						if (len >= coordinate_diff) break;
					}
					else if (breakpoint2read_it->second.cigar_vec[i].second == 'D')
					{
						len += breakpoint2read_it->second.cigar_vec[i].first;
						if (len >= coordinate_diff)
						{
							shift_length -= breakpoint2read_it->second.cigar_vec[i].first;
							break;
						}
						else
							shift_length -= breakpoint2read_it->second.cigar_vec[i].first;
					}
					else if (breakpoint2read_it->second.cigar_vec[i].second == 'I')
					{
						shift_length += breakpoint2read_it->second.cigar_vec[i].first;
					}
				}
			}
			if (shift_length < 0)
			{
				++breakpoint2read_it1;
				continue;
			}

			if (reads_info.seq_left.length() >= shift_length + 5)
			{
				string temp = reads_info.seq_right;
				try {
					reads_info.seq_right = reads_info.seq_left.substr(reads_info.seq_left.length() - shift_length);
				}
				catch (const out_of_range& oor) {
					cerr << "Out of Range error: " << oor.what() << " reads_info.seq_right = reads_info.seq_left.substr(reads_info.seq_left.length() - shift_length);" << endl;
				}
				reads_info.seq_right.append(temp);
				temp = reads_info.qual_right;
				try {
					reads_info.qual_right = reads_info.qual_left.substr(reads_info.qual_left.length() - shift_length);
				}
				catch (const out_of_range& oor) {
					cerr << "Out of Range error: " << oor.what() << " reads_info.qual_right = reads_info.qual_left.substr(reads_info.qual_left.length() - shift_length);" << endl;
				}
				reads_info.qual_right.append(temp);
				try {
					reads_info.seq_left = reads_info.seq_left.substr(0, reads_info.seq_left.length() - shift_length);
				}
				catch (const out_of_range& oor) {
					cerr << "Out of Range error: " << oor.what() << " reads_info.seq_left = reads_info.seq_left.substr(0, reads_info.seq_left.length() - shift_length);" << endl;
				}
				try {
					reads_info.qual_left = reads_info.qual_left.substr(0, reads_info.qual_left.length() - shift_length);
				}
				catch (const out_of_range& oor) {
					cerr << "Out of Range error: " << oor.what() << " reads_info.qual_left = reads_info.qual_left.substr(0, reads_info.qual_left.length() - shift_length);" << endl;
				}
				

				if (CompareStringEndFirst(reads_info.seq_left, breakpoint2read_it->second.seq_left) >= min_match_rate && 
					CompareStringBeginFirst(reads_info.seq_right, breakpoint2read_it->second.seq_right) >= min_match_rate)
				{
					if (breakpoint2read_it->second.seq_left.length() < reads_info.seq_left.length())
					{
						breakpoint2read_it->second.seq_left = reads_info.seq_left;
						breakpoint2read_it->second.qual_left = reads_info.qual_left;
					}
					if (breakpoint2read_it->second.seq_right.length() < reads_info.seq_right.length())
					{
						int add_length = reads_info.seq_right.length() - breakpoint2read_it->second.seq_right.length();
						breakpoint2read_it->second.seq_right.append(reads_info.seq_right, reads_info.seq_right.length()-add_length, add_length);
						breakpoint2read_it->second.qual_right.append(reads_info.qual_right, reads_info.qual_right.length()-add_length, add_length);
						AddCigarRight(breakpoint2read_it->second.cigar_vec, add_length);
					}
					breakpoint2read_it->second.support_read_no += reads_info.support_read_no;
					breakpoint2read.erase(breakpoint2read_it1++);
				}
				else
					++breakpoint2read_it1;
			}
			else
			{
				int overlap_length = breakpoint2read_it->second.seq_right.length() - shift_length;
				if (overlap_length >= 0)
				{
					string seq_left = breakpoint2read_it->second.seq_left;
					seq_left.append(breakpoint2read_it->second.seq_right, 0, shift_length);

					if (CompareStringEndFirst(reads_info.seq_left, seq_left) >= min_match_rate)
					{
						breakpoint2read.erase(breakpoint2read_it1++);
					}
					else
						++breakpoint2read_it1;
				}
				else
					++breakpoint2read_it1;
			}
		}
		++breakpoint2read_it;
	}
}

void GenerateKmer(multimap<pair<string, int>, ReadsInfo> &breakpoint2read, multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator, int> > &kmer2interator, int kmer_len, int left_or_right_clip)
{
	multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = breakpoint2read.begin();

	while (mmap_it != breakpoint2read.end())
	{
		// Add '&& mmap_it->second.support_read_no >= MIN_READ_NO' in revision 37, only soft-clipped reads with support reads number >= MIN_READ_NO(2) will be use to generate kmer.
		if (left_or_right_clip == RIGHT_CLIPPED && mmap_it->second.seq_right.length() >= kmer_len || left_or_right_clip == LEFT_CLIPPED && mmap_it->second.seq_left.length() >= kmer_len && mmap_it->second.support_read_no >= MIN_READ_NO)
		{
			string seq;
			if (left_or_right_clip == RIGHT_CLIPPED)
			{
				seq = mmap_it->second.seq_left;
				string subseq;
				for (int i = 0; i<= 9; i ++)
				{
					int start_pos = (int)seq.length() - kmer_len - i;
					if (start_pos < 0) break;
					try {
						subseq = seq.substr(start_pos, kmer_len);
					}
					catch (const out_of_range& oor) {
						cerr << "Out of Range error: " << oor.what() << " subseq = seq.substr(i, kmer_len) in GenerateKmer" << endl;
					}
					kmer2interator.insert(make_pair(subseq, make_pair(mmap_it, start_pos)));
				}
			}
			else if (left_or_right_clip == LEFT_CLIPPED)
			{
				seq = mmap_it->second.seq_right;
				string subseq;
				for (int i = 0; i + kmer_len <= seq.length() && i<= 9; i ++)
				{
					try {
						subseq = seq.substr(i, kmer_len);
					}
					catch (const out_of_range& oor) {
						cerr << "Out of Range error: " << oor.what() << " subseq = seq.substr(i, kmer_len) in GenerateKmer" << endl;
					}
					kmer2interator.insert(make_pair(subseq, make_pair(mmap_it, i)));
				}
			}
			else
			{
				cerr << "Error of argument int left_or_right_clip in function GenerateKmer()" << endl;
				exit(1);
			}
		}
		++mmap_it;
	}
}


void CompareKmer(ofstream &fout, multimap<pair<string, int>, ReadsInfo> &breakpoint2read, multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator, int> > &kmer2interator, int kmer_len, int type)
{
	double rate_threshold = 0.9;
	pair<multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator,int> >::iterator, multimap<string, pair<multimap<pair<string, int>, ReadsInfo>::iterator, int> >::iterator> iter_pair;

	multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = breakpoint2read.begin(), stored_mmap_it;
	while (mmap_it != breakpoint2read.end())
	{
		if (mmap_it->second.support_read_no < MIN_READ_NO) {
			++mmap_it;
			continue;
		}

		string current_seq_left, current_seq_right;
		vector<pair<int, char> > current_cigar_vec = mmap_it->second.cigar_vec;


		if (type == SAME_STRAND)
		{
			current_seq_left = mmap_it->second.seq_left;
			current_seq_right = mmap_it->second.seq_right;
		}
		else
		{
			current_seq_left = mmap_it->second.seq_right;
			current_seq_right = mmap_it->second.seq_left;
			GetReverseComplementSeq(current_seq_left);	
			GetReverseComplementSeq(current_seq_right);
		}
		string current_seq = current_seq_left;
		current_seq.append(current_seq_right);
		
		string subseq;
		double largest_match_rate = 0, second_largest_match_rate = 0, match_rate = 0;
		int overlap_len = 0, same_len = 0, distance = 0;
		if (type == SAME_STRAND || type == DIFF_STRAND_RIGHT_CLIP)
		{
			if (current_seq_left.length() < kmer_len)
			{
				++mmap_it;
				continue;
			}
			try {
				subseq = current_seq_left.substr(current_seq_left.length() - kmer_len, kmer_len);
			}
			catch (const out_of_range& oor) {
				cerr << "Out of Range error: " << oor.what() << " subseq = current_seq_left.substr(current_seq_left.length() - kmer_len, kmer_len);" << endl;
			}
			iter_pair = kmer2interator.equal_range(subseq);
			while (iter_pair.first != iter_pair.second)
			{
				int pos1 = iter_pair.first->second.second;
				string stored_seq_left = iter_pair.first->second.first->second.seq_left;
				string stored_seq_right = iter_pair.first->second.first->second.seq_right;
				string stored_seq = stored_seq_left;
				stored_seq.append(stored_seq_right);
				int pos2 = current_seq_left.length() - kmer_len;
				if (pos2 >= pos1)
				{
					int dis = pos2 - pos1;
					for (int i = 0; i < stored_seq.length() && i + dis < current_seq.length(); i++)
					{
						if (stored_seq[i] == current_seq[i+dis])
						{
							++same_len;
						}
						++overlap_len;
					}
				}
				else
				{
					int dis = pos1 - pos2;
					for (int i = 0; i + dis < stored_seq.length() && i < current_seq.length(); i++)
					{
						if (stored_seq[i+dis] == current_seq[i])
						{
							++same_len;
						}
						++overlap_len;
					}
				}
				match_rate = (double)same_len/overlap_len;
				if (largest_match_rate <= match_rate)
				{
					second_largest_match_rate = largest_match_rate;
					largest_match_rate = match_rate;
					stored_mmap_it = iter_pair.first->second.first;
					distance = pos2 - pos1;
				}
				else if (second_largest_match_rate < match_rate)
				{
					second_largest_match_rate = match_rate;
				}
				iter_pair.first++;
			}
		}
		else if (type == DIFF_STRAND_LEFT_CLIP)
		{
			if (current_seq_right.length() < kmer_len)
			{
				++mmap_it;
				continue;
			}
			try {
				subseq = current_seq_right.substr(0, kmer_len);
			}
			catch (const out_of_range& oor) {
				cerr << "Out of Range error: " << oor.what() << " subseq = current_seq_left.substr(current_seq_left.length() - kmer_len, kmer_len);" << endl;
			}
			iter_pair = kmer2interator.equal_range(subseq);
			while (iter_pair.first != iter_pair.second)
			{
				int pos1 = iter_pair.first->second.second;
				string stored_seq_left = iter_pair.first->second.first->second.seq_left;
				string stored_seq_right = iter_pair.first->second.first->second.seq_right;
				string stored_seq = stored_seq_left;
				stored_seq.append(stored_seq_right);
				pos1 += stored_seq_left.length();
				int pos2 = current_seq_left.length();
				if (pos2 >= pos1)
				{
					int dis = pos2 - pos1;
					for (int i = 0; i < stored_seq.length() && i + dis < current_seq.length(); i++)
					{
						if (stored_seq[i] == current_seq[i+dis])
						{
							++same_len;
						}
						++overlap_len;
					}
				}
				else
				{
					int dis = pos1 - pos2;
					for (int i = 0; i + dis < stored_seq.length() && i < current_seq.length(); i++)
					{
						if (stored_seq[i+dis] == current_seq[i])
						{
							++same_len;
						}
						++overlap_len;
					}
				}
				match_rate = (double)same_len/overlap_len;
				if (largest_match_rate <= match_rate)
				{
					second_largest_match_rate = largest_match_rate;
					largest_match_rate = match_rate;
					stored_mmap_it = iter_pair.first->second.first;
					distance = pos2 - pos1;
				}
				else if (second_largest_match_rate < match_rate)
				{
					second_largest_match_rate = match_rate;
				}
				iter_pair.first++;
			}

		}
		else
		{
		}

		if (largest_match_rate == second_largest_match_rate || largest_match_rate < 0.9)
		{
			++mmap_it;
			continue;
		}

		int stored_seq_left_len = stored_mmap_it->second.seq_left.length();
		int stored_seq_right_len = stored_mmap_it->second.seq_right.length();
		int current_seq_left_len = current_seq_left.length();
		int current_seq_right_len = current_seq_right.length();

		int microhomology_len = stored_seq_left_len - current_seq_left_len + distance;
		if (type == SAME_STRAND || (type == DIFF_STRAND_RIGHT_CLIP && stored_mmap_it->first < mmap_it->first))
		{

			if (microhomology_len >= 0 && current_seq_right_len >= microhomology_len + 5 && stored_seq_left_len >=  microhomology_len + 5)
			{
				string seq_right = current_seq_right;
				int add_length = stored_seq_right_len + microhomology_len - current_seq_right_len;
				if (add_length > 0)
				{
					seq_right.append(stored_mmap_it->second.seq_right, stored_seq_right_len - add_length, add_length);
					AddCigarRight(current_cigar_vec, add_length);
				}
				vector<pair<int, char> > stored_cigar_vec = stored_mmap_it->second.cigar_vec;
				MinusCigarRight(stored_cigar_vec, microhomology_len);
				string seq_left = current_seq_left;

				add_length = stored_seq_left_len - (int)current_seq_left.length() - microhomology_len;
				if (add_length > 0)
				{
					string temp = seq_left;
					try {
						seq_left = stored_mmap_it->second.seq_left.substr(0, add_length);
					}
					catch (const out_of_range& oor) {
						 cerr << "Out of Range error: " << oor.what() << " seq_left = stored_mmap_it->second.seq_left.substr(0, add_length) in CompareKmer()" << endl;
					}
					seq_left.append(temp);
					//AddCigarLeft(stored_cigar_vec, add_length);
				}
				else
				{
					add_length = -add_length;
					AddCigarLeft(stored_cigar_vec, add_length);
				}
				char up_strand = '+';
				char down_strand = '+';
				if (type == DIFF_STRAND_RIGHT_CLIP) down_strand = '-';
				if (microhomology_len <= 10)
				{
					stored_mmap_it->second.used = 1;
					mmap_it->second.used = 1;
					string sv_type = GetSVType(stored_mmap_it->first.first, stored_mmap_it->first.second - microhomology_len, up_strand, mmap_it->first.first, mmap_it->first.second, down_strand); 
					fout << stored_mmap_it->first.first << '\t' << stored_mmap_it->first.second - microhomology_len << '\t' << up_strand << '\t' << stored_mmap_it->second.support_read_no << '\t' << mmap_it->first.first << '\t' << mmap_it->first.second << '\t' << down_strand << '\t' << mmap_it->second.support_read_no << '\t' << microhomology_len << '\t' << 0 << '\t' << sv_type << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0<< '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t';
					DisplayCigarVector<ofstream> (fout, stored_cigar_vec, 0, 0);
					fout << '\t';
					DisplayCigarVector<ofstream> (fout, current_cigar_vec, 0, 0);
					fout << '\t';
					fout << seq_left << '\t' << seq_right << endl;
				}
			}
		}
		else if (type == DIFF_STRAND_LEFT_CLIP && mmap_it->first < stored_mmap_it->first)
		{
			microhomology_len = -microhomology_len;
			if (microhomology_len >= 0 && current_seq_left_len >= microhomology_len + 5 && stored_seq_right_len >=  microhomology_len + 5)
			{
				string seq_left = current_seq_left;
				int add_length = stored_seq_left_len + microhomology_len - current_seq_left_len;
				if (add_length > 0)
				{
					string temp = seq_left;
					try {
						seq_left = stored_mmap_it->second.seq_left.substr(0, add_length);
					}
					catch (const out_of_range& oor) {
						cerr << "Out of Range error: " << oor.what() << " seq_left = stored_mmap_it->second.seq_left.substr(0, add_length);  in CompareKmer()" << endl;
					}
					seq_left.append(temp);
					AddCigarLeft(current_cigar_vec, add_length);
				}

				vector<pair<int, char> > stored_cigar_vec = stored_mmap_it->second.cigar_vec;
				string seq_right = current_seq_right;
				add_length = stored_seq_right_len - microhomology_len - current_seq_right_len;
				MinusCigarLeft(stored_cigar_vec, microhomology_len);
				if (add_length > 0)
				{
					seq_right.append(stored_mmap_it->second.seq_right, stored_seq_right_len - add_length, add_length);
				}
				else
				{
					add_length = -add_length;
					AddCigarRight(stored_cigar_vec, add_length); 
				}
				char up_strand = '-';
				char down_strand = '+';
				if (microhomology_len <= 10)
				{
					stored_mmap_it->second.used = 1;
					mmap_it->second.used = 1;
					string sv_type = GetSVType(mmap_it->first.first, mmap_it->first.second, up_strand, stored_mmap_it->first.first, stored_mmap_it->first.second + microhomology_len, down_strand);
					fout <<  mmap_it->first.first << '\t' << mmap_it->first.second << '\t' << up_strand << '\t' << mmap_it->second.support_read_no << '\t' << stored_mmap_it->first.first << '\t' << stored_mmap_it->first.second + microhomology_len << '\t' << down_strand << '\t' << stored_mmap_it->second.support_read_no << '\t' << microhomology_len << '\t' << 0 << '\t' << sv_type << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t';
					DisplayCigarVector<ofstream> (fout, current_cigar_vec, 0, 0);
					fout << '\t';
					DisplayCigarVector<ofstream> (fout, stored_cigar_vec, 0, 0);
					fout << '\t';
					fout << seq_left << '\t' << seq_right << endl;
				}
			}
		}
		++mmap_it;
	}
}


int FindDistance(string &seq1, string &seq2, int kmer_len)
{
	multimap<string, int> subseq2tag;
	map<int, int> distance2no;
	multimap<string, int>::iterator subseq2tag_it;
	map<int, int>::iterator distance2no_it;

	string subseq;
	int distance;
	for (int i = 0; i + kmer_len <= seq1.length(); i ++)
	{
		try {
			subseq = seq1.substr(i, kmer_len);
		}
		catch (const out_of_range& oor) {
			cerr << "Out of Range error: " << oor.what()  << " subseq = seq1.substr(i, kmer_len); int FindDistance()" << endl;
		}
		subseq2tag.insert(make_pair(subseq, i));
	}
	for (int i = 0; i + kmer_len <= seq2.length(); i ++)
	{
		try {
			subseq = seq2.substr(i, kmer_len);
		}
		catch (const out_of_range& oor) {
			cerr << "Out of Range error: " << oor.what()  << " subseq = seq2.substr(i, kmer_len); int FindDistance()" << endl;
		}
		if (subseq2tag.count(subseq) == 1)
		{
			subseq2tag_it = subseq2tag.find(subseq);
			distance = i - subseq2tag_it->second;
			distance2no_it = distance2no.find(distance);
			if (distance2no_it != distance2no.end())
				distance2no_it->second = distance2no_it->second + 1;
			else
				distance2no.insert(make_pair(distance, 1));
			
		}
		distance2no_it = distance2no.begin();
		int large_no = 0;
		while (distance2no_it != distance2no.end())
		{
			if (distance2no_it->second > large_no)
			{
				large_no = distance2no_it->second;
				distance = distance2no_it->first;
			}
			++ distance2no_it;
		}
	}
	return distance;
}
