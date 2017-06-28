#include "process_bwasw.h"
#include "clip_reads.h"
#include "getsv.h"

void FindJunction(string file, int min_mapQ, multimap<Junction, OtherInfo> &junction2other)
//void FindJunction(string file, int min_mapQ)
{
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

	g_min_mapQ = min_mapQ;

	string chr, read_id, left_seq, right_seq, left_qual, right_qual;
	int pos, left_seq_len, right_seq_len, map_len_in_ref;
	char clipped_side, strand;
	vector<pair<int, char> > cigar_vec;
	//GenerateCigar(b, map_len_in_ref)

	map<string, Alignment> read_id2align;
	map<string, Alignment>::iterator read_id2align_it;

	//read a line of the bam file, if the line is without soft-clipping read, ignore it
	while ((ret = samread(samfin, b) >= 0))
	{
		// when read failed continue
		if (__g_skip_aln(samfin->header, b)) continue;
		//if the reads is unmap
		if (b->core.flag & BAM_FUNMAP) continue;
		int op1 = bam1_cigar(b)[0] & BAM_CIGAR_MASK, op2 = bam1_cigar(b)[b->core.n_cigar - 1] & BAM_CIGAR_MASK;
		if (op1 == BAM_CHARD_CLIP || op2 == BAM_CHARD_CLIP || (op1 == BAM_CSOFT_CLIP && op2 == BAM_CSOFT_CLIP) || (op1 == BAM_CMATCH && op2 == BAM_CMATCH) || b->core.flag & BAM_FDUP) continue;

		cigar_vec = GenerateCigar(b, map_len_in_ref);
		if (op1 == BAM_CSOFT_CLIP)
		{
			clipped_side = '5';
			left_seq_len = (bam1_cigar(b)[0] >> BAM_CIGAR_SHIFT);
			right_seq_len = b->core.l_qseq - left_seq_len;
			pos = b->core.pos + 1;
		}
		else //op2 == BAM_CSOFT_CLIP
		{
			clipped_side = '3';
			right_seq_len = (bam1_cigar(b)[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT);
			left_seq_len = b->core.l_qseq - right_seq_len;
			pos = b->core.pos + map_len_in_ref;
		}
		//bam1_strand(b) == 0 means strand is '+'
		bam1_strand(b) ? strand = '-' : strand = '+';
		chr = samfin->header->target_name[b->core.tid];
		GetSeq(b, 0, left_seq_len, right_seq_len, left_seq, left_qual, right_seq, right_qual, read_id);

		Alignment alignment, up_alignment, down_alignment;

		if (clipped_side == '5')
		{
			alignment.set_value(chr, pos, left_seq, left_qual, right_seq, right_qual, cigar_vec, '5', strand);
		}
		else
		{
			alignment.set_value(chr, pos, left_seq, left_qual, right_seq, right_qual, cigar_vec, '3', strand);
		}
		
		read_id2align_it = read_id2align.find(read_id);


		if (read_id2align_it != read_id2align.end())
		{
			int microhomology_length = -1;
			if (read_id2align_it->second.strand == strand && read_id2align_it->second.clipped_side != clipped_side || read_id2align_it->second.strand != strand && read_id2align_it->second.clipped_side == clipped_side)
			{
				Junction junction;
				SeqInfo up_seq_info, down_seq_info;
				if (read_id2align_it->second.strand == strand && read_id2align_it->second.clipped_side != clipped_side)
				{

					switch(read_id2align_it->second.clipped_side)
					{
					case '5':
						up_alignment = alignment;
						down_alignment = read_id2align_it->second;
						break;
					case '3':
						up_alignment = read_id2align_it->second;
						down_alignment = alignment;
						break;
					}

					if (up_alignment.left_seq.length() >= down_alignment.left_seq.length())
					{
						microhomology_length = up_alignment.left_seq.length() - down_alignment.left_seq.length();
						junction.set_value(up_alignment.chr, up_alignment.pos - microhomology_length, '+', down_alignment.chr, down_alignment.pos, '+');

						cigar_vec = up_alignment.cigar_vec;
						MinusCigarRight(cigar_vec, microhomology_length);

						up_seq_info.set_value(down_alignment.left_seq, cigar_vec, 0, 0, 0, 2);
						down_seq_info.set_value(down_alignment.right_seq, down_alignment.cigar_vec, 0, 0, 1, 2);

				//		cerr << up_alignment.chr << '\t' << up_alignment.pos - microhomology_length << '\t' << '+' << '\t'
				//	   << down_alignment.chr << '\t' << down_alignment.pos << '\t' << '+' << '\t' << microhomology_length << '\t' << read_id << endl;
					}
					else
					{
						microhomology_length = 0;
						junction.set_value(up_alignment.chr, up_alignment.pos, '+', down_alignment.chr, down_alignment.pos, '+');
						up_seq_info.set_value(down_alignment.left_seq, up_alignment.cigar_vec, 0, down_alignment.left_seq.length() - up_alignment.left_seq.length(), 0, 2);
						down_seq_info.set_value(down_alignment.right_seq, down_alignment.cigar_vec, 0, 0, 1, 2);
				//		cerr << up_alignment.chr << '\t' << up_alignment.pos << '\t' << '+' << '\t'
				//		 << down_alignment.chr << '\t' << down_alignment.pos << '\t' << '+' << '\t' << microhomology_length << '\t' << read_id << endl;
					}

				}
				else if (read_id2align_it->second.strand != strand && read_id2align_it->second.clipped_side == clipped_side)
				{
					if (make_pair(read_id2align_it->second.chr, read_id2align_it->second.pos) < make_pair(chr, pos))
					{
						up_alignment = read_id2align_it->second;
						down_alignment = alignment;
					}
					else
					{
						up_alignment = alignment;
						down_alignment = read_id2align_it->second;
					}

					if (clipped_side == '5')
					{
						if (up_alignment.right_seq.length() >= down_alignment.left_seq.length())
						{
							microhomology_length = up_alignment.right_seq.length() - down_alignment.left_seq.length();
							junction.set_value(up_alignment.chr, up_alignment.pos, '-', down_alignment.chr, down_alignment.pos + microhomology_length, '+');
							GetReverseComplementSeq(up_alignment.left_seq);
							GetReverseComplementSeq(up_alignment.right_seq);
							cigar_vec = down_alignment.cigar_vec;
							AddCigarLeft(cigar_vec, microhomology_length);
							up_seq_info.set_value(up_alignment.right_seq, up_alignment.cigar_vec, 0, 0, 0, 2);
							down_seq_info.set_value(up_alignment.left_seq, cigar_vec, 0, 0, 1, 2);
							//cerr << up_alignment.chr << '\t' << up_alignment.pos << '\t' << '-' << '\t'
						 	//<< down_alignment.chr << '\t' << down_alignment.pos + microhomology_length << '\t' << '+' << '\t' << microhomology_length << '\t' << read_id << endl;
						}
						else
						{
							microhomology_length = 0;
							junction.set_value(up_alignment.chr, up_alignment.pos, '-', down_alignment.chr, down_alignment.pos, '+');
							up_seq_info.set_value(down_alignment.left_seq, up_alignment.cigar_vec, 0, down_alignment.left_seq.length() - up_alignment.right_seq.length(), 0, 2);
							down_seq_info.set_value(down_alignment.right_seq, down_alignment.cigar_vec, 0, 0, 1, 2);
							//cerr << up_alignment.chr << '\t' << up_alignment.pos << '\t' << '-' << '\t'
						 	//<< down_alignment.chr << '\t' << down_alignment.pos << '\t' << '+' << '\t' << microhomology_length << '\t' << read_id << endl;
						}
					}
					else
					{
						if (up_alignment.left_seq.length() >= down_alignment.right_seq.length())
						{
							microhomology_length = up_alignment.left_seq.length() - down_alignment.right_seq.length();
							junction.set_value(up_alignment.chr, up_alignment.pos - microhomology_length, '+', down_alignment.chr, down_alignment.pos, '-');
							GetReverseComplementSeq(down_alignment.left_seq);
							GetReverseComplementSeq(down_alignment.right_seq);
							cigar_vec = up_alignment.cigar_vec;
							MinusCigarRight(cigar_vec, microhomology_length);
							up_seq_info.set_value(down_alignment.right_seq, cigar_vec, 0, 0, 0, 2);
							down_seq_info.set_value(down_alignment.left_seq, down_alignment.cigar_vec, 0, 0, 1, 2);
							//cerr << up_alignment.chr << '\t' << up_alignment.pos - microhomology_length << '\t' << '+' << '\t'
						 //	<< down_alignment.chr << '\t' << down_alignment.pos << '\t' << '-' << '\t' << microhomology_length << '\t' << read_id << endl;
						}
						else
						{
							microhomology_length = 0;
							junction.set_value(up_alignment.chr, up_alignment.pos, '+', down_alignment.chr, down_alignment.pos, '-');
							up_seq_info.set_value(up_alignment.left_seq, up_alignment.cigar_vec, 0, 0, 0, 2);
							down_seq_info.set_value(up_alignment.right_seq, down_alignment.cigar_vec, down_alignment.right_seq.length() - up_alignment.left_seq.length(), 0, 1, 2);

						//	cerr << up_alignment.chr << '\t' << up_alignment.pos << '\t' << '+' << '\t'
						 //	<< down_alignment.chr << '\t' << down_alignment.pos << '\t' << '-' << '\t' << microhomology_length << '\t' << read_id << endl;
						}
					}

				}
				multimap<Junction, OtherInfo>::iterator junction2other_it = junction2other.find(junction);
				if (junction2other_it == junction2other.end())
				{
					OtherInfo other_info(up_seq_info, down_seq_info, microhomology_length, 0);
					junction2other.insert(make_pair(junction, other_info));
				}
				else
				{
					int up_seq_len1 = junction2other_it->second.up_seq_info.seq.length();
					int down_seq_len1 = junction2other_it->second.down_seq_info.seq.length();
					int up_seq_len2 = up_seq_info.seq.length();
					int down_seq_len2 = down_seq_info.seq.length();

					if (up_seq_len1 != up_seq_len2 || down_seq_len1 != down_seq_len2)
					{
						junction2other_it->second.down_seq_info.support_read_no++;
					}
				}
				read_id2align.erase(read_id);
			}

		}
		else
		{
			read_id2align.insert(make_pair(read_id, alignment));
		}
	}
}



//int main(int argc, char *argv[])
//{
//	string file = argv[1];
//	multimap<Junction, OtherInfo> junction2other;
//	FindJunction(file, 1, junction2other);
//	multimap<Junction, OtherInfo>::iterator junction2other_it = junction2other.begin();
//	while (junction2other_it != junction2other.end())
//	{
//		cout << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t'
//			 << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << '\t'
//			 << junction2other_it->second.down_seq_info.support_read_no << '\t' << junction2other_it->second.microhomology_length << '\t';
//		for (int i = 0; i < junction2other_it->second.up_seq_info.cigar_vec.size(); i++)
//			cout << junction2other_it->second.up_seq_info.cigar_vec[i].first << junction2other_it->second.up_seq_info.cigar_vec[i].second;
//		cout << '\t';
//		for (int i = 0; i < junction2other_it->second.down_seq_info.cigar_vec.size(); i++)
//			cout << junction2other_it->second.down_seq_info.cigar_vec[i].first << junction2other_it->second.down_seq_info.cigar_vec[i].second;
//		cout << '\t' << junction2other_it->second.up_seq_info.seq << '\t' << junction2other_it->second.down_seq_info.seq << endl;
//		++junction2other_it;
//	}
//	return 0;
//}
