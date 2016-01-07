/*
 * AUTHOR: Qiu Kunlong
 * EMAIL:  qiukunlong@genomics.cn
 * DATE:   2012/02/11
 *
 */
#include "somatic.h"

typedef pair<multimap<pair<string, int>, Breakpoint>::iterator, multimap<pair<string, int>, Breakpoint>::iterator> iterator_pair;



// read breakpoint information of tumor sample
void ReadTumorFileAndOutputSomaticInfo(string normal_bam_file, string file, string output_file, multimap<pair<string, int>, ReadsInfo> &aligned_pos2reads_3clip, multimap<pair<string, int>, ReadsInfo> &aligned_pos2reads_5clip, int min_mapQ, int offset, double min_map_rate, int mean_insert_size, int deviation, int times)
{
	ofstream fout(output_file.c_str());
	if (!fout) { cerr << "Error: Cannot open output file " << output_file << endl; exit(1); }
	ifstream fin(file.c_str());
	if (!fin) { cerr << "Error: Cannot open output file " << file << endl; exit(1); }
	string up_chr, down_chr, sv_type, up_cigar, down_cigar, up_seq, down_seq, temp;
	int up_pos, down_pos, up_reads_no, down_reads_no, microhomology_length, abnormal_read_pair_no, up_depth, down_depth, up_up_depth, up_down_depth, down_up_depth, down_down_depth;
	double up_clip_rate, down_clip_rate;
	char up_strand, down_strand;
	
	samfile_t *samfin;
	char in_mode[5], *fn_list = 0;
	in_mode[0] = 'r';
	int is_bamin = 0;
	if (normal_bam_file.rfind(".bam") == normal_bam_file.size() - 4)
	{
		//if in_mode[1] == 'b', it will read a bam file
		in_mode[1] = 'b';
		is_bamin = 1;
	}

    if ((samfin = samopen(normal_bam_file.c_str(), in_mode, fn_list)) == 0)
	{
		cerr << "[main_samview] fail to open file for reading." << endl;
   		exit(1);
	}
	if (samfin->header == 0) 
	{
		cerr << "[main_samview] fail to read the header." << endl;
		exit(1);
	}

	bam_index_t *idx = 0;
	if (is_bamin) idx = bam_index_load(normal_bam_file.c_str()); // load BAM index
	if (idx == 0)
	{ // index is unavailable
		cerr << "[main_samview] random alignment retrieval only works for indexed BAM files.\n";
		cerr << normal_bam_file << endl;
		exit(1);
	}
	pair<multimap<pair<string, int>, ReadsInfo>::iterator, multimap<pair<string, int>, ReadsInfo>::iterator> mmap_it_pair;


	while (fin >> up_chr)
	{
		if (up_chr[0] == '@') 
		{
			getline(fin,temp);
			fout << up_chr << temp << '\t' << "left_clip_read_NO_of_control\tright_clip_read_NO_of_control\tabnormal_read_pair_no_of_control"<< endl;
			continue;
		}
		fin >> up_pos >> up_strand >> up_reads_no >> down_chr >> down_pos >> down_strand >> down_reads_no >> microhomology_length >> abnormal_read_pair_no >> sv_type >> up_depth >> down_depth >> up_up_depth >> up_down_depth >> down_up_depth >> down_down_depth >> up_clip_rate >> down_clip_rate >> up_cigar >> down_cigar >> up_seq >> down_seq;
		getline(fin, temp);


		int normal_left_reads = 0, normal_right_reads = 0; // Added in 20121108
		Junction junction(up_chr, up_pos, up_strand, down_chr, down_pos, down_strand);

		if ('+' == up_strand && '+' == down_strand)
		{
			if (microhomology_length != -1)
			{
				int down_smallest_pos = down_pos, up_largest_pos = up_pos + microhomology_length;
				mmap_it_pair = aligned_pos2reads_5clip.equal_range(make_pair(down_chr, down_smallest_pos));


				while (mmap_it_pair.first != mmap_it_pair.second)
				{
					if (CompareStringBeginFirst(down_seq, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq, mmap_it_pair.first->second.seq_left) >= min_map_rate)
					{
						normal_right_reads = mmap_it_pair.first->second.support_read_no;
						break;
					}
					++mmap_it_pair.first;
				}

				mmap_it_pair = aligned_pos2reads_3clip.equal_range(make_pair(up_chr, up_largest_pos));
				string up_seq1 = up_seq;
				if (down_seq.length() >= microhomology_length)
				{
					up_seq1.append(down_seq, 0, microhomology_length);

					string down_seq1 = down_seq.substr(microhomology_length);

					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						if (CompareStringBeginFirst(down_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_left_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}
				}
				int normal_abnormal_read_pair_no = 0;

				normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
				fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
			}
			else
			{
				if (up_reads_no == 0)
				{
					mmap_it_pair = aligned_pos2reads_5clip.equal_range(make_pair(down_chr, down_pos));
					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						if (CompareStringBeginFirst(down_seq, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_right_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}


					multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = aligned_pos2reads_3clip.lower_bound(make_pair(up_chr, up_pos));
					while (mmap_it != aligned_pos2reads_3clip.end() && mmap_it->first.first == up_chr && mmap_it->first.second <= up_pos + offset)
					{
						//seq2 is 3'clipped  seq4 is 3' aligned
						if (Compare(mmap_it->second.seq_left, mmap_it->second.seq_right, up_seq, down_seq, min_map_rate) != -1)
						{
							normal_left_reads = mmap_it->second.support_read_no;
							break;
						}
						++mmap_it;
					}
					int normal_abnormal_read_pair_no = 0;
					normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
					fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
				}
				else if (down_reads_no == 0)
				{
					mmap_it_pair = aligned_pos2reads_3clip.equal_range(make_pair(up_chr, up_pos));
					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						if (CompareStringBeginFirst(down_seq, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_left_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}
					multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = aligned_pos2reads_5clip.lower_bound(make_pair(down_chr, down_pos - offset));
					while (mmap_it != aligned_pos2reads_5clip.end() && mmap_it->first.first == down_chr && mmap_it->first.second <= down_pos)
					{
						//seq2 is 3'clipped  seq4 is 3' aligned
						if (Compare(up_seq, down_seq, mmap_it->second.seq_left, mmap_it->second.seq_right, min_map_rate) != -1)
						{
							normal_right_reads = mmap_it->second.support_read_no;
							break;
						}
						++mmap_it;
					}
					int normal_abnormal_read_pair_no = 0;
					normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
					fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
				}
				else
				{
					cerr << "The tandem repeat length is error in postion: " << up_chr << '\t' << up_pos << '\t' << down_chr << '\t' << down_pos << endl;
				}
			}

		}
		else if ('+' == up_strand && '-' == down_strand)
		{
			if (microhomology_length != -1)
			{
				int up_largest_pos = up_pos + microhomology_length;
				//Here, both up_chr and down_chr are 3 end clipped.
				mmap_it_pair = aligned_pos2reads_3clip.equal_range(make_pair(up_chr, up_largest_pos));
				string up_seq1 = up_seq;
				up_seq1.append(down_seq, 0, microhomology_length);
				string down_seq1 = down_seq.substr(microhomology_length);

				while (mmap_it_pair.first != mmap_it_pair.second)
				{
					if (CompareStringBeginFirst(down_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
					{
						normal_left_reads = mmap_it_pair.first->second.support_read_no;
						break;
					}
					++mmap_it_pair.first;
				}

				mmap_it_pair = aligned_pos2reads_3clip.equal_range(make_pair(down_chr, down_pos));
				up_seq1 = up_seq, down_seq1 = down_seq;
				GetReverseComplementSeq(up_seq1);
				GetReverseComplementSeq(down_seq1);

				while (mmap_it_pair.first != mmap_it_pair.second)
				{
					if (CompareStringBeginFirst(up_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(down_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
					{
						normal_right_reads = mmap_it_pair.first->second.support_read_no;
						break;
					}
					++mmap_it_pair.first;
				}
				int normal_abnormal_read_pair_no = 0;
				normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
				fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
			}
			else
			{
				if (up_reads_no == 0)
				{
					mmap_it_pair = aligned_pos2reads_3clip.equal_range(make_pair(down_chr, down_pos));
					string up_seq1 = up_seq, down_seq1 = down_seq;
					GetReverseComplementSeq(up_seq1);
					GetReverseComplementSeq(down_seq1);

					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						//because up_seq1 and down_seq1 have been reversed, so up_seq1 is compared from begin end.
						if (CompareStringBeginFirst(up_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(down_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_right_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}

					multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = aligned_pos2reads_3clip.lower_bound(make_pair(up_chr, up_pos));
					while (mmap_it != aligned_pos2reads_3clip.end() && mmap_it->first.first == up_chr && mmap_it->first.second <= up_pos + offset)
					{
						//seq2 is 3'clipped  seq4 is 3' aligned
						if (Compare(mmap_it->second.seq_left, mmap_it->second.seq_right, up_seq, down_seq, min_map_rate) != -1)
						{
							normal_left_reads = mmap_it->second.support_read_no;
							break;
						}
						++mmap_it;
					}
					int normal_abnormal_read_pair_no = 0;
					normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
					fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
				}
				else if (down_reads_no == 0)
				{
					mmap_it_pair = aligned_pos2reads_3clip.equal_range(make_pair(up_chr, up_pos));
					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						if (CompareStringBeginFirst(down_seq, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_left_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}


					multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = aligned_pos2reads_3clip.lower_bound(make_pair(down_chr, down_pos));
					string up_seq1 = up_seq, down_seq1 = down_seq;
					GetReverseComplementSeq(up_seq1);
					GetReverseComplementSeq(down_seq1);

					while (mmap_it != aligned_pos2reads_3clip.end() && mmap_it->first.first == down_chr && mmap_it->first.second <= down_pos + offset)
					{
						//seq2 is 3'clipped  seq4 is 3' aligned
						if (Compare(mmap_it->second.seq_left, mmap_it->second.seq_right, down_seq1, up_seq1, min_map_rate) != -1)
						{
							normal_right_reads = mmap_it->second.support_read_no;
							break;
						}
						++mmap_it;
					}
					int normal_abnormal_read_pair_no = 0;
					normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
					fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
				}
				else
				{
					cerr << "The tandem repeat length is error in postion: " << up_chr << '\t' << up_pos << '\t' << down_chr << '\t' << down_pos << endl;
				}
			}
		}
		else if ('-' == up_strand && '+' == down_strand)
		{
			if (microhomology_length != -1)
			{
				string up_seq1 = up_seq, down_seq1 = down_seq;
				GetReverseComplementSeq(up_seq1);
				GetReverseComplementSeq(down_seq1);
				//Here, both up_chr and down_chr are 5 end clipped.
				mmap_it_pair = aligned_pos2reads_5clip.equal_range(make_pair(up_chr, up_pos));
				while (mmap_it_pair.first != mmap_it_pair.second)
				{
					if (CompareStringBeginFirst(up_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(down_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
					{
						normal_left_reads = mmap_it_pair.first->second.support_read_no;
						break;
					}
					++mmap_it_pair.first;
				}
				

				int down_smallest_pos = down_pos - microhomology_length;
				mmap_it_pair = aligned_pos2reads_5clip.equal_range(make_pair(down_chr, down_smallest_pos));
				up_seq1 = up_seq.substr(0, up_seq.length() - microhomology_length);
				down_seq1 = up_seq.substr(up_seq.length() - microhomology_length);
				down_seq1.append(down_seq);
				while (mmap_it_pair.first != mmap_it_pair.second)
				{
					if (CompareStringBeginFirst(down_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
					{
						normal_right_reads = mmap_it_pair.first->second.support_read_no;
						break;
					}
					++mmap_it_pair.first;
				}
				int normal_abnormal_read_pair_no = 0;
				normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
				fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;

			}
			else
			{
				if (up_reads_no == 0)
				{
					mmap_it_pair = aligned_pos2reads_5clip.equal_range(make_pair(down_chr, down_pos));
					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						if (CompareStringBeginFirst(down_seq, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(up_seq, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_right_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}


					multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = aligned_pos2reads_5clip.lower_bound(make_pair(up_chr, up_pos - offset));
					string up_seq1 = up_seq, down_seq1 = down_seq;
					GetReverseComplementSeq(up_seq1);
					GetReverseComplementSeq(down_seq1);
					
					while (mmap_it != aligned_pos2reads_5clip.end() && mmap_it->first.first == up_chr && mmap_it->first.second <= up_pos)
					{
						if (Compare(up_seq1, down_seq1, mmap_it->second.seq_left, mmap_it->second.seq_right, min_map_rate) != -1)
						{
							normal_left_reads = mmap_it->second.support_read_no;
							break;
						}
						++mmap_it;
					}
					int normal_abnormal_read_pair_no = 0;
					normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
					fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
				}
				else if (down_reads_no == 0)
				{
					mmap_it_pair = aligned_pos2reads_5clip.equal_range(make_pair(up_chr, up_pos));
					string up_seq1 = up_seq, down_seq1 = down_seq;
					GetReverseComplementSeq(up_seq1);
					GetReverseComplementSeq(down_seq1);

					while (mmap_it_pair.first != mmap_it_pair.second)
					{
						if (CompareStringBeginFirst(up_seq1, mmap_it_pair.first->second.seq_right) >= min_map_rate && CompareStringEndFirst(down_seq1, mmap_it_pair.first->second.seq_left) >= min_map_rate)
						{
							normal_left_reads = mmap_it_pair.first->second.support_read_no;
							break;
						}
						++mmap_it_pair.first;
					}


					multimap<pair<string, int>, ReadsInfo>::iterator mmap_it = aligned_pos2reads_5clip.lower_bound(make_pair(down_chr, down_pos - offset));
					while (mmap_it != aligned_pos2reads_5clip.end() && mmap_it->first.first == down_chr && mmap_it->first.second <= down_pos)
					{
						if (Compare(up_seq, down_seq, mmap_it->second.seq_left, mmap_it->second.seq_right, min_map_rate) != -1)
						{
							normal_right_reads = mmap_it->second.support_read_no;
							break;
						}
						++mmap_it;
					}
					int normal_abnormal_read_pair_no = 0;
					normal_abnormal_read_pair_no = FindDiscordantReadPairs(samfin, idx, junction, min_mapQ, mean_insert_size, deviation, times);
					fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << abnormal_read_pair_no << '\t' << sv_type << '\t' << up_depth << '\t' << down_depth << '\t' << up_up_depth << '\t' << up_down_depth << '\t' << down_up_depth << '\t' << down_down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << '\t' << normal_left_reads << '\t' << normal_right_reads << '\t' <<  normal_abnormal_read_pair_no << endl;
				}
				else
				{
					cerr << "The tandem repeat length is error in postion: " << up_chr << '\t' << up_pos << '\t' << down_chr << '\t' << down_pos << endl;
				}
			}
			
		}
		
		else
		{
			cerr << "Error: Something error in line " << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << up_reads_no << '\t' << down_chr << '\t' << down_pos << '\t' << down_strand << '\t' << down_reads_no << '\t' << microhomology_length << '\t' << up_depth << '\t' << down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << up_cigar << '\t' << down_cigar << '\t' << up_seq << '\t' << down_seq << endl;
		}

	}
}

/*
void CalculateSomaticBreakpoint(ofstream &fout, multimap<pair<string, int>, Breakpoint> &uppos2breakpoint, string up_chr, int up_pos, char up_strand, char down_strand, int up_depth, int down_depth, double up_clip_rate, double down_clip_rate, Breakpoint &breakpoint, int offset)
{
	multimap<pair<string, int>, Breakpoint>::iterator mmap_it;
	iterator_pair it_pair;
	if (breakpoint.microhomology_length == 0)
	{	
		it_pair = uppos2breakpoint.equal_range(make_pair(up_chr, up_pos));
		while (it_pair.first != it_pair.second)
		{
			if (it_pair.first->second.down_chr == breakpoint.down_chr && it_pair.first->second.down_pos == breakpoint.down_pos)
				break;
			else
				++it_pair.first;
		}
		if (it_pair.first == it_pair.second)
		{
			fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << breakpoint.up_seq_info.support_read_no << '\t' << breakpoint.down_chr << '\t' << breakpoint.down_pos << '\t' << down_strand << '\t' << breakpoint.down_seq_info.support_read_no << '\t' << breakpoint.microhomology_length << '\t' << up_depth << '\t' << down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << breakpoint.up_seq_info.cigar << '\t' << breakpoint.down_seq_info.cigar << '\t' << breakpoint.up_seq_info.seq << '\t' << breakpoint.down_seq_info.seq << endl;
		}
	}
	else if (breakpoint.microhomology_length == -1 || breakpoint.microhomology_length > 0)
	{
		int smallest_up_pos, largest_up_pos;
		if (breakpoint.microhomology_length == -1)
		{
			if (breakpoint.up_seq_info.support_read_no == 0)
			{
				smallest_up_pos = up_pos;
				largest_up_pos = up_pos + offset;
			}
			else
			{
				smallest_up_pos = up_pos - offset;
				largest_up_pos = up_pos;
			}
		}
		else 
		{
			smallest_up_pos = up_pos;
			largest_up_pos = up_pos + breakpoint.microhomology_length;
		}
		mmap_it = uppos2breakpoint.lower_bound(make_pair(up_chr, smallest_up_pos));
		while (mmap_it != uppos2breakpoint.end() && mmap_it->first.first == up_chr && mmap_it->first.second <= largest_up_pos)
		{
			if ('+' == up_strand && '+' == down_strand)
			{
				if (mmap_it->second.down_chr == breakpoint.down_chr && mmap_it->second.down_pos - mmap_it->first.second == breakpoint.down_pos - up_pos)
				{
					break;
				}
				else
				{
					++mmap_it;
				}
			}
			else
			{
				if (mmap_it->second.down_chr == breakpoint.down_chr && mmap_it->second.down_pos + mmap_it->first.second == breakpoint.down_pos + up_pos)
					break;
				else
					++mmap_it;
			}
		}
		if (mmap_it == uppos2breakpoint.end() || mmap_it->first.first != up_chr || mmap_it->first.second > largest_up_pos)
		{
			fout << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << breakpoint.up_seq_info.support_read_no << '\t' << breakpoint.down_chr << '\t' << breakpoint.down_pos << '\t' << down_strand << '\t' << breakpoint.down_seq_info.support_read_no << '\t' << breakpoint.microhomology_length << '\t' << up_depth << '\t' << down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << breakpoint.up_seq_info.cigar << '\t' << breakpoint.down_seq_info.cigar << '\t' << breakpoint.up_seq_info.seq << '\t' << breakpoint.down_seq_info.seq << endl;
		}

	}
	else
	{
		cerr << "Error: Something error in line " << up_chr << '\t' << up_pos << '\t' << up_strand << '\t' << breakpoint.up_seq_info.support_read_no << '\t' << breakpoint.down_chr << '\t' << breakpoint.down_pos << '\t' << down_strand << '\t' << breakpoint.down_seq_info.support_read_no << '\t' << breakpoint.microhomology_length << '\t' << up_depth << '\t' << down_depth << '\t' << up_clip_rate << '\t' << down_clip_rate << '\t' << breakpoint.up_seq_info.cigar << '\t' << breakpoint.down_seq_info.cigar << '\t' << breakpoint.up_seq_info.seq << '\t' << breakpoint.down_seq_info.seq << endl;
	}
}    	
*/
