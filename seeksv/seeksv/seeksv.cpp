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

const char *kVersion = "1.2.2-debug";
//const char *kRevision = "43";
const int kCommandQuantity = 4;
const double kThreshold = 0.85;
const int kNumberOfReadPair = 5000000;

void Usage(char *prog);
void Usage(char *prog, char *command, int i);
char *comarray[kCommandQuantity] = {"getclip", "getsv", "somatic", "cluster"};
void CallGetclip(int argc, char *argv[], int i);
void CallGetsv(int argc, char *argv[], int i);
void CallSomatic(int argc, char *argv[], int i);
void SelectStep(int argc, char *argv[], int i);

int main(int argc, char *argv[])
{
	if (argc == 1)
	{
		//Usage(argv[0]);
		Usage("seeksv");
		//cerr << "Usage: " << argv[0] << " <bamfile> <soft out> <fq out>" << endl;
	
	}
	else
	{
		int i = 0;
		for (; i < kCommandQuantity; i++)
		{
			if (strcmp(comarray[i], argv[1]) == 0)
				break;
		}
		if (i == kCommandQuantity)
		{
			cerr << "[seeksv] unrecognized command '" << argv[1] << "'" << endl;
			exit(1);
		}
		else if (argc == 2)
		{
			//Usage(argv[0], argv[1], i);
			Usage("seeksv", argv[1], i);
		}
		else
		{
			SelectStep(argc, argv, i);
		}
	}
}

void Usage(char *prog)
{
	cerr << "Program: seeksv (software used to dectect structural variations)" << '\n'
		 //<< "Compile Date: " << system("date") << '\n'
		 << "Version: " << kVersion << '\n'
		 << "Contract: " << "Kunlong Qiu(qkl871118@qq.com)\n\n"
		 << "Usage: " << prog << " <command> [option]\n\n"
		 << "Command: " << "getclip\tget soft-clipped reads\n"
		 << "         " << "getsv  \tget final sv\n"
		 << "         " << "somatic\tget somatic sv" << endl;
//		 << "         " << "cluster\tget cluster" << endl;
	exit(1);
}


void Usage(char *prog, char *command, int i)
{
	switch (i)
	{
	case 0:
		//cerr << "Usage: " << prog << " " << command << " [options] <in.sorted.bamfile> <out.soft clipped reads file> <out clipped sequence fastq file>\n\n";
		cerr << "Usage: " << prog << " " << command << " [options] <in.sorted.bamfile>\n\n";
		cerr << "Options: -t <double>           Threshold of match rate while combining two soft-clipped reads [0.85]" << endl;
		cerr << "         -q <int>              Minimum mapping quality of soft-clipped reads [20]" << endl;
		cerr << "         -s                    Save the low quality sequence clipped before alignment by bwa." << endl;
		cerr << "         -o <string>           Prefix of output files [output]" << endl;
		break;
	case 1:
		cerr << "Usage: " << prog << " " << command << " [options] <bamfile of clipped sequence> <orignal sorted bamfile> <soft-clipped reads file> <breakpoint> <unmaped clipped sequence fastq result>\n"
			 << "Options: -F <FILE>             Samfile/Bamfile of connected readthrough reads\n" 
		//	 << "         -B <FILE>             temp breakpoint(sv) result\n"
			 << "         -t <double>           Threshold of match rate while combining two soft-clipped reads [0.85]\n" 
			 << "         -l <int>              Maximum search length to find microhomology[50]\n"
			 << "         -q <int>              Minimum mapping quality of discordant read pair [20]\n"
			 << "         -Q <int>              Minimum mapping quality of clipped sequences [1], if you use bwa samse to align the reads, please set\n"              << "                               this flag to 20\n"
			 << "         -w <int>              Minimum mapping quality connected  readthrough reads [1]\n"
			 << "         -n <int>              Number of segment(read pairs) used to calculate insert size default [5000000], if you donot want to use abnormal read pairs to call sv, set this parameter to 0.\n"
			 << "         -r                    Turn off rescue mode. When rescure mode is on, the a SV with only 1 side with enough soft-clipped \n"
			 << "                               reads is considered as a valid one instead of rejecting it.  Default on.\n"
			 << "         -a <int>              When rescure mode is on, minimum number of soft-clipped reads on one side [5]\n" 
			 << "         -b <int>              Minimum number of soft clipping read,it's the sum of left clipped reads and right clipped reads [3]\n"
			 << "         -d <int>              Minimum distance between the clipped sequence position and the aligned sequence position [50]\n"
			 << "         -D                    Do not calculate depth of the breakpoints and their ajacency regions\n"
			 << "         -e <int>              Minimum number of read pairs which support the junction [0]\n"
			 << "         -R <int>              Minimun number of coverage to be treated as repetitive [500]\n"
			 << "         -f <int>              Minimun mutation frequency(left_pos_clip_percentage >= 0.1 or right_pos_clip_percentage >= 0.1) [0.1].\n"
			 << "                               If you set -D, this value is invalid and set to [0].\n"
			 << "         -T <int>              Maximum length of microhomology, microhomology length longer than [50] will be filtered\n"
			 << "         -m <int>              Minimum length of up_seq or down_seq near the breakpoint when abnormal_read_pair_no == 0 [30]\n"
			 << "         -i <int>              Maximum indel number of up_seq or down_seq near the breakpoint when abnormal_read_pair_no == 0 [1]\n"
			 << "         -L <int>              Calculate average depth  [200] bp upstream or downstream of the breakpoints\n"
			 << endl;
		break;
		
	case 2:
		cerr << "Usage: " << prog << " " << command << " [options] <normal original bam file> <normal soft-clipped reads file> <tumor breakpoint file> <output somatic breakpoint file>\n" << endl;
		cerr << "         -t <int>              Threshold of match rate while comparing two soft-clipped reads [0.85]" << endl;
		cerr << "         -q <int>              Minimum mapping quality of discordant read pair [20]" << endl;
		cerr << "         -l <int>              Maximum search length to find microhomology [30]" << endl;
		cerr << "         -m <int>              Minimum length of the clipped sequence  in normal [10]" << endl;
		cerr << "         -n <int>              Number of segment(read pairs) used to calculate insert size default [5000000], if you donot want to use abnormal read pairs to call sv, set this parameter to 0." << endl;
		break;
//	case 3:
//		cerr << "Usage: " << prog << " " << command << " [options] <orignal sorted bamfile>\n" << endl;
//		cerr << "         -a <int>              maximum insert size [1000]" << endl;
//		cerr << "         -n <int>              number of segment(read pairs) used to calculate insert size [5000000]" << endl;
//		cerr << "         -q <int>              minimum mapping quality [20]" << endl;

	}
	exit(1);
}


void CallGetclip(int argc, char *argv[], int i)
{
	int c, min_mapQ = 20;
	double threshold = kThreshold;
	string prefix = "output";
	bool save_low_quality = 0;
	while ((c = getopt(argc, argv, "t:q:o:s")) >= 0)
	{
		switch (c)
		{
		case 't': threshold = atof(optarg); break;
		case 'q': min_mapQ = atoi(optarg); break;
		case 's': save_low_quality = 1; break;
		case 'o': prefix = optarg; break;
		}
	}

	if (argc != optind + 1)
		Usage(argv[0], argv[1], i);

	string bamfile = argv[optind++];
	
	multimap<pair<string, int>, ReadsInfo> breakpoint2read_l, breakpoint2read_r;
	map<string, pair<pair<string, string>, char> > id2seq_pair;

	InputBamOutputReads<ogzstream> (bamfile, breakpoint2read_l, breakpoint2read_r, id2seq_pair, threshold, min_mapQ, save_low_quality, prefix);

}

void CallGetsv(int argc, char *argv[], int i)
{
	string connect_bam, temp_breakpoint;
	double threshold = kThreshold, frequency = 0.1;
	int c, flank = 50, min_mapQ = 20, min_mapQ1 = 1, min_mapQ2 = 1, read_pair_used = 5000000, baseQ = 0, min_no_one_side_clipped_reads = 5, sum_min_no_both_clipped_reads = 3, min_distance = 50, microhomology_length = 50, times = 4, min_abnormal_read_pair_no = 0, repeat_coverage = 500, flank_length = 200, min_seq_len = 30, max_seq_indel_no = 1;
	bool rescure_mode = 1;
	bool output_depth = 1;
	while ((c = getopt(argc, argv, "F:B:t:l:q:Q:w:n:a:b:d:e:m:i:R:f:T:L:rD")) >= 0)
	{
		switch (c)
		{
		case 'F': connect_bam = optarg; break;
		case 'B': temp_breakpoint = optarg; break;
		case 't': threshold = atof(optarg); break;
		case 'l': flank = atoi(optarg); break;
		case 'q': min_mapQ = atoi(optarg); break;
		case 'Q': min_mapQ1 = atoi(optarg); break;
		case 'w': min_mapQ2 = atoi(optarg); break;
		case 'n': read_pair_used = atoi(optarg); break; 
		case 'r': rescure_mode = 0; break;
		case 'a': min_no_one_side_clipped_reads = atoi(optarg); break;
		case 'b': sum_min_no_both_clipped_reads = atoi(optarg); break;
		case 'd': min_distance = atoi(optarg); break;
		case 'e': min_abnormal_read_pair_no = atoi(optarg); break;
		case 'm': min_seq_len = atoi(optarg); break;
		case 'i': max_seq_indel_no = atoi(optarg); break;
		case 'D': output_depth = 0; break;
		case 'R': repeat_coverage = atoi(optarg); break;
		case 'f': frequency = atof(optarg); break;
		case 'T': microhomology_length = atoi(optarg); break;
		case 'L': flank_length = atoi(optarg); break;
		}
	}
	if (argc != optind + 5)
	{
		Usage(argv[0], argv[1], i);
		//cerr << "aa" << endl;
	}
	else if (flank > 90 || flank < 0 || min_seq_len < 0)
	{
		Usage(argv[0], argv[1], i);
		//cerr << "bb" << endl;
	}

	string file = argv[optind++];
	//char **ba = argv + optind;	
	string original_bam = argv[optind++];
	string clipfile = argv[optind++], breakpoint_file = argv[optind++], clip_unmap_fq_file = argv[optind++];
	multimap<string, AlignInfo> read_id2align_info;
	multimap<pair<string, int>, ClipReads> aligned2clipped;
	map<pair<string, int>, int> pos2depth;
	map<ChrRange, unsigned long> range2depth;
	map<pair<string, int>, int> begin2end;
	map<Junction, pair<pair<ChrRange, ChrRange>, pair<ChrRange, ChrRange> > > junction2range_pair;

	//StoreSamOfClippedSeq(file, read_id2align_info, min_mapQ1);

	//cerr << "'StoreSamOfClippedSeq' finished" << endl;


	multimap<Junction, OtherInfo> junction2other;
	
	if (!temp_breakpoint.empty())
	{
		ReadBreakpoint(temp_breakpoint, junction2other);
		cerr << "[ReadBreakpoint] finish" << endl;
	}

	if (!connect_bam.empty())
	{
		FindJunction(connect_bam, min_mapQ2, junction2other);
		cerr << "'FindJunction' finished" << endl;
	}



	if (clipfile.rfind(".gz") == clipfile.length() - 3)
		//InputSoftInfoStoreBreakpoint<igzstream> (clipfile, read_id2align_info, junction2other, aligned2clipped, threshold);
		InputSoftInfoStoreBreakpoint<igzstream> (clipfile, file, read_id2align_info, junction2other, aligned2clipped, threshold);
	else 
		InputSoftInfoStoreBreakpoint<ifstream> (clipfile, file, read_id2align_info, junction2other, aligned2clipped, threshold);
	cerr << "'InputSoftInfoStoreBreakpoint' finished" << endl;

	MergeJunction(junction2other, flank);
//	cerr << "'MergeJunction' finished" << endl;



//	ChangeSameOrientationJuntion(uppos2breakpoint_same_strand, aligned2clipped, flank, threshold);

//	ChangeDiffOrientationJuntion(uppos2breakpoint_diff_strand_5clip, uppos2breakpoint_diff_strand_3clip, aligned2clipped, flank, threshold);

	int mean_insert_size = 0, deviation = 0;
	
	if (read_pair_used >= 100000) {
		CalculateInsertsizeDeviation(original_bam, min_mapQ, read_pair_used, mean_insert_size, deviation);
		cerr << "'CalculateInsertsizeDeviation' finished" << endl;
		
		samfile_t *samfin;
		char in_mode[5], *fn_list = 0;
		in_mode[0] = 'r';
		int is_bamin = 0;
		if (original_bam.rfind(".bam") == original_bam.size() - 4)
		{
			//if in_mode[1] == 'b', it will read a bam file
			in_mode[1] = 'b';
			is_bamin = 1;
		}
	
	    if ((samfin = samopen(original_bam.c_str(), in_mode, fn_list)) == 0)
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
		if (is_bamin) idx = bam_index_load(original_bam.c_str()); // load BAM index
		if (idx == 0)
		{ // index is unavailable
			cerr << "[main_samview] random alignment retrieval only works for indexed BAM files.\n" << endl;
			cerr << original_bam << endl;
			exit(1);
		}
		FindDiscordantReadPairs(samfin, idx, junction2other, min_mapQ, mean_insert_size, deviation, times);
		bam_index_destroy(idx);
		cerr << "'FindDiscordantReadPairs' finished" << endl;
	}
	else {
		min_abnormal_read_pair_no = 0;
	}

	if (output_depth == 1) 
	{
		GetBreak(junction2other, pos2depth, range2depth, junction2range_pair, flank_length);
		GetBreak(aligned2clipped, pos2depth);
		MergeOverlap(range2depth, begin2end);
		cerr << "'MergeOverlap' finished" << endl;
	
		main_depth(pos2depth, range2depth, begin2end, original_bam, 0, min_mapQ, 0);
		cerr << "'main_depth' finished" << endl;
	}
	else
	{
		frequency = 0;
	}


	ofstream fout(breakpoint_file.c_str());

	if (!fout) { cerr << "Cannot open file " << breakpoint_file << endl; exit(1); }

	fout << "@left_chr\tleft_pos\tleft_strand\tleft_clip_read_NO\tright_chr\tright_pos\tright_strand\tright_clip_read_NO\tmicrohomology_length\tabnormal_readpair_NO\tsvtype\tleft_pos_depth\tright_pos_depth\taverage_depth_of_left_pos_5end\taverage_depth_of_left_pos_3end\taverage_depth_of_right_pos_5end\taverage_depth_of_right_pos_3end\tleft_pos_clip_percentage\tright_pos_clip_percentage\tleft_seq_cigar\tright_seq_cigar\tleft_seq\tright_seq" << endl;
	/*
	map<Junction, OtherInfo>::iterator junction2other_it = junction2other.begin();
	while (junction2other_it != junction2other.end())
	{
		fout << junction2other_it->first.up_chr << '\t' << junction2other_it->first.up_pos << '\t' << junction2other_it->first.up_strand << '\t' << junction2other_it->second.up_seq_info.support_read_no << '\t'
		     << junction2other_it->first.down_chr << '\t' << junction2other_it->first.down_pos << '\t' << junction2other_it->first.down_strand << '\t' << junction2other_it->second.down_seq_info.support_read_no << '\t' << junction2other_it->second.microhomology_length << '\t';
			 DisplayCigarVector<ostream> (fout, junction2other_it->second.up_seq_info.cigar_vec);
			 fout << '\t';
			 DisplayCigarVector<ostream> (fout, junction2other_it->second.down_seq_info.cigar_vec);
			 fout << '\t';
			 fout << junction2other_it->second.up_seq_info.left_clipped_seq_length << '\t' << junction2other_it->second.up_seq_info.right_clipped_seq_length << '\t' << junction2other_it->second.down_seq_info.left_clipped_seq_length << '\t' << junction2other_it->second.down_seq_info.right_clipped_seq_length << '\t' << junction2other_it->second.up_seq_info.seq << '\t' << junction2other_it->second.down_seq_info.seq << endl;
			 
		++junction2other_it;
	}
	*/
	OutputBreakpoint(fout, junction2other, pos2depth, range2depth, junction2range_pair, rescure_mode, min_no_one_side_clipped_reads, sum_min_no_both_clipped_reads, min_abnormal_read_pair_no, repeat_coverage, frequency, min_distance, microhomology_length, min_seq_len, max_seq_indel_no);

	ofstream foutuq(clip_unmap_fq_file.c_str());
	if (!foutuq) { cerr << "Cannot open file " << clip_unmap_fq_file << endl; exit(1); }

	OutputOneendUnmapBreakpoint(fout, foutuq, aligned2clipped, pos2depth);
	




/*	map<pair<pair<string, int>, pair<string, int> >, pair<int, int> > samestrand_breakpoint2support_count;
	map<pair<pair<string, int>, pair<string, int> >, pair<int, int> > diffstrand_breakpoint2support_count;
	multimap<pair<string, int>, BreakpointInfoAndSupportReadsCount> samestrand_breakpoint2breakpoint;
	multimap<pair<string, int>, BreakpointInfoAndSupportReadsCount> diffstrand_breakpoint2breakpoint;
	map<pair<string, pair<int, int> >, SupportCounts> inversion_area2support_count;
	multimap<pair<string, int>, InsertSeqInfo> insert_pos2insert_seq;
	map<pair<string, pair<int, int> >, SupportReadsCountInfo> deletion2support_count;


//	SamInfo sam_info(argv[2]);
	map<string, AlignInfo> read_id2align_info;
	map<string, string> seq2qual;
	multimap<pair<string, int>, SeqQualEnd> mapped_end2unmapped_end;
	ofstream fout(argv[5]);
	string samfile = argv[2];
	StoreContigSam(samfile, read_id2align_info, seq2qual);
	InputSoftInfoStoreBreakpoint(argv[3], read_id2align_info, samestrand_breakpoint2support_count, diffstrand_breakpoint2support_count, seq2qual, mapped_end2unmapped_end);
	ProcessOneEndMppedClippedreads(argv[6], mapped_end2unmapped_end);
	ProcessBreakpoint(argv[4], samestrand_breakpoint2support_count, diffstrand_breakpoint2support_count, samestrand_breakpoint2breakpoint);
	GetInsertion(samestrand_breakpoint2breakpoint, insert_pos2insert_seq, 90);
	DisplayInsertion(fout, insert_pos2insert_seq);

	GetDeletion(samestrand_breakpoint2breakpoint, deletion2support_count);
	DisplayDeletion(fout, deletion2support_count, 90);

	GetInversion(diffstrand_breakpoint2support_count, inversion_area2support_count);

	DisplayInversion(fout, inversion_area2support_count, 90);
	*/
}

void CallSomatic(int argc, char *argv[], int i)
{
	int offset = 30, c, min_len_of_clipped_seq = 10, read_pair_used = 5000000, min_mapQ = 20;
	double min_map_rate = 0.85;
	while ((c = getopt(argc, argv, "t:q:l:m:n:")) >= 0)
	{
		switch (c)
		{
		case 't': min_map_rate = atof(optarg); break;
		case 'q': min_mapQ = atoi(optarg); break;
		case 'l': offset = atoi(optarg); break;
		case 'm': min_len_of_clipped_seq = atoi(optarg); break;
		case 'n': read_pair_used = atoi(optarg); break; 
		}
	}
	if (argc != optind + 4)
	{
		cerr << argc << '\t' << optind << endl;
		Usage(argv[0], argv[1], i);
	}
	else if (offset >= 90 || offset < 0)
	{
		cerr << "Error: value of -l must in range [0, 90) " << endl;
		Usage(argv[0], argv[1], i);
	}
	string normal_bam_file = argv[optind++];
	string clipped_file = argv[optind++];
	string tumor_file = argv[optind++];
	string somatic_file = argv[optind++];

	multimap<pair<string, int>, ReadsInfo> aligned_pos2reads_3clip, aligned_pos2reads_5clip;
	if (clipped_file.rfind(".gz") == clipped_file.length() - 3)
	{
		ReadsClipReads<igzstream> (clipped_file, aligned_pos2reads_3clip, aligned_pos2reads_5clip, min_len_of_clipped_seq);
	}
	else
	{
		ReadsClipReads<ifstream> (clipped_file, aligned_pos2reads_3clip, aligned_pos2reads_5clip, min_len_of_clipped_seq);
	}
	int mean_insert_size = 0, deviation = 0;
	if (read_pair_used >= 100000) {
		CalculateInsertsizeDeviation(normal_bam_file, min_mapQ, read_pair_used, mean_insert_size, deviation);
	}
	ReadTumorFileAndOutputSomaticInfo(normal_bam_file, tumor_file, somatic_file, aligned_pos2reads_3clip, aligned_pos2reads_5clip, min_mapQ, offset, min_map_rate, mean_insert_size, deviation, 4);
}


//Added on 2012-1-7
//Author : Qiu Kunlong
void CallCluster(int argc, char *argv[], int i)
{
	int number_of_read_pair = kNumberOfReadPair;
	int c;
	int min_mapQ = 20, maximum_insert_size = 1000;

	while ((c = getopt(argc, argv, "a:n:q:")) >= 0)
	{
		switch (c)
		{
		case 'a': maximum_insert_size = atoi(optarg); break;
		case 'n': number_of_read_pair = atoi(optarg); break;
		case 'q': min_mapQ = atoi(optarg); break;
		}
	}

	if (argc != optind + 1)
		Usage(argv[0], argv[1], i);

	string bamfile = argv[optind++];
	
	//CalculateInsertDeviation(bamfile);
	int mean_insert_size = 0, deviation = 0;
	CalculateInsertsizeDeviation(bamfile, min_mapQ, number_of_read_pair, mean_insert_size, deviation);
//	DivideCluster(bamfile, min_mapQ, mean_insert_size, deviation);
	//RandomFetch(bamfile, min_mapQ, mean_insert_size, deviation);
	
}

void SelectStep(int argc, char *argv[], int i)
{
	switch (i)
	{
	case 0:
		CallGetclip(argc - 1, argv + 1, i); break;
	case 1:
		CallGetsv(argc - 1, argv + 1, i); break;
	case 2:
		CallSomatic(argc - 1, argv + 1, i); break;
//	case 3:
//		CallCluster(argc - 1, argv + 1, i); break;
	}
}
