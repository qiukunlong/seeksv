/*
 ***************************************
 * cluster.cpp                         *
 *                                     *
 *  Created on:   2013-1-7             *
 *  Updated on:   2016-7-7             *
 *  Author:       Qiu Kunlong          *
 *                                     *
 ***************************************
 */
#include "cluster.h"
#include "head.h"
#include "clip_reads.h"

bool CalculateInsertsizeDeviation(string file, int min_mapQ, int read_pair_used, int &mean_insert_size, int &deviation)
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
		cerr << "[main_samview] fail to open file " << file << "for reading." << endl;
   		exit(1);
	}
	if (samfin->header == 0) 
	{
		cerr << "[main_samview] fail to read the header of file " << file  << endl;
		exit(1);
	}

	b = bam_init1();
	int ret = 0;

	g_min_mapQ = min_mapQ;

	unsigned long total_insert_size = 0;
	int insert_size = 0, read_pair_number = 0;
	vector<int> vec_insert_sizes;

	//read a line of the bam file, if the line is without soft-clipping read, ignore it
	while ((ret = samread(samfin, b) >= 0))
	{
		// when read failed continue
		if (__g_skip_aln(samfin->header, b)) continue;

		insert_size = b->core.isize;
//		if (!(b->core.flag&BAM_FREVERSE) && (b->core.flag&BAM_FMREVERSE) && !(b->core.flag&BAM_FDUP) && insert_size > 0 && insert_size <= maximum_insert_size)
//		{
//			total_insert_size += insert_size;
//			vec_insert_sizes.push_back(insert_size);
//			++read_pair_number;
//		}
//      ignore hard clipped reads
		if (IsHardClip(b)) continue;
		if ((b->core.flag&BAM_FPAIRED) && (b->core.flag&BAM_FPROPER_PAIR) && !(b->core.flag&BAM_FDUP) && insert_size > 0)
		{
			total_insert_size += insert_size;
			vec_insert_sizes.push_back(insert_size);
			++read_pair_number;
		}
		if (read_pair_number == read_pair_used) break;
	}

	if (read_pair_number == 0) return 1;
	mean_insert_size = total_insert_size / read_pair_number;
	double d_deviation = 0;
	vector<int>::iterator vec_it = vec_insert_sizes.begin();
	while (vec_it != vec_insert_sizes.end())
	{
		d_deviation += ((*vec_it) - mean_insert_size) * ((*vec_it) - mean_insert_size);
		++ vec_it;
	}
	deviation = (int)sqrt(d_deviation/read_pair_number);
	cerr << "Bam/sam " << file << "    Mean insert size : " << mean_insert_size << "\n" << "Mean deviation: " << deviation << endl;
	return 0;
}


bool DivideCluster(string file, int min_mapQ, int &mean_insert_size, int &deviation)
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

	unsigned long total_insert_size = 0;
	int insert_size = 0, read_pair_number = 0;
	vector<int> vec_insert_sizes;
	char *buffer;

	//read a line of the bam file, if the line is without soft-clipping read, ignore it
	while ((ret = samread(samfin, b) >= 0))
	{
		// when read failed continue
		if (__g_skip_aln(samfin->header, b)) continue;

		insert_size = b->core.isize;
		if (IsConcordant(samfin->header, b, mean_insert_size, deviation, 5) || (b->core.flag&BAM_FDUP) || (b->core.flag&BAM_FREAD2)) continue;
		buffer = bam_format1(samfin->header, b); // read the buffer, and store it into buffer
		cout << buffer << endl;
		free(buffer);
	}

	return 0;
}


bool IsConcordant(const bam_header_t *header, const bam1_t *b, int mean_insert_size, int deviation, int times)
{
	int min_insert_size = mean_insert_size - deviation * times, max_insert_size = mean_insert_size + deviation * times, insert_size = b->core.isize;
	if (!(b->core.flag&BAM_FREVERSE) && (b->core.flag&BAM_FMREVERSE) && min_insert_size <= insert_size && insert_size <= max_insert_size) { return 1; }
	else if ((b->core.flag&BAM_FREVERSE) && !(b->core.flag&BAM_FMREVERSE) && insert_size < 0)
	{
		insert_size = abs(insert_size);
		if (min_insert_size <= insert_size && insert_size <= max_insert_size) return 1;
		else return 0;
	}
	else { return 0; }
}

bool RandomFetch(string file, int min_mapQ, int &mean_insert_size, int &deviation)
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

	int ret = 0;
	bam_index_t *idx = 0;
	idx = bam_index_load(file.c_str()); // load BAM index
	if (idx == 0)
	{ // index is unavailable
		cerr << "[main_samview] random alignment retrieval only works for indexed BAM files.\n" << endl;
		exit(1);
	}


	g_min_mapQ = min_mapQ;

	map<string, int> seq_name2tid;
	StoreSeqName2Tid(samfin->header, seq_name2tid);
	int tid = BamGetTid(seq_name2tid, "chr17");
	//beg is 0-base and end is 1-base
	int beg = 199999, end = 20000;
	//bam_parse_region(samfin->header, "chr17:20000-20000", &tid, &beg, &end);
	//cout << tid << '\t' << beg << '\t' << end << endl;

	bam_iter_t iter = bam_iter_query(idx, tid, beg, end);
	bamFile fp = samfin->x.bam;
	b = bam_init1();
	g_min_mapQ = 60;
	while ((ret = bam_iter_read(fp, iter, b)) >= 0)
	{
		if (__g_skip_aln(samfin->header, b)) continue;
		cout << b->core.tid << '\t' << b->core.pos << endl;
	}


	return 0;
}

int StoreSeqName2Tid(const bam_header_t *header, map<string, int> &seq_name2tid)
{
	int i = 0;	
	string seq_name;
	for (i = 0; i < header->n_targets; i++)
	{
		seq_name = header->target_name[i];
		seq_name2tid.insert(make_pair(seq_name, i));
	}
	return 0;
}

int BamGetTid(map<string, int> &seq_name2tid, const string seq_name)
{
	map<string, int>::iterator map_it = seq_name2tid.find(seq_name);
	if (map_it == seq_name2tid.end())
	{
		return -1;
	}
	else
	{
		return map_it->second;
	}
}
