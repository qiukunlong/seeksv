/*
 ***************************************
 * cluster.h                           *
 *                                     *
 *  Created on:   2013-3-19            *
 *  Author:       Qiu Kunlong          *
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
#include "head.h"
#include "gzstream.h"
#include "getsv.h"
#include "clip_reads.h"

using namespace std;


class Alignment 
{
public:
	string chr;
	int pos;
	string left_seq;
	string left_qual;
	string right_seq;
	string right_qual;
	vector<pair<int, char> > cigar_vec;
	char clipped_side;
	char strand;
	
	Alignment() { }
	Alignment(string chr1, int pos1, string left_seq1, string left_qual1, string right_seq1, string right_qual1, vector<pair<int, char> > &cigar_vec1, char clipped_side1, char strand1)
	{
		chr = chr1;
		pos = pos1;
		left_seq = left_seq1;
		left_qual = left_qual1;
		right_seq = right_seq1;
		right_qual = right_qual1;
		cigar_vec = cigar_vec1;
		clipped_side = clipped_side1;
		strand = strand1;
	}
	Alignment& operator= (const Alignment &alignment)
	{
		chr = alignment.chr;
		pos = alignment.pos;
		left_seq = alignment.left_seq;
		left_qual = alignment.left_qual;
		right_seq = alignment.right_seq;
		right_qual = alignment.right_qual;
		cigar_vec = alignment.cigar_vec;
		clipped_side = alignment.clipped_side;
		strand = alignment.strand;
		return *this;
	}

	void set_value(string chr1, int pos1, string left_seq1, string left_qual1, string right_seq1, string right_qual1, vector<pair<int, char> > &cigar_vec1, char clipped_side1, char strand1)
	{
		chr = chr1;
		pos = pos1;
		left_seq = left_seq1;
		left_qual = left_qual1;
		right_seq = right_seq1;
		right_qual = right_qual1;
		cigar_vec = cigar_vec1;
		clipped_side = clipped_side1;
		strand = strand1;
	}

};


void FindJunction(string file, int min_mapQ, multimap<Junction, OtherInfo> &junction2other);
