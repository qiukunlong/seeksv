/*
 ***************************************
 * cluster.h                           *
 *                                     *
 *  Created on:   2013-9-9             *
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

using namespace std;


bool CalculateInsertsizeDeviation(string file, int min_mapQ, int read_pair_used, int &mean_insert_size, int &deviation);

bool DivideCluster(string file, int min_mapQ, int &mean_insert_size, int &deviation);

bool IsConcordant(const bam_header_t *header, const bam1_t *b, int mean_insert_size, int deviation, int times);

bool RandomFetch(string file, int min_mapQ, int &mean_insert_size, int &deviation);

//int BamGetTid(const bam_header_t *header, const char *seq_name);
int StoreSeqName2Tid(const bam_header_t *header, map<string, int> &seq_name2tid);

int BamGetTid(map<string, int> &seq_name2tid, const string seq_name);
