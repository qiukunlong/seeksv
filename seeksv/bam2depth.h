/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */
#ifndef BAM2DEPTH_H_
#define BAM2DEPTH_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "head.h"
#include <string>
#include <map>
#include "getsv.h"


using namespace std;

typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ filter
} aux_t;


// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}

//int main_depth(map<pair<string, int>, int> &pos2depth, string bam_file, int baseQ, int mapQ, char *reg);
int main_depth(map<pair<string, int>, int> &pos2depth, map<ChrRange, unsigned long> &range2depth, map<pair<string, int>, int> &begin2end, string bam_file, int baseQ, int mapQ, char *reg);
#endif
