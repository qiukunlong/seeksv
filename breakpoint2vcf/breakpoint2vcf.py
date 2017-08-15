#!/usr/bin/python
import sys
import argparse
import vcf


def parser():
	parser = argparse.ArgumentParser(prog="breakpoint2vcf.py", description="Change seeksv results to vcf")
	parser.add_argument("breakpoint", help="Input breakpoint file")
	parser.add_argument("vcf_file", help="Output vcf file")
	return parser

base_reverse_complementary = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'a' : 'T', 't' : 'A', 'c' : 'G', 'g' : 'C'}


def breakpoint2vcfRecord(breakpoint_dict, breakpoint_id):	
	if breakpoint_dict['left_strand'] == '+' and breakpoint_dict['right_strand'] == '+':
		ref1 = breakpoint_dict['left_seq'][-1]
		alt1 = vcf.model._Breakend(breakpoint_dict['right_chr'], breakpoint_dict['right_pos'], False, True, ref1, None)
		ref2 = breakpoint_dict['right_seq'][0]
		alt2 = vcf.model._Breakend(breakpoint_dict['left_chr'], breakpoint_dict['left_pos'], True, False, ref2, None)
	elif breakpoint_dict['left_strand'] == '+' and breakpoint_dict['right_strand'] == '-':
		ref1 = breakpoint_dict['left_seq'][-1]
		alt1 = vcf.model._Breakend(breakpoint_dict['right_chr'], breakpoint_dict['right_pos'], True, True, ref1, None)
		ref2 = base_reverse_complementary[breakpoint_dict['right_seq'][0]]
		alt2 = vcf.model._Breakend(breakpoint_dict['left_chr'], breakpoint_dict['left_pos'], True, False, ref2, None)
		
	record1 = vcf.model._Record(\
			CHROM=breakpoint_dict['left_chr'], \
			ID="bnt" + str(breakpoint_id) + "_U", \
			POS=breakpoint_dict['left_pos'],  \
			REF=ref1, \
			ALT=[alt1], \
			QUAL='.', \
			FILTER='PASS', \
			INFO={}, \
			FORMAT='.', \
			sample_indexes='.')

	record2 = vcf.model._Record(\
			CHROM=breakpoint_dict['right_chr'], \
			ID="bnt" + str(breakpoint_id) + "_D", \
			POS=breakpoint_dict['right_pos'],  \
			REF=ref2, \
			ALT=[alt2], \
			QUAL='.', \
			FILTER='PASS', \
			INFO={}, \
			FORMAT='.', \
			sample_indexes='.')
	return record1, record2

def main():
	args = parser().parse_args(sys.argv[1:])
	in_vcf = args.breakpoint
	out_vcf = args.vcf_file
	vcf_writer = vcf.Writer(open(out_vcf, 'w'), vcf.Reader(filename=in_vcf))
	record = vcf.model._Record(CHROM=20, ID='.', POS=1110696, REF='A', ALT=[vcf.model._Substitution('G'), vcf.model._Substitution('T')], QUAL='.', FILTER='.', INFO={}, FORMAT='.', sample_indexes='.')
	vcf_writer.write_record(record)

	return 0

if __name__ == '__main__':
	main()
