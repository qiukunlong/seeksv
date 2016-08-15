seeksv: bam2depth.cpp seeksv.cpp bam2depth.cpp cluster.cpp gzstream.C seeksv.cpp clip_reads.cpp  getsv.cpp somatic.cpp process_bwasw.cpp clip_reads.h getsv.h somatic.h process_bwasw.h
	g++ bam2depth.cpp cluster.cpp gzstream.C seeksv.cpp clip_reads.cpp  getsv.cpp somatic.cpp process_bwasw.cpp -o seeksv  -lz -lpthread -lm -L/samtools/0.1.18/src/samtools -lbam

clean:
	rm *.o seeksv
