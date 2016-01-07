seeksv: bam2depth.cpp seeksv.cpp bam2depth.cpp cluster.cpp gzstream.C seeksv.cpp clip_reads.cpp  getsv.cpp somatic.cpp process_bwasw.cpp clip_reads.h getsv.h somatic.h process_bwasw.h
	#g++ bam2depth.cpp gzstream.C seeksv.cpp clip_reads.cpp  getsv.cpp somatic.cpp -o seeksv -lz -lpthread -lm -L/share/backup/qiukl/lib -lbam
	#g++ -g -rdynamic   bam2depth.cpp cluster.cpp gzstream.C seeksv.cpp clip_reads.cpp  getsv.cpp somatic.cpp process_bwasw.cpp -o seeksv  -lz -lpthread -lm -L/ifs4/BC_CANCER/02usr/hej/src/samtools/samtools-0.1.18 -lbam
	g++ bam2depth.cpp cluster.cpp gzstream.C seeksv.cpp clip_reads.cpp  getsv.cpp somatic.cpp process_bwasw.cpp -o seeksv  -lz -lpthread -lm -L/ifs4/BC_CANCER/02usr/qiukl/src/samtools-0.1.18 -lbam

clean:
	rm *.o seeksv
