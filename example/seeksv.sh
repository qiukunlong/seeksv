./bin/seeksv getclip -o normal normal.sort.bam 
### bwa aln + bwa sampe and bwa bwasw also OK
./bin/bwa mem ./reference/example.fa normal.clip.fq.gz |./bin/samtools view -Sb - >normal.clip.bam
./bin/seeksv getsv normal.clip.bam normal.sort.bam normal.clip.gz normal.sv normal.clipunmap.gz
