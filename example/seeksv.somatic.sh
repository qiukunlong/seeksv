./bin/seeksv getclip -o cancer cancer.sort.bam 
### bwa aln + bwa sampe and bwa bwasw also OK
./bin/bwa mem ./reference/example.fa cancer.clip.fq.gz |./bin/samtools view -Sb - >cancer.clip.bam
./bin/seeksv getsv cancer.clip.bam cancer.sort.bam cancer.clip.gz cancer.sv cancer.clipunmap.gz
./bin/seeksv somatic normal.sort.bam normal.clip.gz cancer.sv cancer.somatic.temp.sv
awk '$24 == 0 && $25 == 0 && $26 == 0 || /^@/' cancer.somatic.temp.sv >cancer.somatic.sv
