#!/bin/bash

##kallisto index -i trinity_out_dir/kallisto.idx Trinity.fasta

for i in {1..16}
do
	SAMPLE=Sample_$i
	kallisto quant -i trinity_reverse/kallisto.idx -o kallisto/kallisto_$SAMPLE -l 200 -s 1 -b 15 -t 5 --single $SAMPLE\/$SAMPLE\_trim.fastq.gz

done