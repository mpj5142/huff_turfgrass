#!/bin/bash

#kallisto index -i kallisto_maize.idx transps_final.fa

for i in {1..16}
do
	SAMPLE=Sample_$i
	kallisto quant -i assembly/transps_maize/kallisto_maize.idx -o kallisto_transps/kallisto_$SAMPLE -l 200 -s 1 -b 15 -t 5 --single samples/$SAMPLE\/$SAMPLE\_trim.fastq.gz

done