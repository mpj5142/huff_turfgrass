#!/bin/bash

for i in {2..16}
do
	RAWPATH=~/work/data/turfgrass/david-turfgrass-RNAseq
	ADAPATH=~/work/src/Trimmomatic-0.35/adapters
	DIR=Sample_$i
	mkdir $DIR
	cd $DIR
	trimmomatic SE -threads 8 -phred33 -trimlog $DIR\.log $RAWPATH/$DIR/$i\_R1.fastq.gz $DIR\_trim.fastq.gz ILLUMINACLIP:$ADAPATH/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:18 MINLEN:50
	fastqc --nogroup $DIR\_trim.fastq.gz
	cd ..
done

