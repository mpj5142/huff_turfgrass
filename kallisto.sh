#!/bin/bash

kallisto index -i trinity_out_dir/kallisto.idx Trinity.fasta

kallisto quant -i trinity_out_dir/kallisto.idx -o kallisto_out/ -l 200 -s 1  -b 15 -t 5 --single Sample_1/Sample_1_trim.fastq.gz Sample_2/Sample_2_trim.fastq.gz Sample_3/Sample_3_trim.fastq.gz Sample_4/Sample_4_trim.fastq.gz Sample_5/Sample_5_trim.fastq.gz Sample_6/Sample_6_trim.fastq.gz Sample_7/Sample_7_trim.fastq.gz Sample_8/Sample_8_trim.fastq.gz Sample_9/Sample_9_trim.fastq.gz Sample_10/Sample_10_trim.fastq.gz Sample_11/Sample_11_trim.fastq.gz Sample_12/Sample_12_trim.fastq.gz Sample_13/Sample_13_trim.fastq.gz Sample_14/Sample_14_trim.fastq.gz Sample_15/Sample_15_trim.fastq.gz Sample_16/Sample_16_trim.fastq.gz