#!/bin/bash

makeblastdb -in transps/maize_protein.fa -dbtype prot

blastx -db transps/maize_protein.fa -query trinity_out_dir/Trinity.fasta -out transps/blastx.out -outfmt 6 -evalue 0.01 -max_target_seqs 20

perl /home/matt/work/src/TransPS1.1.0/transps.pl -t trinity_out_dir/Trinity.fasta -b transps/blastx.out