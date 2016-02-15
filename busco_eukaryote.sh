#!/bin/bash

## Note: Trial frun for BUSCO on the turfgrass transcriptome assembly. For a reference ortholog file, I'm using a eukaryote-specific file;
## BUSCO may release a plantonly file at some point in the future that may be better for the turfgrass.

python ~/work/src/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o turfgrass_out -in trinity_out_dir/Trinity.fasta -l eukaryota/ -m trans