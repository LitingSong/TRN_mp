#!/bin/bash

# Genome scan of TF motifs

## for motifs in JASPAR 2020
~/miniconda3/bin/fimo --text  --skip-matched-sequence ./JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt  /home1/songlt/miniconda3/share/homer-4.10-0/data/genomes/hg38/genome.fa > JASPAR.tsv

### Here, to speed up the scaning process, we divided the JASPAR2020 meme file into 12 sub-files.
~/miniconda3/bin/fimo --text  --skip-matched-sequence ./J1-j12.txt  /home1/songlt/miniconda3/share/homer-4.10-0/data/genomes/hg38/genome.fa > j1-12.tsv


## for motifs in HOCOMOCO
~/miniconda3/bin/fimo --text --skip-matched-sequence ~/GTRD/ho_addi_motif.txt  /home1/songlt/miniconda3/share/homer-4.10-0/data/genomes/hg38/genome.fa > h1.tsv

## for motifs in SwissRegulon
~/miniconda3/bin/fimo --text --skip-matched-sequence ~/GTRD/swiss_addi_motif.txt   /home1/songlt/miniconda3/share/homer-4.10-0/data/genomes/hg38/genome.fa > s1.tsv

 
