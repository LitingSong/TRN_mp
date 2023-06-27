# Disrupted long-range gene regulations elucidate shared tissue-specific mechanisms of neuropsychiatric disorders


## Step1: Genome scan of TF motifs.

s1.1 determine the genome (GECh38) positions for TF motifs in JASPAR CORE Vertebrate 2020, HOMOCOMO and Swiss Regulon. To speed up the scaning process, we divided the JASPAR2020 meme file into 12 sub-files (i.e., j1-j12.txt).
./s1/s1_motif_genome_scan.sh

## Step2: get the DNase footprints of each tissue from the GTRD (the Gene Transcription Regulation Database).

s2.1: get peak ids for each tissue.
./s2/s2.1_dnase_seq_file_cellline.sh

s2.2: combine footprint bed files for each tissue.
./s2/s2.2_dnase_seq_file_cellline.sh

## Step3: Get TF-gene pairs for each tissue.

s3.1 overlap the predicted TF binding regions (outputs of step1) with the DNase footprints of each tissue (outputs of step2).

s3.2 overlap predicted motif region of a TF with the +/-10kb region surrounding a gene's transcriptional starting site.

./s3/*tissues*/tf1-12.R: correspond to predicted TF binding regions using 12 JASPAR motif subfiles.
./s3/*tissues*/tf_h1.R and tf_s1.R correspond to HOMOCOMO and Swiss Regulon motif file, seperately.

## Step4: filter TF-gene pair for each tissue.

4.1 get expressed genes (TPM >0.1 in at least 25% samples) in each tissue using GTEx data.
./s4/s4.1_gtex_exp.R 

4.2 filter TF-gene pair based on expressed genes and Pearson correlation and Lasso regression.

./s4/s4.2_*tissue*_tf_combine.R; 

## Step5: WGCNA using gtex data

Weighted correlation netwok analysis (WGCNA) for each tissue
./s5/wgcna_expgene.R

enrichment of each tissue in neuropsychiatric disorders pathways
./s5/plot_tissue_disease.R



## Contact
Dr. Jingqi Chen (jingqichen@fudan.edu.cn);Dr. Liting Song (ltsong18@fudan.edu.cn)


