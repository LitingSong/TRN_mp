library(limma)
library(stringr)
### filter expressed genes in each tissue based on the GTEx.
### gtex https://xenabrowser.net/datapages/?dataset=gtex_RSEM_gene_tpm&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443


### gtex 
gtex_pheno <- read.table('~/GTRD/GTEX/GTEX_phenotype', sep='\t', header=T, row.names = 1 )
rownames(gtex_pheno) <- gsub('-','.', rownames(gtex_pheno) )
#gtex_gene <- read.table('~/GTRD/GTEX/gencode.v23.annotation.gene.probemap', header=T, row.names = 1 )

grch38 <- read.table('~/GTRD/Homo_sapiens.GRCh38.87.txt', header = T)
rownames(grch38) <- grch38$ensg

gtex_exp <- read.table('~/GTRD/GTEX/gtex_RSEM_gene_tpm',sep='\t',header = T, row.names = 1)

rownames(gtex_exp) <- str_split_fixed(rownames(gtex_exp), '\\.', 2)[,1]

gtex_exp <- gtex_exp[ intersect(rownames(grch38), rownames(gtex_exp)), ]
gtex_exp$gene_symbol <- grch38[rownames(gtex_exp),'symbol']
gtex_exp <- gtex_exp[ !duplicated(gtex_exp$gene_symbol) , ]
rownames(gtex_exp) <- gtex_exp$gene_symbol

gtex_exp <- gtex_exp[, -7863]
gtex_exp <- round(2^gtex_exp,3)-0.001


#  TPM >0.1 in at least 25% samples ；Genes were filtered to include only those with transcripts per million reads (TPM) > 0.1 in at least 25% of samples,
gtex_brain <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Brain'], colnames(gtex_exp)) ]
gtex_brain_exp <- gtex_brain[apply(gtex_brain, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]

gtex_liver <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Liver'], colnames(gtex_exp)) ]
gtex_liver_exp <- gtex_liver[apply(gtex_liver, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]


gtex_lung <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Lung'], colnames(gtex_exp)) ]
gtex_lung_exp <- gtex_lung[apply(gtex_lung, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]


gtex_colon <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Colon'], colnames(gtex_exp)) ]
gtex_colon_exp <- gtex_colon[apply(gtex_colon, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]


gtex_stomach <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Stomach'], colnames(gtex_exp)) ]
gtex_stomach_exp <- gtex_stomach[apply(gtex_stomach, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]

gtex_intestine <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$X_primary_site=='Small Intestine'], colnames(gtex_exp)) ]
gtex_intestine_exp <- gtex_intestine[apply(gtex_intestine, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]

gtex_brain_cortex <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.%in%c('Brain - Cortex','Brain - Frontal Cortex (BA9)')], colnames(gtex_exp)) ]
gtex_brain_cortex_exp <- gtex_brain_cortex[apply(gtex_brain_cortex, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]

gtex_hippocampus <- gtex_exp[, intersect(rownames(gtex_pheno)[gtex_pheno$body_site_detail..SMTSD.%in%c('Brain - Hippocampus')], colnames(gtex_exp)) ]
gtex_hippocampus_exp <- gtex_hippocampus[apply(gtex_hippocampus, 1, function(x){ length(which(x>0.1)) > length(x)*0.25 }),]

save(gtex_exp, gtex_pheno, file='~/GTRD/GTEX/gtex_processed_tpm_symbol.RData')
save(gtex_brain_exp, gtex_liver_exp, gtex_lung_exp, gtex_colon_exp, gtex_stomach_exp, gtex_intestine_exp,gtex_brain_cortex_exp, gtex_hippocampus_exp, file='~/GTRD/GTEX/gtex_tissue_exp.RData')


