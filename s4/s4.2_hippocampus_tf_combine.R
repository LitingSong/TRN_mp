library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
#library(zFPKM)
#library(lars)
library(org.Hs.eg.db)
library('glmnet')
library(limma)

# combine
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts1.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts2.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts3.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts4.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts5.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts6.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts7.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts8.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts9.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts10.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts11.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts12.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts_h1.RData')
load('~/GTRD/final_gtrd/hippocampus/hippocampus_hts_s1.RData')


hippocampus_hts <- do.call('rbind', list(hippocampus_hts1 ,hippocampus_hts2,  hippocampus_hts3, hippocampus_hts4,hippocampus_hts5, hippocampus_hts6, hippocampus_hts7, hippocampus_hts8,hippocampus_hts9, hippocampus_hts10, hippocampus_hts11, hippocampus_hts12, hippocampus_hts_h1, hippocampus_hts_s1))
hippocampus_hts$queryHits <- toupper(hippocampus_hts$queryHits)
hippocampus_hts$queryHits <- str_split_fixed(hippocampus_hts$queryHits, pattern = '_|\\.|\\{|\\,|\\(',2)[,1]
hippocampus_hts <- unique(hippocampus_hts)

#### alias to gene symbol
TF <- unique(hippocampus_hts$queryHits)
TF <- as.data.frame(TF)

for(i in 1:nrow(TF)){
  TF[i,'sym'] <-  alias2Symbol(TF[i, 'TF'],  species = "Hs",expand.symbols = F)[1]
}
rownames(TF) <- TF$TF

hippocampus_hts$queryHits <- TF[hippocampus_hts$queryHits,'sym']
#hippocampus_hts$subjectHits <- GENES[hippocampus_hts$subjectHits , 'sym']
hippocampus_hts <- subset(hippocampus_hts, !is.na(hippocampus_hts$queryHits) & !is.na(hippocampus_hts$subjectHits))
hippocampus_hts$subjectHits <- as.vector(hippocampus_hts$subjectHits)
hippocampus_hts_bu <- hippocampus_hts
hippocampus_hts <- hippocampus_hts_bu


#### GTEX
##Genes were filtered to include only those with transcripts per million reads (TPM) > 0.1 in at least 25% of samples,
load('~/GTRD/GTEX/gtex_tissue_exp.RData')
hippocampus_hts <- subset(hippocampus_hts,queryHits%in% rownames(gtex_hippocampus_exp) & subjectHits%in% rownames(gtex_hippocampus_exp))

# Pearson correlation |r|>0.25

func_pearson <- function(i,hippocampus_hts,gtex_hippocampus_exp){
  
  pear_cor <-  cor(as.numeric(gtex_hippocampus_exp[hippocampus_hts[i,1],]), as.numeric(gtex_hippocampus_exp[hippocampus_hts[i,2],]), method=c("pearson"))
  res <- c(hippocampus_hts[i,],pear_cor)
  return(res)
}


cl <- makeCluster(12) 
results_pearson <- parLapply(cl,1:nrow(hippocampus_hts),hippocampus_hts=hippocampus_hts,gtex_hippocampus_exp=gtex_hippocampus_exp, func_pearson) # lapply的并行版本
stopCluster(cl) 

hippocampus_hts_p <- as.data.frame(matrix(unlist(results_pearson),ncol=3, byrow = T))
colnames(hippocampus_hts_p) <- c('queryHits','subjectHits','cor')
hippocampus_hts_p$cor <- as.numeric(as.vector(hippocampus_hts_p$cor))
hippocampus_hts_p <- subset(hippocampus_hts_p, abs(cor)>0.25)
hippocampus_hts_p$queryHits<- as.character(hippocampus_hts_p$queryHits)
hippocampus_hts_p$subjectHits<- as.character(hippocampus_hts_p$subjectHits)

# lasso regression: 
# independent variable: TFs; dependent variable: genes

t_matrix <- t(gtex_hippocampus_exp)
hippocampus_hts_l <- NULL
r_sq <- c()

for (genes in intersect(names(table(hippocampus_hts_p$subjectHits))[table(hippocampus_hts_p$subjectHits)!=1], colnames(t_matrix) )){
  
  tfs <- hippocampus_hts_p$queryHits[  hippocampus_hts_p$subjectHits==genes ]
  x <- as.matrix(t_matrix[, tfs])
  y <- as.matrix(t_matrix[, genes])
  
  cv_glmnet <- cv.glmnet(x,y,alpha=1, family='gaussian', type.measure='mse')
  lam <- cv_glmnet$lambda.1se
  
  r_sq <-  c(r_sq, c(genes,  1 - (cv_glmnet$cvm[which(cv_glmnet$lambda == lam)])/var(y)) )
  
  fit = glmnet(x, y, intercept = F, standardize = F, lambda = lam)
  
  l <- subset(hippocampus_hts, subjectHits== genes & queryHits %in% tfs[as.numeric(coef(fit))!=0] )
  
  hippocampus_hts_l <- rbind(hippocampus_hts_l,l)
  
}

r_sq <- as.data.frame(matrix(r_sq, ncol=2, byrow=T))

r_sq$V1 <- as.character(as.vector(r_sq$V1))
r_sq$V2 <- as.numeric(as.vector(r_sq$V2))
hippocampus_r_sq <- r_sq


# lasso r2>0.5

hippocampus_hts_l2 <- rbind(hippocampus_hts_l, subset(hippocampus_hts_p[,1:2], subjectHits%in%names(table(hippocampus_hts_p$subjectHits))[table(hippocampus_hts_p$subjectHits)==1] ))

hippocampus_hts_p_l <- merge(hippocampus_hts_p,hippocampus_hts_l,all.y=T)
hippocampus_hts_p_r2 <- merge(hippocampus_hts_p_l, hippocampus_r_sq, by.x='subjectHits',by.y='V1',all.x=T)
colnames(hippocampus_hts_p_r2)[4] <- 'r2'

hippocampus_hts_p_r2_l <-  subset(hippocampus_hts_p_r2, r2>0.5 )
