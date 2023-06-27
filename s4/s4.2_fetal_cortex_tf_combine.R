library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
#library(zFPKM)
#library(lars)
library(org.Hs.eg.db)
library('glmnet')
library(limma)


# combine
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts1.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts2.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts3.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts4.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts5.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts6.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts7.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts8.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts9.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts10.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts11.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts12.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts_h1.RData')
load('~/GTRD/final_gtrd/fetal_brain/fetal_brain_hts_s1.RData')


fetal_brain_hts <- do.call('rbind', list(fetal_brain_hts1 ,fetal_brain_hts2,  fetal_brain_hts3, fetal_brain_hts4,fetal_brain_hts5, fetal_brain_hts6, fetal_brain_hts7, fetal_brain_hts8,fetal_brain_hts9, fetal_brain_hts10, fetal_brain_hts11, fetal_brain_hts12, fetal_brain_hts_h1, fetal_brain_hts_s1))
fetal_brain_hts$queryHits <- toupper(fetal_brain_hts$queryHits)
fetal_brain_hts$queryHits <- str_split_fixed(fetal_brain_hts$queryHits, pattern = '_|\\.|\\{|\\,|\\(',2)[,1]
fetal_brain_hts <- unique(fetal_brain_hts)

#### alias to gene symbol
TF <- unique(fetal_brain_hts$queryHits)
TF <- as.data.frame(TF)
for(i in 1:nrow(TF)){
  TF[i,'sym'] <-  alias2Symbol(TF[i, 'TF'],  species = "Hs",expand.symbols = F)[1]
}
rownames(TF) <- TF$TF

fetal_brain_hts$queryHits <- TF[fetal_brain_hts$queryHits,'sym']
#fetal_brain_hts$subjectHits <- GENES[fetal_brain_hts$subjectHits , 'sym']
fetal_brain_hts <- subset(fetal_brain_hts, !is.na(fetal_brain_hts$queryHits) & !is.na(fetal_brain_hts$subjectHits))
fetal_brain_hts$subjectHits <- as.vector(fetal_brain_hts$subjectHits)
fetal_brain_hts_bu <- fetal_brain_hts
fetal_brain_hts <- fetal_brain_hts_bu


#### GTEX
## Genes were filtered to include only those with transcripts per million reads (TPM) > 0.1 in at least 25% of samples,
load('~/GTRD/GTEX/gtex_tissue_exp.RData')
fetal_brain_hts <- subset(fetal_brain_hts,queryHits%in% rownames(gtex_fetal_brain_exp) & subjectHits%in% rownames(gtex_fetal_brain_exp))

# Pearson correlation |r|>0.25

func_pearson <- function(i,fetal_brain_hts,gtex_fetal_brain_exp){
  
  pear_cor <-  cor(as.numeric(gtex_fetal_brain_exp[fetal_brain_hts[i,1],]), as.numeric(gtex_fetal_brain_exp[fetal_brain_hts[i,2],]), method=c("pearson"))
  res <- c(fetal_brain_hts[i,],pear_cor)
  return(res)
}

#detectCores()
cl <- makeCluster(12) # 
results_pearson <- parLapply(cl,1:nrow(fetal_brain_hts),fetal_brain_hts=fetal_brain_hts,gtex_fetal_brain_exp=gtex_fetal_brain_exp, func_pearson) # lapply的并行版本
stopCluster(cl) 


fetal_brain_hts_p <- as.data.frame(matrix(unlist(results_pearson),ncol=3, byrow = T))
colnames(fetal_brain_hts_p) <- c('queryHits','subjectHits','cor')
fetal_brain_hts_p$cor <- as.numeric(as.vector(fetal_brain_hts_p$cor))
fetal_brain_hts_p <- subset(fetal_brain_hts_p, abs(cor)>0.25)
fetal_brain_hts_p$queryHits<- as.character(fetal_brain_hts_p$queryHits)
fetal_brain_hts_p$subjectHits<- as.character(fetal_brain_hts_p$subjectHits)

# lasso regression: 
# independent variable: TFs; dependent variable: genes

t_matrix <- t(gtex_fetal_brain_exp)
fetal_brain_hts_l <- NULL
r_sq <- c()

for (genes in names(table(fetal_brain_hts_p$subjectHits))[table(fetal_brain_hts_p$subjectHits)>1] ){
  print(genes)
  tfs <- fetal_brain_hts_p$queryHits[  fetal_brain_hts_p$subjectHits==genes ]
  x <- as.matrix(t_matrix[, tfs])
  y <- as.matrix(t_matrix[, genes])
  
  cv_glmnet <- cv.glmnet(x,y,alpha=1, family='gaussian', type.measure='mse')
  lam <- cv_glmnet$lambda.1se
  
  r_sq <-  c(r_sq, c(genes,  1 - (cv_glmnet$cvm[which(cv_glmnet$lambda == lam)])/var(y)) )
  
  fit = glmnet(x, y, intercept = F, standardize = F, lambda = lam)
  
  l <- subset(fetal_brain_hts, subjectHits== genes & queryHits %in% tfs[as.numeric(coef(fit))!=0] )
  
  fetal_brain_hts_l <- rbind(fetal_brain_hts_l,l)
  
}

r_sq <- as.data.frame(matrix(r_sq, ncol=2, byrow=T))

r_sq$V1 <- as.character(as.vector(r_sq$V1))
r_sq$V2 <- as.numeric(as.vector(r_sq$V2))
fetal_brain_r_sq <- r_sq

# lasso r2>0.5
fetal_brain_hts_l2 <- rbind(fetal_brain_hts_l, subset(fetal_brain_hts_p[,1:2], subjectHits%in%names(table(fetal_brain_hts_p$subjectHits))[table(fetal_brain_hts_p$subjectHits)==1] ))
fetal_brain_hts_p_l <- merge(fetal_brain_hts_p,fetal_brain_hts_l,all.y=T)
fetal_brain_hts_p_r2 <- merge(fetal_brain_hts_p_l, fetal_brain_r_sq, by.x='subjectHits',by.y='V1',all.x=T)
colnames(fetal_brain_hts_p_r2)[4] <- 'r2'


fetal_brain_hts_p_r2_l <-  subset(fetal_brain_hts_p_r2, r2>0.5 )


