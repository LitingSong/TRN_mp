library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
#library(zFPKM)
#library(lars)
library(org.Hs.eg.db)
library('glmnet')
library(limma)

# combine
load('~/GTRD/final_gtrd/intestine/intestine_hts1.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts2.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts3.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts4.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts5.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts6.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts7.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts8.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts9.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts10.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts11.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts12.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts_h1.RData')
load('~/GTRD/final_gtrd/intestine/intestine_hts_s1.RData')


intestine_hts <- do.call('rbind', list(intestine_hts1 ,intestine_hts2,  intestine_hts3, intestine_hts4,intestine_hts5, intestine_hts6, intestine_hts7, intestine_hts8,intestine_hts9, intestine_hts10, intestine_hts11, intestine_hts12, intestine_hts_h1, intestine_hts_s1))
intestine_hts$queryHits <- toupper(intestine_hts$queryHits)
intestine_hts$queryHits <- str_split_fixed(intestine_hts$queryHits, pattern = '_|\\.|\\{|\\,|\\(',2)[,1]
intestine_hts <- unique(intestine_hts)

#### alias to gene symbol
TF <- unique(intestine_hts$queryHits)
TF <- as.data.frame(TF)

for(i in 1:nrow(TF)){
  TF[i,'sym'] <-  alias2Symbol(TF[i, 'TF'],  species = "Hs",expand.symbols = F)[1]
}
rownames(TF) <- TF$TF

intestine_hts$queryHits <- TF[intestine_hts$queryHits,'sym']
#intestine_hts$subjectHits <- GENES[intestine_hts$subjectHits , 'sym']
intestine_hts <- subset(intestine_hts, !is.na(intestine_hts$queryHits) & !is.na(intestine_hts$subjectHits))
intestine_hts$subjectHits <- as.vector(intestine_hts$subjectHits)
intestine_hts_bu <- intestine_hts
intestine_hts <- intestine_hts_bu


#### GTEX
##  至少25%的样本中>0.1 TPM ；Genes were filtered to include only those with transcripts per million reads (TPM) > 0.1 in at least 25% of samples,
load('~/GTRD/GTEX/gtex_tissue_exp.RData')
intestine_hts <- subset(intestine_hts,queryHits%in% rownames(gtex_intestine_exp) & subjectHits%in% rownames(gtex_intestine_exp))

# Pearson correlation |r|>0.25

func_pearson <- function(i,intestine_hts,gtex_intestine_exp){
  
  pear_cor <-  cor(as.numeric(gtex_intestine_exp[intestine_hts[i,1],]), as.numeric(gtex_intestine_exp[intestine_hts[i,2],]), method=c("pearson"))
  res <- c(intestine_hts[i,],pear_cor)
  return(res)
}

#detectCores()
cl <- makeCluster(12) 
results_pearson <- parLapply(cl,1:nrow(intestine_hts),intestine_hts=intestine_hts,gtex_intestine_exp=gtex_intestine_exp, func_pearson) # lapply的并行版本
stopCluster(cl) 

intestine_hts_p <- as.data.frame(matrix(unlist(results_pearson),ncol=3, byrow = T))
colnames(intestine_hts_p) <- c('queryHits','subjectHits','cor')
intestine_hts_p$cor <- as.numeric(as.vector(intestine_hts_p$cor))
intestine_hts_p <- subset(intestine_hts_p, abs(cor)>0.25)
intestine_hts_p$queryHits<- as.character(intestine_hts_p$queryHits)
intestine_hts_p$subjectHits<- as.character(intestine_hts_p$subjectHits)

# lasso regression: 
# independent variable: TFs; dependent variable: genes
t_matrix <- t(gtex_intestine_exp)
intestine_hts_l <- NULL
r_sq <- c()

for (genes in intersect(names(table(intestine_hts_p$subjectHits))[table(intestine_hts_p$subjectHits)!=1], colnames(t_matrix) )){
  
  tfs <- intestine_hts_p$queryHits[  intestine_hts_p$subjectHits==genes ]
  x <- as.matrix(t_matrix[, tfs])
  y <- as.matrix(t_matrix[, genes])#响应变量
  
  cv_glmnet <- cv.glmnet(x,y,alpha=1, family='gaussian', type.measure='mse')
  lam <- cv_glmnet$lambda.1se
  
  r_sq <-  c(r_sq, c(genes,  1 - (cv_glmnet$cvm[which(cv_glmnet$lambda == lam)])/var(y)) )
  
  fit = glmnet(x, y, intercept = F, standardize = F, lambda = lam)
  
  l <- subset(intestine_hts, subjectHits== genes & queryHits %in% tfs[as.numeric(coef(fit))!=0] )
  
  intestine_hts_l <- rbind(intestine_hts_l,l)
  
}

r_sq <- as.data.frame(matrix(r_sq, ncol=2, byrow=T))

r_sq$V1 <- as.character(as.vector(r_sq$V1))
r_sq$V2 <- as.numeric(as.vector(r_sq$V2))
intestine_r_sq <- r_sq


# lasso r2>0.5
intestine_hts_l2 <- rbind(intestine_hts_l, subset(intestine_hts_p[,1:2], subjectHits%in%names(table(intestine_hts_p$subjectHits))[table(intestine_hts_p$subjectHits)==1] ))
intestine_hts_p_l <- merge(intestine_hts_p,intestine_hts_l,all.y=T)
intestine_hts_p_r2 <- merge(intestine_hts_p_l, intestine_r_sq, by.x='subjectHits',by.y='V1',all.x=T)
colnames(intestine_hts_p_r2)[4] <- 'r2'

intestine_hts_p_r2_l <-  subset(intestine_hts_p_r2, r2>0.5 )
