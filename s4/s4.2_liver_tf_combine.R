library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(org.Hs.eg.db)
library(glmnet)
library(limma)

# combine
load('~/GTRD/final_gtrd/liver/liver_hts1.RData')
load('~/GTRD/final_gtrd/liver/liver_hts2.RData')
load('~/GTRD/final_gtrd/liver/liver_hts3.RData')
load('~/GTRD/final_gtrd/liver/liver_hts4.RData')
load('~/GTRD/final_gtrd/liver/liver_hts5.RData')
load('~/GTRD/final_gtrd/liver/liver_hts6.RData')
load('~/GTRD/final_gtrd/liver/liver_hts7.RData')
load('~/GTRD/final_gtrd/liver/liver_hts8.RData')
load('~/GTRD/final_gtrd/liver/liver_hts9.RData')
load('~/GTRD/final_gtrd/liver/liver_hts10.RData')
load('~/GTRD/final_gtrd/liver/liver_hts11.RData')
load('~/GTRD/final_gtrd/liver/liver_hts12.RData')
load('~/GTRD/final_gtrd/liver/liver_hts_h1.RData')
load('~/GTRD/final_gtrd/liver/liver_hts_s1.RData')


liver_hts <- do.call('rbind', list(liver_hts1 ,liver_hts2,  liver_hts3, liver_hts4,liver_hts5, liver_hts6, liver_hts7, liver_hts8,liver_hts9, liver_hts10, liver_hts11, liver_hts12, liver_hts_h1, liver_hts_s1))
liver_hts$queryHits <- toupper(liver_hts$queryHits)
liver_hts$queryHits <- str_split_fixed(liver_hts$queryHits, pattern = '_|\\.|\\{|\\,|\\(',2)[,1]
liver_hts <- unique(liver_hts)

#### alias to gene symbol
TF <- unique(liver_hts$queryHits)
TF <- as.data.frame(TF)
for(i in 1:nrow(TF)){
  TF[i,'sym'] <-  alias2Symbol(TF[i, 'TF'],  species = "Hs",expand.symbols = F)[1]
}
rownames(TF) <- TF$TF

liver_hts$queryHits <- TF[liver_hts$queryHits,'sym']
#liver_hts$subjectHits <- GENES[liver_hts$subjectHits , 'sym']
liver_hts <- subset( liver_hts, !is.na(liver_hts$queryHits) & !is.na(liver_hts$subjectHits))
liver_hts$subjectHits <- as.vector(liver_hts$subjectHits)
liver_hts_bu <- liver_hts
liver_hts <- liver_hts_bu


#### GTEX
#Genes were filtered to include only those with transcripts per million reads (TPM) > 0.1 in at least 25% of samples,
load('~/GTRD/GTEX/gtex_tissue_exp.RData')
liver_hts <- subset(liver_hts,queryHits%in% rownames(gtex_liver_exp) & subjectHits%in% rownames(gtex_liver_exp))

# Pearson correlation |r|>0.25

func_pearson <- function(i,liver_hts,gtex_liver_exp){
  
  pear_cor <-  cor(as.numeric(gtex_liver_exp[liver_hts[i,1],]), as.numeric(gtex_liver_exp[liver_hts[i,2],]), method=c("pearson"))
  res <- c(liver_hts[i,],pear_cor)
  return(res)
}

#detectCores()
cl <- makeCluster(12) # 初始化16核心集群
#results_pearson <- parLapply(cl,1:200,liver_hts=liver_hts,gtex_liver_exp=gtex_liver_exp, func_pearson) # lapply的并行版本
results_pearson <- parLapply(cl,1:nrow(liver_hts),liver_hts=liver_hts,gtex_liver_exp=gtex_liver_exp, func_pearson) # lapply的并行版本
stopCluster(cl) # 关闭集群'

liver_hts_p <- as.data.frame(matrix(unlist(results_pearson),ncol=3, byrow = T))
colnames(liver_hts_p) <- c('queryHits','subjectHits','cor')
liver_hts_p$cor <- as.numeric(as.vector(liver_hts_p$cor))
liver_hts_p <- subset(liver_hts_p, abs(cor)>0.25)
liver_hts_p$queryHits<- as.character(liver_hts_p$queryHits)
liver_hts_p$subjectHits<- as.character(liver_hts_p$subjectHits)


# lasso regression: 
# independent variable: TFs; dependent variable: genes
t_matrix <- t(gtex_liver_exp)
liver_hts_l <- NULL
r_sq <- c()

for (genes in intersect(names(table(liver_hts_p$subjectHits))[table(liver_hts_p$subjectHits)!=1], colnames(t_matrix) )){
  
  tfs <- liver_hts_p$queryHits[  liver_hts_p$subjectHits==genes ]
  x <- as.matrix(t_matrix[, tfs])
  y <- as.matrix(t_matrix[, genes])#响应变量
  
  cv_glmnet <- cv.glmnet(x,y,alpha=1, family='gaussian', type.measure='mse')
  lam <- cv_glmnet$lambda.1se
  
  r_sq <-  c(r_sq, c(genes,  1 - (cv_glmnet$cvm[which(cv_glmnet$lambda == lam)])/var(y)) )
  
  fit = glmnet(x, y, intercept = F, standardize = F, lambda = lam)
  
  l <- subset(liver_hts, subjectHits== genes & queryHits %in% tfs[as.numeric(coef(fit))!=0] )
  
  liver_hts_l <- rbind(liver_hts_l,l)
  
}

r_sq <- as.data.frame(matrix(r_sq, ncol=2, byrow=T))

r_sq$V1 <- as.character(as.vector(r_sq$V1))
r_sq$V2 <- as.numeric(as.vector(r_sq$V2))
liver_r_sq <- r_sq

# lasso r2>0.5
liver_hts_l2 <- rbind(liver_hts_l, subset(liver_hts_p[,1:2], subjectHits%in%names(table(liver_hts_p$subjectHits))[table(liver_hts_p$subjectHits)==1] ))
liver_hts_p_l <- merge(liver_hts_p,liver_hts_l,all.y=T)
liver_hts_p_r2 <- merge(liver_hts_p_l, liver_r_sq, by.x='subjectHits',by.y='V1',all.x=T)
colnames(liver_hts_p_r2)[4] <- 'r2'

liver_hts_p_r2_l <-  subset(liver_hts_p_r2, r2>0.5 )

