library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(org.Hs.eg.db)
library(glmnet)
library(limma)


# combine
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts1.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts2.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts3.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts4.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts5.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts6.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts7.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts8.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts9.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts10.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts11.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts12.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts_h1.RData')
load('~/GTRD/final_gtrd/adult_brain/adult_brain_hts_s1.RData')


adult_brain_hts <- do.call('rbind', list(adult_brain_hts1 ,adult_brain_hts2,  adult_brain_hts3, adult_brain_hts4,adult_brain_hts5, adult_brain_hts6, adult_brain_hts7, adult_brain_hts8,adult_brain_hts9, adult_brain_hts10, adult_brain_hts11, adult_brain_hts12, adult_brain_hts_h1, adult_brain_hts_s1))
adult_brain_hts$queryHits <- toupper(adult_brain_hts$queryHits)
adult_brain_hts$queryHits <- str_split_fixed(adult_brain_hts$queryHits, pattern = '_|\\.|\\{|\\,|\\(',2)[,1]
adult_brain_hts <- unique(adult_brain_hts)

#### alias to gene symbol
TF <- unique(adult_brain_hts$queryHits)
TF <- as.data.frame(TF)
for(i in 1:nrow(TF)){
  TF[i,'sym'] <-  alias2Symbol(TF[i, 'TF'],  species = "Hs",expand.symbols = F)[1]
}
rownames(TF) <- TF$TF


adult_brain_hts$queryHits <- TF[adult_brain_hts$queryHits,'sym']
#adult_brain_hts$subjectHits <- GENES[adult_brain_hts$subjectHits , 'sym']
adult_brain_hts <- subset( adult_brain_hts, !is.na(adult_brain_hts$queryHits) & !is.na(adult_brain_hts$subjectHits))
adult_brain_hts$subjectHits <- as.vector(adult_brain_hts$subjectHits)
adult_brain_hts_bu <- adult_brain_hts
adult_brain_hts <- adult_brain_hts_bu


#### GTEX
##  Genes were filtered to include only those with transcripts per million reads (TPM) > 0.1 in at least 25% of samples,
load('~/GTRD/GTEX/gtex_tissue_exp.RData')
adult_brain_hts <- subset(adult_brain_hts,queryHits%in% rownames(gtex_adult_brain_exp) & subjectHits%in% rownames(gtex_adult_brain_exp))


# Pearson correlation |r|>0.25

func_pearson <- function(i,adult_brain_hts,gtex_adult_brain_exp){
  
  pear_cor <-  cor(as.numeric(gtex_adult_brain_exp[adult_brain_hts[i,1],]), as.numeric(gtex_adult_brain_exp[adult_brain_hts[i,2],]), method=c("pearson"))
  res <- c(adult_brain_hts[i,],pear_cor)
  return(res)
}

#detectCores()
cl <- makeCluster(12) 
results_pearson <- parLapply(cl,1:nrow(adult_brain_hts),adult_brain_hts=adult_brain_hts,gtex_adult_brain_exp=gtex_adult_brain_exp, func_pearson) # lapply的并行版本
stopCluster(cl) 

adult_brain_hts_p <- as.data.frame(matrix(unlist(results_pearson),ncol=3, byrow = T))
colnames(adult_brain_hts_p) <- c('queryHits','subjectHits','cor')
adult_brain_hts_p$cor <- as.numeric(as.vector(adult_brain_hts_p$cor))

adult_brain_hts_p <- subset(adult_brain_hts_p, abs(cor)>0.25)
adult_brain_hts_p$queryHits<- as.character(adult_brain_hts_p$queryHits)
adult_brain_hts_p$subjectHits<- as.character(adult_brain_hts_p$subjectHits)

# lasso regression: 
# independent variable: TFs; dependent variable: genes
t_matrix <- t(gtex_adult_brain_exp)
adult_brain_hts_l <- NULL
r_sq <- c()

for (genes in intersect(names(table(adult_brain_hts_p$subjectHits))[table(adult_brain_hts_p$subjectHits)!=1], colnames(t_matrix) )){

  tfs <- adult_brain_hts_p$queryHits[  adult_brain_hts_p$subjectHits==genes ]
  x <- as.matrix(t_matrix[, tfs])
  y <- as.matrix(t_matrix[, genes])

  cv_glmnet <- cv.glmnet(x,y,alpha=1, family='gaussian', type.measure='mse')
  lam <- cv_glmnet$lambda.1se

  r_sq <-  c(r_sq, c(genes,  1 - (cv_glmnet$cvm[which(cv_glmnet$lambda == lam)])/var(y)) )

  fit = glmnet(x, y, intercept = F, standardize = F, lambda = lam)

  l <- subset(adult_brain_hts, subjectHits== genes & queryHits %in% tfs[as.numeric(coef(fit))!=0] )

  adult_brain_hts_l <- rbind(adult_brain_hts_l,l)

}

r_sq <- as.data.frame(matrix(r_sq, ncol=2, byrow=T))

r_sq$V1 <- as.character(as.vector(r_sq$V1))
r_sq$V2 <- as.numeric(as.vector(r_sq$V2))
adult_brain_r_sq <- r_sq

# lasso r2>0.5
adult_brain_hts_l2 <- rbind(adult_brain_hts_l, subset(adult_brain_hts_p[,1:2], subjectHits%in%names(table(adult_brain_hts_p$subjectHits))[table(adult_brain_hts_p$subjectHits)==1] ))
adult_brain_hts_p_l <- merge(adult_brain_hts_p,adult_brain_hts_l,all.y=T)
adult_brain_hts_p_r2 <- merge(adult_brain_hts_p_l, adult_brain_r_sq, by.x='subjectHits',by.y='V1',all.x=T)
colnames(adult_brain_hts_p_r2)[4] <- 'r2'

adult_brain_hts_p_r2_l <-  subset(adult_brain_hts_p_r2, r2>0.5 )


