library(Hmisc)
library(ggplot2)

# =pathwaw enrichment
# 2021.01.16

setwd('~/Desktop/tissue_disease_enrich_go/14 disorders tissue GO/')
options(stringsAsFactors = F)

for (ff in dir('/Users/songlt/Desktop/tissue_disease_enrich_go/14 disorders tissue GO/')){
  
  go <- read.csv(ff,stringsAsFactors = F)
  go$X <- capitalize(gsub('^go ','',gsub('_', ' ', tolower(go$X) )))
  mgo <- melt(go)
  colnames(mgo) <- c('Pathway','Tissue','pvalue')
  mgo$Tissue <- gsub('_',' ',mgo$Tissue)
  mgo$`Tissue type`  <- 'Brain'
  mgo$`Tissue type`[mgo$Tissue%in%c( 'Liver' ,'Lung', 'Small bowel')] <- 'Other tissues'
  
  
  mgo <- mgo[order(mgo$pvalue, mgo$Tissue),] #yes
  #mgo <- mgo[order( mgo$Tissue, mgo$pvalue),]
  
  mgo$Pathway <- factor(mgo$Pathway, levels = unique(mgo$Pathway) )
  
  p1 <- ggplot(data=mgo)+
    geom_point(aes( x = Pathway,y = -log10(pvalue), color = Tissue,shape=`Tissue type`,size=-log10(pvalue) )) +
    coord_flip()  +theme_bw()+
    scale_color_manual(values=c("#9900FF",'#D55E00',"#0033FF","#66CCFF",'#999900','#dd9966')) +
    xlab('') + geom_hline(yintercept=-log10(0.05),colour = 'red', linetype = "dashed")+
    scale_x_discrete(labels=function(x) str_wrap(x, width=60))+
    guides(size = guide_legend(order = 3),color = guide_legend(order = 1), shape = guide_legend(order = 2))
    
  tissue <- unlist(strsplit(split = '\\.', ff))[3]
  
  ggsave(p1, file=paste0('~/Desktop/tissue_disease_enrich_go/',tissue,'.pdf'))
  
}





