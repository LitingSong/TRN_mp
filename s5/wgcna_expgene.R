# WGCNA

library(WGCNA)
setwd('~/Desktop/WGCNA/exp_gene_100/')
load('~/Desktop/WGCNA/gtex_wgcna.RData')
allowWGCNAThreads()

summary <- c()

for (tissue in c('PFC', 'liver','lung','intestine','hippa')){
  for(cutHeight in c(100,150)){
    
  if (tissue=='PFC'){expro <- gtex_brain_exp;  }
  if (tissue=='liver'){expro <- gtex_Liver_exp;  }
  if (tissue=='lung'){expro <- gtex_Lung_exp;  }
  if (tissue=='intestine'){expro <- gtex_Intestine_exp;  }
  if (tissue=='hippa'){expro <- gtex_hippa_exp;  }
  
  expro <- log2(expro+1)
  #datExpr = t(expro[order(apply(expro,1,mad), decreasing = T)[1:5000],])
  datExpr = t(expro)
  
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  gsg = goodSamplesGenes(datExpr, verbose = 0)
  
  ##  Flagging genes and samples with too many missing values...
  ##   ..step 1
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:",
                       paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  dim(datExpr)
  
  
  sampleTree = hclust(dist(datExpr), method = "average")

  
  clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
  print(table(clust))
  keepSamples = (clust== 1)
  datExpr = datExpr[keepSamples, ]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  
  powers = c(c(4:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers,networkType='signed',RsquaredCut = 0.8)
  sft1 = pickSoftThreshold(datExpr, powerVector = powers,networkType='signed',RsquaredCut = 0.75)
  
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.8;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.8,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  
  #dev.print(pdf, file=paste(tissue, 'sft.pdf',sep=''))
  
  
  if(is.na(sft$powerEstimate) & is.na(sft1$powerEstimate)){power <- 12}
  if(is.na(sft$powerEstimate) & !is.na(sft1$powerEstimate)){power <- sft1$powerEstimate}
  if(!is.na(sft$powerEstimate) & !is.na(sft1$powerEstimate)){power <- sft$powerEstimate}
  
  
  minModuleSize = min(20, ncol(datExpr)/2 )
  cor <- WGCNA::cor

  # step3:Automatic network construction and module detection
  for (mergeCutHeight in c(0.15)){
    for (deepSplit in c(0,1)){
   
      net <-  blockwiseModules(
        datExpr,networkType = "signed",
        power = power,deepSplit=deepSplit,
        detectCutHeight = 0.995, 
        maxBlockSize=2000,
        minModuleSize = min(20, ncol(datExpr)/2 ),
        
        # Gene reassignment, module trimming, and module "significance" criteria
        
        reassignThreshold = 1e-6,
        minCoreKME = 0.5, 
        minCoreKMESize = minModuleSize/3,
        minKMEtoStay = 0.3,
        
        # Module merging options
        
        mergeCutHeight = mergeCutHeight, 
        impute = TRUE, 
        trapErrors = FALSE, 
        
        # Output options
        
        numericLabels = T,
        
        # Options controlling behaviour
        
        useInternalMatrixAlgebra = FALSE,
        useCorOptionsThroughout = TRUE,
        verbose = 0, indent = 0,
      )
      
      table(net$colors)
      
      moduleLabels = net$colors
      moduleColors = labels2colors(moduleLabels)

      modulec <- as.data.frame(cbind(colnames(datExpr),moduleColors ))
      modulec <- modulec[order(modulec$moduleColors),]
      colnames(modulec) <-c('gene','module')
      
      
      write.csv(modulec, file=paste0("MergeCutHeight",mergeCutHeight,"_","deepSplit",deepSplit,'/',tissue,"_",'cutHeight',cutHeight,"_","mergeCutHeight",mergeCutHeight,"_","deepSplit",deepSplit,'_module.csv'),row.names = F)
      
      summary0 <- c(tissue,cutHeight, mergeCutHeight, deepSplit, paste(table(modulec$module),collapse  ='_') )
      summary <- rbind(summary0,summary)
      
        }}
  }
}


summary_m <- as.data.frame(summary)
colnames(summary_m) <- c('tissue','cutHeight', 'mergeCutHeight', 'deepSplit','module_size')
summary_m <- summary_m[order(summary_m$mergeCutHeight,summary_m$deepSplit),]

write.csv(summary_m, file='./exp_gene_100_summary_module_size.csv')




