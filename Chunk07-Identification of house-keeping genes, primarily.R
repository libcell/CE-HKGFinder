
################################################################################
#    &&&....&&&    % Project: Identification of HKG candidates in C. elegans   #
#  &&&&&&..&&&&&&  % Author: Bo Li, Jingxin Tao, Youjin Hao                    #
#  &&&&&&&&&&&&&&  % Date: Dec. 24th, 2019                                     #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
#        &         %                                                           #
################################################################################

###==========================================================================### 
### ****************************************************************************
### code chunk number 07: Identification of house-keeping genes, SGLs, primarily.
### ****************************************************************************

setwd("D:\\03-重要数据备份\\26-Tao's HKG results")

exprs.mat2 <- get(load("HKs_149ds.RData")) # the new data. 

l <- length(exprs.mat2) 

#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
## Adding a previouse filtering for all genes, according to the CV.
# Loading the all packages needed in this step. 

library(RankAggreg)

outlier.up <- 2
outlier.down <- 0.5

top5000list <- list()

second.list <- list()

for (ds in 1:l) {
  
  data_nor.method <- exprs.mat2[[ds]]
  
  #----------------------------------------------------------------------------#
  
  hk.ds.each <- list()
  
  top5000 <- NULL
  
  for (d in 1:length(data_nor.method)) {
    
    ### Three rules: 
    # 1) expression observed in all tissues and various conditions. 
    
    eset_d <- data_nor.method[[d]]
    
    # all_genes <- rownames(eset_d)
    
    zero.detect <- function(x) {
      
      res <-  any(x == 0)
      
      return(res)
      
    }
    
    pos_0 <- apply(eset_d, 1, zero.detect)
    
    eset_d <- eset_d[!pos_0, ]
    
    # 2) low variance over tissues: standard-deviation [log 2 (RPKM)]<1. With sd-based ranking. 
    
    # sd_seq <- apply(eset_d, 1, sd)
    # hk.sd <- eset_d[sd_seq < sd_cutoff, ]
    
    a <- apply(eset_d, 1, sd)
    
    hk.sd.ranked <- eset_d[order(a), ]
    
    # 3) no exceptional expression in any single tissue; that is, no log-expression value differed from the averaged log 2 (expression value) by two (four-fold) or more. 
    
    outlier.det <- function(x) {
      
      res <- (max(x) / mean(x) > outlier.up | min(x) / mean(x) < outlier.down)
      
      return(res)
      
    }
    
    pos_outlier <- apply(hk.sd.ranked, 1, outlier.det)
    
    # table(pos_outlier)
    
    hk.sd.ranked <- hk.sd.ranked[!pos_outlier, ]
    
    ranked.hk <- rownames(hk.sd.ranked)
    
    hk.ds.each[[d]] <- ranked.hk
    
    top5000 <- rbind(top5000, ranked.hk[1:5000])
    
  }
  
  # dim(top5000)
  
  # rownames(top5000) <- names(hk.ds.each) <- names(rma)
  
  rownames(top5000) <- names(hk.ds.each) <- names(data_nor.method)
  
  colnames(top5000) <- paste("Gene", 1:5000, sep = "-")
  
  # Removing the wrong datasets from 149 datasets !!!. 
  
  DS.tobe.remove <- c("GSE81854", "GSE95603", "GSE41486", "GSE7354")
  
  pos.remove <- match(DS.tobe.remove, rownames(top5000))
  
  top5000.e <- top5000[-pos.remove, ]
  
  # top5000list <- rbind(top5000list, rownames(top5000.e))
  
  # dim(top5000.e) 
  
  top5000list[[ds]] <- top5000.e
  
  # ranking
  
  # more complex example (to get a better solution, increase maxIter)
  
  Topgenelist <- RankAggreg(top5000.e, 30, method = "CE", distance = "Spearman", 
                            maxIter = 1000)
  plot(Topgenelist)
  
  # Topgenelist
  
  toplist10 <- Topgenelist$top.list
  
  #----------------------------------------------------------------------------#
  
  second.list[[ds]] <- toplist10
  
}

names(exprs.mat2)[[4]] <- "Li-Wong"

names(top5000list) <- gsub("_exprs", "", names(exprs.mat2))

names(second.list) <- gsub("_exprs", "", names(exprs.mat2))

save(top5000list, file = "top5000list.RData")

library(openxlsx)

write.xlsx(top5000list, "S4-top5000.xlsx", row.names = TRUE)

###==========================================================================###

# End of this line. 

