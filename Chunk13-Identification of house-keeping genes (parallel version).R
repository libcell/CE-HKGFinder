
################################################################################
#    &&&....&&&    % Project: Identification of HKG candidates in C. elegans   #
#  &&&&&&..&&&&&&  % Author: Bo Li, Jingxin Tao, Youjin Hao                    #
#  &&&&&&&&&&&&&&  % Date: Dec. 24th, 2019                                     #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
#        &         %                                                           #
################################################################################

### ****************************************************************************
### code chunk number 13: Paralle computing for RankAggregation, about chunk 07.
### ****************************************************************************

# RMA

# rma <- get(load("rma_top5000.RData"))
# exprs.mat2

# setwd("/Users/libo/Desktop/")

exprs.mat2 <- get(load("HKs_149ds.RData"))

library(parallel)

library(RankAggreg)

hk.finder <- function(x, top = 5000, num = 100) {
  
  outlier.up <- 2
  outlier.down <- 0.5
  
  hk.ds.each <- list()
  
  top5000 <- NULL
  
  for (d in 1:length(x)) {
    
    ### Three rules: 
    # 1) expression observed in all tissues. 
    
    eset_d <- x[[d]]
    
    all_genes <- rownames(eset_d)
    
    pos_0 <- NULL
    
    for (i in all_genes) {
      
      res <-  any(eset_d[i, ] == 0)
      
      pos_0 <- c(pos_0, res)
      
    }
    
    # 2) low variance over tissues: standard-deviation [log 2 (RPKM)]<1. With sd-based ranking. 
    
    sd_seq <- apply(eset_d, 1, sd)
    
    hk.sd <- eset_d
    
    a <- apply(hk.sd, 1, sd)
    
    hk.sd.ranked <- hk.sd[order(a), ]
    
    # 3) no exceptional expression in any single tissue; that is, no log-expression value differed from the averaged log 2 (RPKM) by two (four-fold) or more. 
    
    outlier.det <- function(x) {
      
      res <- (max(x) / mean(x) > outlier.up | min(x) / mean(x) < outlier.down)
      
      return(res)
      
    }
    
    pos_outlier <- apply(hk.sd.ranked, 1, outlier.det)
    
    # table(pos_outlier)
    
    hk.sd.ranked <- hk.sd.ranked[!pos_outlier, ]
    
    ranked.hk <- rownames(hk.sd.ranked)
    
    hk.ds.each[[d]] <- ranked.hk
    
    top5000 <- rbind(top5000, ranked.hk[1:top])
    
  }
  
  # dim(top5000)
  
  rownames(top5000) <- names(hk.ds.each) <- names(data_nor.method)
  
  colnames(top5000) <- paste("Gene", 1:top, sep = "-")
  
  # Removing the wrong datasets. 
  
  DS.tobe.remove <- c("GSE81854", "GSE95603", "GSE41486", "GSE7354")
  
  pos.remove <- match(DS.tobe.remove, rownames(top5000))
  
  top5000.e <- top5000[-pos.remove, ]
  
  # dim(top5000.e)
  
  # ranking
  
  # more complex example (to get a better solution, increase maxIter)
  
  Topgenelist <- RankAggreg::RankAggreg(top5000.e, num, method = "CE", distance = "Spearman", 
                                        maxIter = 1000)
  # plot(Topgenelist)
  # Topgenelist
  
  toplist10 <- Topgenelist$top.list
  
  return(toplist10)
  
}

hk.finder(exprs.mat2, top = 100, num = 10)

# no_cores <- detectCores() - 1

# cl <- makeCluster(no_cores)

# result.rank <- parLapply(cl, exprs.mat2, hk.finder)

# stopCluster(cl)

# determine a function for identifying the hkgs. 

# End of this line. 
