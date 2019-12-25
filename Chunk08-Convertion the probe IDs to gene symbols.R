
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
### code chunk number 08: Convertion the probe IDs to gene symbols in top 5000.
### ****************************************************************************

top5000list <- get(load("J:/top5000list.RData"))

gpl200anno <- get(load("J:/gpl200anno.RData"))

for (n in 1:length(top5000list)) {
  
  for (i in 1:ncol(top5000list[[n]])) {
    
    probe.len <- length(top5000list[[n]][, i]) 
    
    at.pos <- grep("_at", top5000list[[n]][, i])
    
    at2sym <- NULL
    
    for (g in at.pos) {
      
      tp <- gpl200anno$`Gene Symbol`[match(top5000list[[n]][g, i], gpl200anno$ID)]
      
      print(g)
      
      print(tp)
      
      tmp <- strsplit(tp, " /// ")[[1]][1]
      
      if (length(grep("WBGene", tmp)) != 0) {
        
        x <- strsplit(tp, " /// ")[[1]][2] 
        
        tmp <- x
        
      }
      
      print(tmp)
      
      at2sym <- c(at2sym, tmp)
      
    }
    
    top5000list[[n]][at.pos, i] <- at2sym
    
  }
  
}

save(top5000list, file = "top5000list_symbol.RData")

library(openxlsx)

write.xlsx(top5000list, "S4-top5000_gene_symbol.xlsx", row.names = TRUE)

# End of this line. 
