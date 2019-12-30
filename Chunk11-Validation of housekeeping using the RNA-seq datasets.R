
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
### code chunk number 11: Validation of HKGs using RNA-sequencing datasets.
### ****************************************************************************

### Step-01. Comparison of gene expression levels using RNA-sequencing datasets

setwd("J:/")

dir.create("HKG_test")

setwd("HKG_test")

# devtools::install_github("markziemann/dee2/getDEE2")

library('getDEE2')

mdat <- getDee2Metadata("celegans")

expre_gse <- c(# "GSE52861",
               # "GSE63528",
               # "GSE94879",
               # "GSE87528",
               # "GSE102680",
               # "GSE124049", 
               # "GSE87078", 
  
               # "GSE111364", 
               # "GSE49043", 
               # "GSE60063", 
               # "GSE108263",
               "GSE49043",
               "GSE98915"
               
               # "GSE50548", 
               # "GSE60755"
               )

for (i in expre_gse) {
 
 mdat1 <- mdat[which(mdat$GSE_accession %in% i),] 

 SRRlist <- as.vector(mdat1$SRR_accession)

 dat <- getDEE2("celegans", SRRlist)
 
 gene.count <- dat$GeneCounts
 
 gene.info <- dat$GeneInfo
  
 expr <- cbind(gene.info, gene.info)
  
 write.csv(expr, paste0(i, "_count.csv")) 
 
 }


library(edgeR)

library(openxlsx)

fileName <- dir()[grep("_count.csv", dir())]

aaa <- function(x) {mean(x, na.rm = TRUE)}

genes_26 <- c("cyc-1", "tba-1", "atp-3", "mdh-1", "gpd-2", 
              "eif-3.C", "act-1", "cdc-42", "pmp-3", "act-2", 
              "csq-1", "ama-1", "rbd-1", "rps-23", "rps-27", 
              "rps-26", "rps-4", "rps-2", "rps-16", "rps-17", 
              "rpl-24.1", "rpl-15", "rpl-35", "rpl-36", "rpl-33","rpl-27")

for (f in fileName) { 
  
  data <- read.csv(f) 
  
  data_norm <- cpm.default(data[, -(1:6)]) 
  
  sym <- as.character(data$GeneSymbol) 
  
  rownames(data_norm) <- sym 
  
  cpm_data <- data_norm 
  
  mat.26 <- NULL 
  
  for (g in genes_26) { 
    
    pos <- which(rownames(cpm_data) == g)
    
    tmp <- cpm_data[pos, ]
    
    if (!is.null(nrow(tmp))) tmp1 <- apply(tmp, 2, aaa) else tmp1 <- tmp
    
    mat.26 <- rbind(mat.26, tmp1)  
    
    } 
  
  rownames(mat.26) <- genes_26 
  
  mat.26 <- log(mat.26 + 1, 2) 
  
  #Gini rank  
  
  tmp <- apply(mat.26, 1, ineq::Gini) 
  
  GC <- sort(tmp) 
  
  GC <- as.matrix(GC)
  
  print("############################ Start ##################################")
  
  print(f)
  
  #print(as.data.frame(GC[1:10])) 
  
  # write.csv(GC[1:10], paste(f, "csv", sep = "."))
  
  #SD rank 
  
  tmp1 <- apply(mat.26, 1, sd) 
  
  SD <- sort(tmp1) 
  
  SD <- as.matrix(SD)
  
  #print(as.data.frame(SD[1:10])) 
  
  res <- list(SD = SD, GC = GC)
  
  write.xlsx(res, file = paste(f, "xlsx", sep = "."), col.names = TRUE, row.names = TRUE)
  
  Sys.sleep(3)
  
  print("############################# End ####################################")

  }

### End of Step-01.

