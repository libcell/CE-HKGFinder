
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
### code chunk number 12: Validation of HKGs using DNA microarray datasets.
### ****************************************************************************

### ======================================================================== ###
###--For DNA microarray dataset, GSE76380.-----------------------------------### 

library(GEOquery)

dat <- getGEO("GSE76380", GSEMatrix = TRUE)

dat <- dat[[1]]

anno.tab <- dat@featureData@data; dim(anno.tab)

eset <- exprs(dat); dim(eset) 

table(anno.tab$ID == rownames(eset))

probe.num <- nrow(anno.tab)


all.gene <- c("cyc-1", "tba-1", "atp-3", "mdh-1", "gpd-2", 
              "eif-3.C", "act-1", "cdc-42", "pmp-3", "act-2", 
              "csq-1", "ama-1", "rbd-1", "rps-23", "rps-27", 
              "rps-26", "rps-4", "rps-2", "rps-16", "rps-17", 
              "rpl-24.1", "rpl-15", "rpl-35", "rpl-36", "rpl-33", "rpl-27")

eset.mat <- NULL

for (g in all.gene) {
  
  gene.pos <- NULL
  
  for (i in 1:probe.num) {
    
    gene.hit <- which(strsplit(anno.tab$GENE_SYMBOL[i], " /// ")[[1]] == g)
    
    if (length(gene.hit) != 0) gene.pos <- c(gene.pos, i)
    
  }
  
  # gene.pos
  
  tmp <- eset[gene.pos, ]
  
  # print(tmp)
  
  if (length(gene.pos) > 1) gene.eset <- apply(as.matrix(tmp), 2, mean) else gene.eset <- tmp
  
  print(g)
  
  print(gene.eset) 
  
  if (length(gene.eset) == 0) gene.eset <- 0
  
  print("##################################################################")
  
  eset.mat <- rbind(eset.mat, gene.eset)
  
}

rownames(eset.mat) <- all.gene

boxplot(t(eset.mat))

# SD. 

sort(apply(eset.mat, 1, sd))

# GC. 

sort(apply(eset.mat, 1, ineq::Gini))

### End of GSE76380 =========================================================###

### ======================================================================== ###
###--For DNA microarray dataset, GSE118294.----------------------------------### 

library(GEOquery)

dat <- getGEO("GSE118294", GSEMatrix = TRUE)

dat <- dat[[1]]

anno.tab <- dat@featureData@data; dim(anno.tab)

eset <- exprs(dat); dim(eset) 

table(anno.tab$ID == rownames(eset))

probe.num <- nrow(anno.tab)

all.gene <- c("cyc-1", "tba-1", "atp-3", "mdh-1", "gpd-2", 
              "eif-3.C", "act-1", "cdc-42", "pmp-3", "act-2", 
              "csq-1", "ama-1", "rbd-1", "rps-23", "rps-27", 
              "rps-26", "rps-4", "rps-2", "rps-16", "rps-17", 
              "rpl-24.1", "rpl-15", "rpl-35", "rpl-36", "rpl-33", "rpl-27")

eset.mat <- NULL

for (g in all.gene) {
  
  gene.pos <- NULL
  
  for (i in 1:probe.num) {
    
    gene.hit <- which(strsplit(anno.tab$`Gene Symbol`[i], " /// ")[[1]] == g)
    
    if (length(gene.hit) != 0) gene.pos <- c(gene.pos, i)
    
  }
  
  # gene.pos
  
  tmp <- eset[gene.pos, ]
  
  # print(tmp)
  
  if (length(gene.pos) > 1) gene.eset <- apply(as.matrix(tmp), 2, mean) else gene.eset <- tmp

  print(g)
  
  print(gene.eset)
  
  print("##################################################################")
  
  eset.mat <- rbind(eset.mat, gene.eset)
  
}

rownames(eset.mat) <- all.gene

boxplot(t(eset.mat))

# SD. 

sort(apply(eset.mat, 1, sd))

# GC. 

sort(apply(eset.mat, 1, ineq::Gini))

### End of GSE118294 ======================================================= ###




### ======================================================================== ###
###--For DNA microarray dataset, GSE108968.----------------------------------### 

library(GEOquery)

dat <- getGEO("GSE108968", GSEMatrix = TRUE)

dat <- dat[[1]]

anno.tab <- dat@featureData@data; dim(anno.tab)

eset <- exprs(dat); dim(eset) 

table(anno.tab$ID == rownames(eset))

probe.num <- nrow(anno.tab)


all.gene <- c("cyc-1", "tba-1", "atp-3", "mdh-1", "gpd-2", 
              "eif-3.C", "act-1", "cdc-42", "pmp-3", "act-2", 
              "csq-1", "ama-1", "rbd-1", "rps-23", "rps-27", 
              "rps-26", "rps-4", "rps-2", "rps-16", "rps-17", 
              "rpl-24.1", "rpl-15", "rpl-35", "rpl-36", "rpl-33", "rpl-27")

eset.mat <- NULL

for (g in all.gene) {
  
  gene.pos <- NULL
  
  for (i in 1:probe.num) {
    
    gene.hit <- which(strsplit(anno.tab$GENE_SYMBOL[i], " /// ")[[1]] == g)
    
    if (length(gene.hit) != 0) gene.pos <- c(gene.pos, i)
    
  }
  
  # gene.pos
  
  tmp <- eset[gene.pos, ]
  
  # print(tmp)
  
  if (length(gene.pos) > 1) gene.eset <- apply(as.matrix(tmp), 2, mean) else gene.eset <- tmp
  
  print(g)
  
  print(gene.eset) 
  
  if (length(gene.eset) == 0) gene.eset <- 0
  
  print("##################################################################")
  
  eset.mat <- rbind(eset.mat, gene.eset)
  
}

rownames(eset.mat) <- all.gene

boxplot(t(eset.mat))

# SD. 

sort(apply(eset.mat, 1, sd))

# GC. 

sort(apply(eset.mat, 1, ineq::Gini))

### End of GSE108968. =======================================================###


sort(apply(eset.mat, 1, sd))

library(ineq)

sort(apply(eset.mat, 1, Gini))


str(dat)

class(dat)

dim(exprs(dat))

mdat<-getDee2Metadata("celegans")

mdat[which(mdat$GSE_accession %in% "GSE94879"), ]


BiocManager::install("scater")




data("sc_example_counts")
data("sc_example_cell_info")
example_sce <- SingleCellExperiment(
  list(counts = sc_example_counts), 
  colData = sc_example_cell_info)

cpm(example_sce) <- calculateCPM(example_sce, use_size_factors = FALSE)

y <- matrix(rnbinom(20,size=1,mu=10),5,4)
cpm(y)

d <- DGEList(counts=y, lib.size=1001:1004)
cpm(d)
cpm(d,log=TRUE)

d$genes <- data.frame(Length=c(1000,2000,500,1500,3000))
rpkm(d)

cpmByGroup(d, group=c(1,1,2,2))

rpkmByGroup(d, group=c(1,1,2,2))

library(SingleCellExperiment)

