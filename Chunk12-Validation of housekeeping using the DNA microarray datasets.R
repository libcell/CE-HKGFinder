
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
### code chunk number 12: Validation of HKGs using DNA microarray datasets.
### ****************************************************************************
###==========================================================================### 

### Common variation. 
genes_26 <- c("cyc-1", "tba-1", "atp-3", "mdh-1", "gpd-2", 
              "eif-3.C", "act-1", "cdc-42", "pmp-3", "act-2", 
              "csq-1", "ama-1", "rbd-1", "rps-23", "rps-27", 
              "rps-26", "rps-4", "rps-2", "rps-16", "rps-17", 
              "rpl-24.1", "rpl-15", "rpl-35", "rpl-36", "rpl-33","rpl-27")

###--For DNA microarray dataset, GSE118294.----------------------------------### 

library(affy)

setwd("GSE118294_RAW/")

dat <- ReadAffy()

eset <- rma(dat)

es <- exprs(eset)

setwd("..")

library(GEOquery)

dat <- getGEO("GSE118294", GSEMatrix = TRUE)

anno.tab <- dat[[1]]@featureData@data; dim(anno.tab)

eset <- es; dim(eset) 

table(anno.tab$ID == rownames(eset))

probe.num <- nrow(anno.tab)

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


### ======================================================================== ###
###--For DNA microarray dataset, GSE76380.-----------------------------------### 

library(limma)

#. dat <- read.maimages(files = dir("GSE76380_RAW"),
#.                      source="agilent", 
#.                      columns = list(G = "gMedianSignal", 
#.                                     Gb = "gBGMedianSignal", 
#.                                     R = "gMedianSignal", 
#.                                     Rb = "gBGMedianSignal"), 
#.                      annotation = c("Row", 
#.                                     "Col",
#.                                     "FeatureNum",
#.                                     "ControlType",
#.                                     "ProbeName",
#.                                     "GeneName",
#.                                     "SystematicName"))

dat <- get(load("GSE76380.RData"))

dat2 <- backgroundCorrect(dat, "normexp", offset = 50)

eset <- normalizeBetweenArrays(dat2$R, method = "quantile")

eset2 <- log(eset, 2)

rownames(eset2) <- dat$genes$ProbeName

eset2 <- as.data.frame(eset2)

# setwd("..")

# write.csv(eset2,"GSE108968_expr.csv")

dat.tmp <- eset2

########################################################

gse76380 <- getGEO("GSE76380", GSEMatrix = TRUE)

str(gse76380)

anno.tab <- gse76380[[1]]@featureData@data; dim(anno.tab)

rownames(eset2) <- NULL

eset <- eset2; dim(eset) 

# table(anno.tab$NAME == rownames(eset))

nr <- nrow(anno.tab)

eset.mat <- NULL

for (g in all.gene) {
  
  gene.pos <- NULL
  
  for (i in 1:nr) {
    
    gene.hit <- which(strsplit(anno.tab$GENE_SYMBOL[i], " /// ")[[1]] == g)
    
    if (length(gene.hit) != 0) gene.pos <- c(gene.pos, i)
    
  }
  
  # gene.pos
  
  tmp <- eset[gene.pos, ]
  
  # print(tmp)
  
  if (length(gene.pos) > 1) gene.eset <- apply(as.matrix(tmp), 2, mean) else 
    
    if (length(gene.pos) == 0) gene.eset <- 0 
  
  else gene.eset <- tmp
  
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

### End of the chunk 12. ====================================================###
### ======================================================================== ###
