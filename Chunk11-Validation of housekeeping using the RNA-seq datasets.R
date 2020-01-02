
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
### code chunk number 11: Validation of HKGs using RNA-sequencing datasets.
### ****************************************************************************

###==========================================================================### 

### Common variation. 
genes_26 <- c("cyc-1", "tba-1", "atp-3", "mdh-1", "gpd-2", 
              "eif-3.C", "act-1", "cdc-42", "pmp-3", "act-2", 
              "csq-1", "ama-1", "rbd-1", "rps-23", "rps-27", 
              "rps-26", "rps-4", "rps-2", "rps-16", "rps-17", 
              "rpl-24.1", "rpl-15", "rpl-35", "rpl-36", "rpl-33","rpl-27")

col_26 <- rep(c("#FF4500", "#00A600FF"), times = c(13, 13))

op <- par(mfcol = c(3, 2)) # width:length = 12:8 

setwd("J:/HKG_test/")


###----------------------- (1) Microarray: GSE11894 -------------------------###
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

for (g in genes_26) {
  
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
  
  print("######################################################333############")
  
  eset.mat <- rbind(eset.mat, gene.eset)
  
}

rownames(eset.mat) <- genes_26

# boxplot(t(eset.mat))

sd.seq <- apply(eset.mat, 1, sd)

x <- t(eset.mat[sd.seq != 0, ])

y <- reshape2::melt(x)

vioplot(value ~ Var2, data = y, col = col_26[sd.seq != 0], las = 3, 
        main = "GSE118294", xlab = NULL, ylab = "expression level") 

# SD. 

sort(sd.seq)

# GC. 

sort(apply(eset.mat, 1, ineq::Gini))

### End of GSE118294 ======================================================= ###
###--------------------------------------------------------------------------###


###----------------------- (2) Microarray: GSE108968 ------------------------###
###--For DNA microarray dataset, GSE108968.----------------------------------### 

library(limma)

#. dat <- read.maimages(files = dir("GSE108968_RAW"),
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

dat <- get(load("GSE108968.RData"))

dat2 <- backgroundCorrect(dat, "normexp", offset = 50)

eset <- normalizeBetweenArrays(dat2$R, method = "quantile")

eset2 <- log(eset, 2)

rownames(eset2) <- dat$genes$ProbeName

eset2 <- as.data.frame(eset2)

# write.csv(eset2,"GSE108968_expr.csv")

########################################################

gse108968 <- getGEO("GSE108968", GSEMatrix = TRUE)

str(gse108968)

anno.tab <- gse108968[[1]]@featureData@data; dim(anno.tab)

rownames(eset2) <- NULL

eset <- eset2; dim(eset) 

# table(anno.tab$NAME == rownames(eset))

nr <- nrow(anno.tab)

eset.mat <- NULL

for (g in genes_26) {
  
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

rownames(eset.mat) <- genes_26

sd.seq <- apply(eset.mat, 1, sd)

x <- t(eset.mat[sd.seq != 0, ]); 
y <- reshape2::melt(x)

vioplot(value ~ Var2, data = y, col = col_26[sd.seq != 0], las = 3, 
        main = "GSE108968", xlab = NULL, ylab = "expression level") 

# SD. 

sort(sd.seq)

# GC. 

sort(apply(eset.mat, 1, ineq::Gini))

### End of GSE108968. =======================================================###
###--------------------------------------------------------------------------###


###----------------------- (3) Microarray: GSE76380 -------------------------###
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

########################################################

gse76380 <- getGEO("GSE76380", GSEMatrix = TRUE)

str(gse76380)

anno.tab <- gse76380[[1]]@featureData@data; dim(anno.tab)

rownames(eset2) <- NULL

eset <- eset2; dim(eset) 

# table(anno.tab$NAME == rownames(eset))

nr <- nrow(anno.tab)

eset.mat <- NULL

for (g in genes_26) {
  
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

rownames(eset.mat) <- genes_26

sd.seq <- apply(eset.mat, 1, sd)

x <- t(eset.mat[sd.seq != 0, ]); 
y <- reshape2::melt(x)

vioplot(value ~ Var2, data = y, col = col_26[sd.seq != 0], las = 3, 
        main = "GSE76380", xlab = NULL, ylab = "expression level") 

# SD. 

sort(sd.seq)

# GC. 

sort(apply(eset.mat, 1, ineq::Gini))

### End of GSE76380. ========================================================###
###--------------------------------------------------------------------------###


###----------------------- (4) Microarray: GSE76380 -------------------------###
###--For RNA-sequencing dataset, GSE63528, GSE60755 and GSE98919-------------### 

### Step-01. Comparison of gene expression levels using RNA-sequencing datasets

# devtools::install_github("markziemann/dee2/getDEE2")

library('getDEE2')

mdat <- getDee2Metadata("celegans")

expre_gse <- c("GSE63528", 
               "GSE60755", 
               "GSE98919")

for (i in expre_gse) {
 
 mdat1 <- mdat[which(mdat$GSE_accession %in% i),] 

 SRRlist <- as.vector(mdat1$SRR_accession)

 dat <- getDEE2("celegans", SRRlist)
 
 gene.count <- dat$GeneCounts
 
 gene.info <- dat$GeneInfo
  
 expr <- cbind(gene.info, gene.count)
  
 write.csv(expr, paste0(i, "_count.csv")) 
 
 }



library(edgeR)

library(openxlsx)

fileName <- dir()[grep("_count.csv", dir())]

ask.mean <- function(x) {mean(x, na.rm = TRUE)}

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
    
    if (!is.null(nrow(tmp))) tmp1 <- apply(tmp, 2, ask.mean) else tmp1 <- tmp
    
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
  
  gse.nam <- strsplit(f, "_")[[1]][1]
  
  xlsx.file <- paste(gse.nam, "xlsx", sep = ".")
  
  write.xlsx(res, file = xlsx.file, col.names = TRUE, row.names = TRUE)
  
  Sys.sleep(20)
  
  # Boxplot 
  # print(mat.26[tmp1 != 0, ])
  print(tmp1)
  
  # boxplot(t(mat.26[tmp1 != 0, ]), col = col_26[tmp1 != 0], las = 3)
  
  x <- t(mat.26[tmp1 != 0, ]); 
  y <- reshape2::melt(x)

  vioplot(value ~ Var2, data = y, col = col_26[tmp1 != 0], las = 3, 
          main = gse.nam, xlab = NULL, ylab = "expression level") 
  
  print("############################# End ####################################")

  }

par(op)

### End of the chunk 11. 
###==========================================================================### 
