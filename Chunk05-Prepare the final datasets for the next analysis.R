
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
### code chunk number 05: Prepare the final datasets for the next analysis.
### ****************************************************************************

### Step-01. Load the saved datasets with RData format.

setwd("00-Preprocessed data/old/")

setwd("../new/")

data_rma <- get(load("rma_exprs.RData"))
data_mas <- get(load("mas5log2_exprs.RData"))
data_gcrma <- get(load("gcrma_exprs.RData"))
data_dchip <- get(load("dchip_log2_exprs.RData"))
data_plier <- get(load("plier_exprs.RData"))
data_vsn <- get(load("vsn_exprs.RData"))

Series.name <- get(load("Series.name.RData"))

# save(Series.name, file = "Series.name.RData")

error.gse <- c("GSE81854", "GSE95603", "GSE41486", "GSE7354")

### Step-02. Load the currently saved datasets with RData format.

names(data_rma) <- Series.name
names(data_mas) <- Series.name
names(data_gcrma) <- Series.name
names(data_dchip) <- Series.name
names(data_plier) <- Series.name
names(data_vsn) <- Series.name

# save(error.gse, file = "error.gse.RData")

error.pos <- match(error.gse, Series.name)

table(data_rma[[123]] == data_rma[["GSE81854"]])

# Prepare the final datasets for next analysis. 

final.DS <- list(
  data_rma = data_rma, 
  data_mas = data_mas, 
  data_gcrma = data_gcrma, 
  data_dchip = data_dchip, 
  data_plier = data_plier, 
  data_vsn = data_vsn
)

save(final.DS, file = "finalDS2.RData")

#. library(qs)
#. qsave(final.DS, file = "finalDS3.RData")
#. help(package = "qs")

### Step-02. Data integrity check. 

for (i in 1:length(data_rma)) {
  
  check.mas <- table((data_rma[[i]] == data_mas[[i]]))
  print(check.mas)
  
  check.gcrma <- table((data_rma[[i]] == data_gcrma[[i]]))
  print(check.gcrma)
  
  check.dchip <- table((data_rma[[i]] == data_dchip[[i]]))
  print(check.dchip)
  
  check.plier <- table((data_rma[[i]] == data_plier[[i]]))
  print(check.plier)
  
  check.vsn <- table((data_rma[[i]] == data_vsn[[i]]))
  print(check.vsn)
  
}

save(data_rma, file = "rma_exprs.RData")
save(data_mas, file = "mas5log2_exprs.RData")
save(data_gcrma, file = "gcrma_exprs.RData")
save(data_dchip, file = "dchip_log2_exprs.RData")
save(data_plier, file = "plier_exprs.RData")
save(data_vsn, file = "vsn_exprs.RData")

### Step-03. Annotation and construction of gene expression matrix. 

library(annotate)
library(celegans.db)

dim(data_rma[[1]])

### S1. Checking the probesets for all subdatasets. 

probe.list <- rownames(data_rma[[1]])

for (m in 1:length(data_rma)) {
  
  check.res <- (rownames(data_rma[[m]]) == probe.list)
  check.res <- (rownames(data_mas[[m]]) == probe.list)
  check.res <- (rownames(data_gcrma[[m]]) == probe.list)
  check.res <- (rownames(data_dchip[[m]]) == probe.list)
  check.res <- (rownames(data_plier[[m]]) == probe.list)
  check.res <- (rownames(data_vsn[[m]]) == probe.list)
  print(check.res)
  
}

table(rownames(data_rma[[1]]) == rownames(data_plier[[141]]))


### Step-1. Parse probesets to gene symbols. 

### S1. get the probeset-gene mapping matrix.

# geneID <- select(celegans.db, keys = probe.list, "ENTREZID", "PROBEID")

gene.smb <- select(celegans.db, keys = probe.list, "SYMBOL", "PROBEID")

length(unique(gene.smb$PROBEID))

nrow(gene.smb)

### S2. remove the rows with NAs.

gene.symbol <- na.omit(gene.smb)

nrow(gene.symbol)

#. table(is.na(gene.symbol))
#. length(unique(gene.symbol$PROBEID))
#. length(table(gene.symbol$PROBEID))

### S3. remain the probe which mapped only one gene. 

# sum(table(table(gene.symbol$PROBEID)))

one2any <- names(table(gene.symbol$PROBEID)[table(gene.symbol$PROBEID) == 1]) # 

# Method-1. 
one2any.probe <- gene.symbol[match(one2any, gene.symbol$PROBEID), ]

# Method-2. 
p.pos <- NULL

for (p in one2any) {
  
  p.pos <- c(p.pos, which(gene.symbol$PROBEID == p))
  
}

p.pos

one2any.probe1 <- gene.symbol[p.pos, ]

table(one2any.probe == one2any.probe1)

nrow(one2any.probe)

# table(is.na(one2any.probe))

length(unique(one2any.probe$PROBEID))


# Investigate the detail information of probesets. 

# gene number in total. 

length(table(one2any.probe$SYMBOL))

# one2one.

one2one <- names(table(one2any.probe$SYMBOL)[table(one2any.probe$SYMBOL) == 1])

# sum(table(one2any.probe$SYMBOL))

# Method-1
one2one.probe <- one2any.probe[match(one2one, one2any.probe$SYMBOL), ]

# Method-2. 
p.pos <- NULL

for (p in one2one) {
  
  p.pos <- c(p.pos, which(one2any.probe$SYMBOL == p))
  
}

p.pos

one2one.probe1 <- one2any.probe[p.pos, ]

table(one2one.probe == one2one.probe1)

nrow(one2one.probe)

# length(match(one2one, one2any.probe$SYMBOL))

nrow(one2one.probe)
length(unique(one2one.probe$PROBEID))
length(unique(one2one.probe$SYMBOL))

# more2one.

# compute the number of genes which count is more than 1. 

sum(table(one2any.probe$SYMBOL) > 1)

# the number of probesets with more2one style.  
sum(table(one2any.probe$SYMBOL)[table(one2any.probe$SYMBOL) > 1]) 

# the gene symbol with more2one style. 
more2one <- names(table(one2any.probe$SYMBOL)[table(one2any.probe$SYMBOL) > 1])

length(more2one)

# End of this line. 

