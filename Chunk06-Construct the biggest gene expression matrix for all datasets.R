
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
### code chunk number 06: Construct the gene expression matrix.
### ****************************************************************************

###- RMA Expression Matrix --------------------------------------------------###

len.DS <- length(data_rma)

RMA_exprs <- list()

for (i in 1:len.DS) {
  
  ### S1. One2one gene expression matrix. 
  
  one2one.matrix <- data_rma[[i]][one2one.probe$PROBEID, ]
  
  rownames(one2one.matrix) == one2one.probe$PROBEID
  
  ### S2. More2one gene expression matrix. 
  
  more2one.matrix <- NULL
  
  for (s in 1:length(more2one)) {
    
    pos <- which(one2any.probe$SYMBOL == more2one[s])
    
    probe.pos <- one2any.probe$PROBEID[pos]
    
    tmp <- data_rma[[i]][probe.pos, ]
    
    exp <- apply(tmp, 2, median)
    
    more2one.matrix <- rbind(more2one.matrix, exp)
    
  }
  
  rownames(more2one.matrix) <- more2one
  
  dim(more2one.matrix)
  
  ### S3. Whole gene expression matrix. 
  
  whole.matrix <- rbind(one2one.matrix, more2one.matrix)
  
  ### S4. Gene expression matrix.
  
  RMA_exprs[[i]] <- whole.matrix
  
}



###- MAS Expression Matrix --------------------------------------------------###

len.DS <- length(data_mas)

MAS_exprs <- list()

for (i in 1:len.DS) {
  
  ### S1. One2one gene expression matrix. 
  
  one2one.matrix <- data_mas[[i]][one2one.probe$PROBEID, ]
  
  rownames(one2one.matrix) == one2one.probe$PROBEID
  
  ### S2. More2one gene expression matrix. 
  
  more2one.matrix <- NULL
  
  for (s in 1:length(more2one)) {
    
    pos <- which(one2any.probe$SYMBOL == more2one[s])
    
    probe.pos <- one2any.probe$PROBEID[pos]
    
    tmp <- data_mas[[i]][probe.pos, ]
    
    exp <- apply(tmp, 2, median)
    
    more2one.matrix <- rbind(more2one.matrix, exp)
    
  }
  
  rownames(more2one.matrix) <- more2one
  
  dim(more2one.matrix)
  
  ### S3. Whole gene expression matrix. 
  
  whole.matrix <- rbind(one2one.matrix, more2one.matrix)
  
  ### S4. Gene expression matrix.
  
  MAS_exprs[[i]] <- whole.matrix
  
}



###- GCRMA Expression Matrix ------------------------------------------------###

len.DS <- length(data_gcrma)

GCRMA_exprs <- list()

for (i in 1:len.DS) {
  
  ### S1. One2one gene expression matrix. 
  
  one2one.matrix <- data_gcrma[[i]][one2one.probe$PROBEID, ]
  
  rownames(one2one.matrix) == one2one.probe$PROBEID
  
  ### S2. More2one gene expression matrix. 
  
  more2one.matrix <- NULL
  
  for (s in 1:length(more2one)) {
    
    pos <- which(one2any.probe$SYMBOL == more2one[s])
    
    probe.pos <- one2any.probe$PROBEID[pos]
    
    tmp <- data_gcrma[[i]][probe.pos, ]
    
    exp <- apply(tmp, 2, median)
    
    more2one.matrix <- rbind(more2one.matrix, exp)
    
  }
  
  rownames(more2one.matrix) <- more2one
  
  dim(more2one.matrix)
  
  ### S3. Whole gene expression matrix. 
  
  whole.matrix <- rbind(one2one.matrix, more2one.matrix)
  
  ### S4. Gene expression matrix.
  
  GCRMA_exprs[[i]] <- whole.matrix
  
}


###- DCHIP Expression Matrix ------------------------------------------------###

len.DS <- length(data_dchip)

DCHIP_exprs <- list()

for (i in 1:len.DS) {
  
  ### S1. One2one gene expression matrix. 
  
  one2one.matrix <- data_dchip[[i]][one2one.probe$PROBEID, ]
  
  rownames(one2one.matrix) == one2one.probe$PROBEID
  
  ### S2. More2one gene expression matrix. 
  
  more2one.matrix <- NULL
  
  for (s in 1:length(more2one)) {
    
    pos <- which(one2any.probe$SYMBOL == more2one[s])
    
    probe.pos <- one2any.probe$PROBEID[pos]
    
    tmp <- data_dchip[[i]][probe.pos, ]
    
    exp <- apply(tmp, 2, median)
    
    more2one.matrix <- rbind(more2one.matrix, exp)
    
  }
  
  rownames(more2one.matrix) <- more2one
  
  dim(more2one.matrix)
  
  ### S3. Whole gene expression matrix. 
  
  whole.matrix <- rbind(one2one.matrix, more2one.matrix)
  
  ### S4. Gene expression matrix.
  
  DCHIP_exprs[[i]] <- whole.matrix
  
}

###- PLIER Expression Matrix ------------------------------------------------###

len.DS <- length(data_plier)

PLIER_exprs <- list()

for (i in 1:len.DS) {
  
  ### S1. One2one gene expression matrix. 
  
  one2one.matrix <- data_plier[[i]][one2one.probe$PROBEID, ]
  
  rownames(one2one.matrix) == one2one.probe$PROBEID
  
  ### S2. More2one gene expression matrix. 
  
  more2one.matrix <- NULL
  
  for (s in 1:length(more2one)) {
    
    pos <- which(one2any.probe$SYMBOL == more2one[s])
    
    probe.pos <- one2any.probe$PROBEID[pos]
    
    tmp <- data_plier[[i]][probe.pos, ]
    
    exp <- apply(tmp, 2, median)
    
    more2one.matrix <- rbind(more2one.matrix, exp)
    
  }
  
  rownames(more2one.matrix) <- more2one
  
  dim(more2one.matrix)
  
  ### S3. Whole gene expression matrix. 
  
  whole.matrix <- rbind(one2one.matrix, more2one.matrix)
  
  ### S4. Gene expression matrix.
  
  PLIER_exprs[[i]] <- whole.matrix
  
}


###- VSN Expression Matrix ------------------------------------------------###

len.DS <- length(data_vsn)

VSN_exprs <- list()

for (i in 1:len.DS) {
  
  ### S1. One2one gene expression matrix. 
  
  one2one.matrix <- data_vsn[[i]][one2one.probe$PROBEID, ]
  
  rownames(one2one.matrix) == one2one.probe$PROBEID
  
  ### S2. More2one gene expression matrix. 
  
  more2one.matrix <- NULL
  
  for (s in 1:length(more2one)) {
    
    pos <- which(one2any.probe$SYMBOL == more2one[s])
    
    probe.pos <- one2any.probe$PROBEID[pos]
    
    tmp <- data_vsn[[i]][probe.pos, ]
    
    exp <- apply(tmp, 2, median)
    
    more2one.matrix <- rbind(more2one.matrix, exp)
    
  }
  
  rownames(more2one.matrix) <- more2one
  
  dim(more2one.matrix)
  
  ### S3. Whole gene expression matrix. 
  
  whole.matrix <- rbind(one2one.matrix, more2one.matrix)
  
  ### S4. Gene expression matrix.
  
  VSN_exprs[[i]] <- whole.matrix
  
}

### Checking the dimension of each gene expression

dim(RMA_exprs[[1]])
dim(MAS_exprs[[1]])
dim(GCRMA_exprs[[1]])
dim(DCHIP_exprs[[1]])
dim(PLIER_exprs[[1]])
dim(VSN_exprs[[1]])

### renaming the six gene expression lists.
#
names(RMA_exprs) <- names(MAS_exprs) <- names(GCRMA_exprs) <- names(DCHIP_exprs) <- names(PLIER_exprs) <- names(VSN_exprs) <- Series.name

### Merging six gene expression lists into an large list consisting of all gene expression matrices.  

exprs.mat2 <- list(RMA_exprs = RMA_exprs, 
                   MAS_exprs = MAS_exprs, 
                   GCRMA_exprs = GCRMA_exprs, 
                   DCHIP_exprs = DCHIP_exprs, 
                   PLIER_exprs = PLIER_exprs,
                   VSN_exprs = VSN_exprs)

### Cleaning the menmory of computer. 

rm(RMA_exprs)
rm(MAS_exprs)
rm(GCRMA_exprs)
rm(DCHIP_exprs)
rm(PLIER_exprs)
rm(VSN_exprs)

### Checking the data integrity

rma <- exprs.mat2$RMA_exprs
mas <- exprs.mat2$MAS_exprs
gcrma <- exprs.mat2$GCRMA_exprs
dchip <- exprs.mat2$DCHIP_exprs
plier <- exprs.mat2$PLIER_exprs
vsn <- exprs.mat2$VSN_exprs

names(rma) == names(mas)

setwd("/Users/libo/Desktop/")

save(exprs.mat2, file = "HKs_149ds.RData")

###--------------------------------------------------------------------------###

