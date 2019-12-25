
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
### code chunk number 04: Preprocessing and integrating the individual datasets.
### ****************************************************************************

################################################################################

###--------------------------- S1. RMA algorithm ----------------------------###

setwd("GSE datasets about C.elegans (GPL200)")

cel_dir <- dir()[grep("GSE", dir())]

for (i in cel_dir) {
  
  setwd(i)
  
  nd <- nd + 1
  
  dat <- affy::ReadAffy()
  eset <- affy::rma(dat)
  expre_mat <- affy::exprs(eset)
  
  # probe_rank <- rownames(expre_mat)[order(apply(expre_mat, 1, f))]
  # result_list[[nd]] <- probe_rank
  
  # tmp <- apply(expre_mat, 1, f)
  
  result_list[[nd]] <- expre_mat
  
  setwd("..")
  
}

# names(result_list) <- cel_dir

eset_rma <- result_list

save(eset_rma, file = "rma_exprs.RData")

###--------------------------- S2. MAS algorithm ----------------------------###

setwd("GSE datasets about C.elegans (GPL200)")

result_list <- list()

nd <- 0

cel_dir <- dir()[grep("GSE", dir())]

for (i in cel_dir) {
  
  setwd(i)
  
  nd <- nd + 1
  
  dat <- affy::ReadAffy()
  eset <- affy::mas5(dat)
  expre_mat <- log2(affy::exprs(eset))
  
  result_list[[nd]] <- expre_mat
  
  print(nd)
  
  setwd("..")
  
}

result_list

# names(result_list) <- cel_dir

eset_mas <- result_list

save(eset_mas, file = "mas5log2_exprs.RData")

###-------------------------- S3. dChip algorithm ---------------------------###

setwd("GSE datasets about C.elegans (GPL200)")

result_list <- list()

nd <- 0

cel_dir <- dir()[grep("GSE", dir())]

for (i in cel_dir) {
  
  setwd(i)
  
  nd <- nd + 1
  
  if (nd > 120) {
    
    dat <- affy::ReadAffy()
    
    # eset <- affy::mas5(dat)
    
    eset <- affy::expresso(dat, 
                           normalize.method = "invariantset", 
                           bg.correct = FALSE, 
                           pmcorrect.method = "pmonly", 
                           summary.method = "liwong")
    
    
    expre_mat <- log2(affy::exprs(eset))
    
    result_list[[nd]] <- expre_mat
    
    print(nd)
    
    
  }
  
  setwd("..")
  
}

result_list

# names(result_list) <- cel_dir

eset_dchip <- result_list

save(eset_dchip, file = "dchip_log2_exprs2.RData")


###-------------------------- S4. gcRMA algorithm ---------------------------###

setwd("GSE datasets about C.elegans (GPL200)")

result_list <- list()

nd <- 0

cel_dir <- dir()[grep("GSE", dir())]

for (i in cel_dir) {
  
  setwd(i)
  
  nd <- nd + 1
  
  dat <- affy::ReadAffy()
  
  eset <- gcrma::gcrma(dat)
  
  expre_mat <- exprs(eset)
  
  result_list[[nd]] <- expre_mat
  
  setwd("..")
  
}

result_list

# names(result_list) <- cel_dir

eset_gcrma <- result_list

save(eset_gcrma, file = "gcrma_exprs.RData")


###-------------------------- S5. PLIER algorithm ---------------------------###

setwd("GSE datasets about C.elegans (GPL200)")

result_list <- list()

nd <- 0

cel_dir <- dir()[grep("GSE", dir())]

for (i in cel_dir) {
  
  setwd(i)
  
  nd <- nd + 1
  
  dat <- affy::ReadAffy()
  
  eset <- plier::justPlier(dat)
  
  expre_mat <- affy::exprs(eset)
  
  result_list[[nd]] <- expre_mat
  
  setwd("..")
  
}

result_list

# names(result_list) <- cel_dir

eset_plier <- result_list

save(eset_plier, file = "plier_exprs.RData")

###--------------------------- S6. VSN algorithm ----------------------------###

setwd("GSE datasets about C.elegans (GPL200)")

result_list <- list()

nd <- 0

cel_dir <- dir()[grep("GSE", dir())]

for (i in cel_dir) {
  
  setwd(i)
  
  nd <- nd + 1
  
  dat <- affy::ReadAffy()
  
  eset <- vsn::vsnrma(dat)
  
  expre_mat <- affy::exprs(eset)
  
  result_list[[nd]] <- expre_mat
  
  setwd("..")
  
}

result_list

eset_vsn <- result_list

save(eset_vsn, file = "vsn_exprs.RData")

# End of this line. 


