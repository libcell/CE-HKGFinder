
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
### code chunk number 02: Download the datasets, GPL200.
### ****************************************************************************

### Step-01. Obtain the all GSE datasets with GPL200.

library(GEOquery)
dir.create("GSE datasets about C.elegans (GPL200)")
setwd("GSE datasets about C.elegans (GPL200)")

### get the GSE Series Number. 

gpl <- getGEO('GPL200') # C.elegans. 

gpl_gse <- gpl@header$series_id

### download the GSE Series datasets. 

gse_download <- gpl_gse

gse_non <- NULL

for (i in gse_download) {
  
  test <- try(getGEOSuppFiles(i, makeDirectory = TRUE, baseDir = getwd(),
                              fetch_files = TRUE, filter_regex = NULL), 
              silent=TRUE)
  
  if (class(test) == "NULL") {
    
    gse_non <- c(gse_non, i)
    
    next
    
  }
  
}

gse_downloaded <- dir()

gse_non

# save(gse_non, file = "gse_non.RData")

write.csv(gse_non, file = "undownloaded_datasets.csv")

gpl_gse[-match(gse_downloaded, gpl_gse)]

dir.create("Other_processed_data")
setwd("Other_processed_data")

dir.create("Expression datasets")
setwd("Expression datasets")

library(ArrayExpress)

gse_non2 <- NULL

gse_non_AE <- paste("E-GEOD", (gsub("GSE", "", gse_non)), sep = "-") 

for (m in gse_non_AE) {
  
  dat <-  try(getAE(m, type = "processed"), silent = TRUE)
  
  if (class(dat) == "NULL") {
    
    gse_non2 <- c(gse_non2, i)
    
    next
    
  }
}

gse_non2

del_no <- setdiff(1:length(dir()), grep("processed-data", dir()))

file.remove(dir()[del_no])


setwd("..")

dir.create("array_design")

setwd("array_design")

getAE(gse_non_AE[1], type = "processed")

del_no <- setdiff(1:length(dir()), grep("A-AFFY-", dir()))

file.remove(dir()[del_no])

dir()

write.csv(gse_non, file = "gse_non.csv")

setwd("..")

setwd("..")

# source("http://bioconductor.org/biocLite.R")
# biocLite("ArrayExpress")

# End. 

