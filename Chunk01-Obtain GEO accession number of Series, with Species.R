
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
### code chunk number 01: Obtain GEO accession number of Series, with Species.
### ****************************************************************************

### Step-01. Obtain the all GSE datasets with Classic Affymetrix DNA Microarray.

# For windows, work directory. 
# setwd("F:/House-keeping genes/00-Preprocessed data/new/")

# For Mac OS, work directory. 

setwd("/Volumes/TOSHIBA EXT/03-重要数据备份/16-House-keeping genes/")

# Loading a R script, to install packages in bioconductor. 

source("http://bioconductor.org/biocLite.R")

# biocLite("limma")

library(GEOquery)

### get the GSE Series number for four Species. 

sc <- c("Homo_sapiens", "Mus_musculus", "Rattus_norvegicus", "Macaca_mulatta")

pf <- c("GPL570", "GPL1261", "GPL1355", "GPL3535")

GSEInfo <- list()

p <- 0

for (s in pf) {
  
  gpl <- getGEO(s)
  
  gpl_gse <- gpl@header$series_id
  
  gpl_sam <- gpl@header$sample_id
  
  sc_data <- list(GSE_number = gpl_gse, Sample_number = gpl_sam)
  
  # assign(sc[1], sc_data)
  
  p <- p + 1
  
  GSEInfo[[p]] <- sc_data
  
}

names(GSEInfo) <- sc

num_series <- c(length(GSEInfo$Homo_sapiens$GSE_number), 
                length(GSEInfo$Mus_musculus$GSE_number), 
                length(GSEInfo$Rattus_norvegicus$GSE_number), 
                length(GSEInfo$Macaca_mulatta$GSE_number))

num_sample <- c(length(GSEInfo$Homo_sapiens$Sample_number), 
                length(GSEInfo$Mus_musculus$Sample_number), 
                length(GSEInfo$Rattus_norvegicus$Sample_number), 
                length(GSEInfo$Macaca_mulatta$Sample_number))


id <- function(x) {
  y <- getGEO(x)
  y <- nrow(y@dataTable@table)
  return(y)
}

num_probe <- as.numeric(sapply(pf, id))

GSEtable <- data.frame(GEO_accession = pf, 
                       Number_of_Probes = num_probe, 
                       Number_of_Series = num_series, 
                       Number_of_Sample = num_sample)

rownames(GSEtable) <- sc


### Step-02. Save the all GSE datasets with Classic Affymetrix DNA Microarray.

print(GSEtable)

save(GSEInfo, file = paste(Sys.Date(), "Datasets.Rdata", sep = "_"))

### End of code chunk number 01. 
