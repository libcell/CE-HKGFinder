
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
### code chunk number 03: Obtain the statistics informations about GPL200.
### ****************************************************************************

### Step-01. Get the all GSE datasets with GPL200.

library(RCurl)
library(XML)

gpl_id <- "200"
gse_list <- NULL

for (i in 1:5) {
  
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/browse/?view=series&platform=", 
                gpl_id, 
                "&zsort=date&display=500&page=", 
                i)
  
  wp <- getURL(url)
  doc <- htmlParse(wp, asText = TRUE)
  tables <- readHTMLTable(doc)
  record <- tables[[1]]
  record <- record[-1, ]
  # DT::datatable(record[1:6, ])
  # print(head(record))
  
  gse_list <- rbind(gse_list, record)
  
}

gse_list <- as.data.frame(gse_list)

# DT::datatable(gse_list)

save(gse_list, file = "GSE_GPL200.RData")

openxlsx::write.xlsx(gse_list, file = "GEO-Series_C.elegans.xlsx")


### Step-02. Get the all GSE samples with GPL200.

gpl_id <- "200"
sam_list <- NULL

for (i in 1:5) {
  
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&platform=", 
                gpl_id, 
                "&zsort=date&display=500&page=", 
                i)
  
  wp <- getURL(url)
  doc <- htmlParse(wp, asText = TRUE)
  tables <- readHTMLTable(doc)
  record <- tables[[1]]
  record <- record[-1, ]
  # DT::datatable(record[1:6, ])
  # print(head(record))
  
  sam_list <- rbind(sam_list, record)
  
}

sam_list <- as.data.frame(sam_list)

# DT::datatable(sam_list)

# save(sam_list, file = "Sample_GPL200.RData")

openxlsx::write.xlsx(sam_list, file = "S1-GEO-Sample_C.elegans.xlsx")

# End of this line. 

