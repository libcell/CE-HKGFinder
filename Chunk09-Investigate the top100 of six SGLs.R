
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
### code chunk number 09: Investigation of top 100 genes of SGLs.
### ****************************************************************************

top50HKG <- get(load("J:/rank50_top5000_v2.RData"))

top100HKG <- get(load("J:/top5000_rank100.RData"))

gpl200anno <- get(load("J:/gpl200anno.RData"))

for (n in 1:length(top100HKG)) {
  
  probe.len <- length(top100HKG[[n]]) 
  
  at.pos <- grep("_at", top100HKG[[n]])
  
  at2sym <- NULL
  
  for (g in at.pos) {
    
    tp <- gpl200anno$`Gene Symbol`[match(top100HKG[[n]][g], gpl200anno$ID)]
    
    print(g)
    
    print(tp)
    
    tmp <- strsplit(tp, " /// ")[[1]][1]
    
    if (length(grep("WBGene", tmp)) != 0) {
      
      x <- strsplit(tp, " /// ")[[1]][2] 
      
      tmp <- x
      
    }
    
    print(tmp)
    
    at2sym <- c(at2sym, tmp)
    
  }
  
  top100HKG[[n]][at.pos] <- at2sym
  
}

names(top50HKG)[4] <- names(top100HKG)[4] <- "Li-Wong_exprs" 

save(top50HKG, file = "rank50_top5000_v2.RData")
save(top100HKG, file = "rank100_top5000_v2.RData")


library(UpSetR)

# help(package = "UpSetR")

## ---------------- Venn plot for top 50 genes in six SGLs ------------------ ##

listInput <- list(RMA = top50HKG$RMA_exprs, 
                  MAS5 = top50HKG$MAS_exprs, 
                  gcRMA = top50HKG$GCRMA_exprs, 
                  Li_Wong = top50HKG$`Li-Wong_exprs`, 
                  PLIER = top50HKG$PLIER_exprs, 
                  VSN = top50HKG$VSN_exprs)

df.listInput <- fromList(listInput)

col_seq <- c("#376092","#77933C","#993735","#604A7B","#31859C","#E46C0A" )

# #F08080 standing for the color named by LightCoral

hkgs50 <- upset(df.listInput, 
          sets = c ("RMA", "MAS5","gcRMA","Li_Wong", "PLIER","VSN"), 
          main.bar.color = "grey20", 
          queries = list(list(query = intersects, 
                              params = list("PLIER", "Li_Wong", "MAS5", "RMA"), 
                              color = "#F08080", active = TRUE), 
                         list(query = intersects, 
                              params = list("PLIER", "Li_Wong", "RMA"), 
                              color = "#F08080", active = TRUE)
                         ), 
          sets.bar.color = col_seq,
          order.by = "degree")


hk4.top50 <- Reduce(intersect, 
                    list(top50HKG$PLIER_exprs, 
                         top50HKG$`Li-Wong_exprs`, 
                         top50HKG$MAS_exprs, 
                         top50HKG$RMA_exprs)
                   )

## ---------------- Venn plot for top 100 genes in six SGLs ----------------- ##

listInput <- list(RMA = top100HKG$RMA_exprs, 
                  MAS5 = top100HKG$MAS_exprs, 
                  gcRMA = top100HKG$GCRMA_exprs, 
                  Li_Wong = top100HKG$`Li-Wong_exprs`, 
                  PLIER = top100HKG$PLIER_exprs, 
                  VSN = top100HKG$VSN_exprs)

df.listInput <- fromList(listInput)

col_seq <- c( "#376092","#77933C","#993735","#604A7B","#31859C","#E46C0A" )

upset(df.listInput, 
      sets = c ("RMA", "MAS5","gcRMA","Li_Wong", "PLIER","VSN"), 
      main.bar.color = "grey20", 
      queries = list(list(query = intersects, 
                          params = list("PLIER", "Li_Wong", "MAS5", "RMA"), 
                          color = "#F08080", active = TRUE), 
                     list(query = intersects, 
                          params = list("VSN", "Li_Wong", "RMA"), 
                          color = "#F08080", active = TRUE), 
                     list(query = intersects, 
                          params = list("PLIER", "Li_Wong", "RMA"), 
                          color = "#F08080", active = TRUE), 
                     list(query = intersects, 
                          params = list("PLIER", "MAS5", "RMA"), 
                          color = "#F08080", active = TRUE), # coral red
                     list(query = intersects, 
                          params = list("PLIER", "RMA"), 
                          color = "violet", active = TRUE) # violet
                     ), 
      sets.bar.color = col_seq,
      order.by = "degree")

hk4.top100 <- Reduce(intersect, 
                     list(top100HKG$PLIER_exprs, 
                          top100HKG$`Li-Wong_exprs`, 
                          top100HKG$MAS_exprs, 
                          top100HKG$RMA_exprs)
)

hk3.1.top100 <- Reduce(intersect, 
                       list(top100HKG$VSN_exprs, 
                            top100HKG$`Li-Wong_exprs`, 
                            top100HKG$RMA_exprs)
)

hk3.2.top100 <- Reduce(intersect, 
                       list(top100HKG$PLIER_exprs, 
                            top100HKG$`Li-Wong_exprs`, 
                            top100HKG$RMA_exprs)
)

hk3.3.top100 <- Reduce(intersect, 
                       list(top100HKG$PLIER_exprs, 
                            top100HKG$MAS_exprs, 
                            top100HKG$RMA_exprs)
)

hk2.top100 <- Reduce(intersect, 
                     list(top100HKG$PLIER_exprs, 
                          top100HKG$RMA_exprs)
)

hk2.top100 # rps-17

hk.3_4.100 <- unique(c(hk4.top100, hk3.1.top100, hk3.2.top100, hk3.3.top100))

# End of this line. 
