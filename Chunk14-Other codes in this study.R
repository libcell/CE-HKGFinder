
################################################################################
#    &&&....&&&    % Project: Identification of HKG candidates in C. elegans   #
#  &&&&&&..&&&&&&  % Author: Bo Li, Jingxin Tao, Youjin Hao                    #
#  &&&&&&&&&&&&&&  % Date: Dec. 24th, 2019                                     #
#   &&&&&&&&&&&&   %                                                           #
#     &&&&&&&&     % Environment: R version 3.5.3;                             #
#       &&&&       % Platform: x86_64-pc-linux-gnu (64-bit)                    #
#        &         %                                                           #
################################################################################

library(GEOquery)

gse <- get(load("gse.RData"))

gse <- as.character(gse)

gse.asthma.m <- list()

error.gse <- NULL

for (i in 1:length(gse)) {
  
  tmp <- try(getGEO(gse[i]), silent = TRUE)
  
  if (class(tmp) == "try-error") {error.gse <- c(error.gse, gse[i])}
  
  gse.asthma.m[[i]] <- tmp
  
  Sys.sleep(5)
  
}

names(gse.asthma.m) <- gse

save(gse.asthma.m, file = "gse_asthma_m.RData")



tpm <- read.delim("clipboard", header = TRUE)

x <-c("cyc-1",
      "tba-1",
      "atp-3",
      "mdh-1",
      "gpd-2",
      "eif-3.C",
      "act-1",
      "daf-16",
      "cdc-42",
      "pmp-3",
      "act-2",
      "csq-1",
      "ama-1",
      "rbd-1", 
      "rps-4", 
      "rpl-15"
)

sink("bb.txt")

f <- function(x) {
  
  mean(x, na.rm = TRUE)
  
}


for (i in x) {
  
  tmp <- which(as.character(tpm$Gene.name) == i)
  
  tmp <-  tpm[tmp, 2:184]
  
  print("##########################################################################")
  
  print(i)
  
  # print(tmp)
  
  s <- sd(as.numeric(tmp), na.rm = TRUE)
  
  m <- mean(as.numeric(tmp), na.rm = TRUE)
  
  # print(s); print(m)
  
  print(s/m)
  
}

sink()



#1.读取文件
#1.1设置路径

path <- "C:/Users/Li Bo/Desktop/srr"

setwd(path)

getwd()

fileNames <- dir(path)

#1.2解压文件
library(R.utils)

for (i in fileNames) {
  
  gunzip(i, remove = T) 
  
}

#1.3拼接成表格


# data1 <- cbind(rownames(data1), data1)

fileNames <- dir()

# data1 <- read.table("SRR962090.se.tsv")

data1 <- NULL

for (i in fileNames) {
  
  tmp <- read.table(i)[, 1]
  
  data1 <- cbind(data1, tmp)
  
}

rownames(data1) <- rownames(read.table(i))

colnames(data1) <- fileNames


#2.注释文件

library(biomaRt)

mart <- useMart("parasite_mart", dataset = "wbps_gene", 
                host = "https://parasite.wormbase.org", port = 443)

genes <- getBM(mart = mart, 
               filters = c("species_id_1010", "only_elegaprjna13758_homologue"),
               value = list("mansoprjea36577", TRUE),
               attributes = c("wbps_gene_id", "elegaprjna13758_gene", "elegaprjna13758_gene_name"))
head(genes)



library(org.Ce.eg.db)

y <-toTable(org.Ce.egWORMBASE) 


dat1 <-merge(x = data1, y = y, by.x= "rownames(data1)",by.y="wormbase_id",all.x= TRUE ) 

require(clusterProfiler)
eg = bitr(dat1$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Ce.eg.db")


dat2 <- merge(x = dat1,y = eg,by.x="gene_id",by.y="ENTREZID",all.x=TRUE)



#3.保存
setwd("..")
write.csv(dat2, "GSE50548counts.csv")



###############################################################################################################################
#4.验证
#4.2读取表格
setwd("/Users/touyasushishin/early C. elegans Embryo/")

dat <-read.csv("GSE50548counts.csv")

colnames(dat)

genename <- as.character(dat[,188])

expr1 <- dat[,4:187]

tpm = t(t(expr1)/colSums(expr1)) * 10^6

expr1 <- tpm

dat <- cbind(expr1,genename)

write.csv(dat,"tpm.csv")

#4.3对应基因的位置提取出
x1<-c("cyc-1","tba-1","atp-3","mdh-1","gpd-2","eif-3.C","act-1","daf-16","cdc-42","pmp-3","act-2","csq-1","ama-1","rbd-1")
x2<-c("rps-23","rps-27","rps-26","rps-4","rps-2","rps-16","rps-17","rpl-24.1","rpl-15","rpl-35","rpl-36","rpl-33","rpl-27")


dat <- read.csv("tpm.csv")
dat <- dat[,-1]

a1  <- NULL

aaa <- function(x){mean(x, na.rm = TRUE)}

for (i in x1) {
  a <- dat[which(dat[,"genename"]== i),-185]
  
  tmp1 <- apply(a, 2, aaa)
  
  a1 <-rbind(a1,tmp1)
}

rownames(a1) <- x1



a2 <- NULL

for (i in x2) {
  a<-dat[which(dat[,"genename"]== i),-185]
  
  tmp1 <- apply(a, 2, aaa)
  
  a2<-rbind(a2,tmp1)
}


rownames(a2) <- x2



#4.4相应提出来做成表格
a<-rbind(a1,a2)



############################cv#############################################
cv <-function(x){
  cv=sd(x)/mean(x)
  return(cv)
}

tmp <- apply(a, 1, cv) 

sort(tmp)


############################Gini##########################################
#BiocManager::install("ineq")

library(ineq)

tmp1 <- apply(a, 1, ineq::Gini)

sort(tmp1)
uresult <- cbind( tmp1, tmp)



##########5.标准化####################################################################
dat <-read.csv("GSE50548counts.csv")

colnames(dat)

genename <- as.character(dat[,188])

expr1 <- dat[,4:187]

library(edgeR)

y <- DGEList(counts=expr1)

y <- cpm(y)

dat <- cbind(y,genename)

write.csv(dat,"cpm.csv")



############################cv#############################################
x1<-c("cyc-1","tba-1","atp-3","mdh-1","gpd-2","eif-3.C","act-1","daf-16","cdc-42","pmp-3","act-2","csq-1","ama-1","rbd-1")
x2<-c("rps-23","rps-27","rps-26","rps-4","rps-2","rps-16","rps-17","rpl-24.1","rpl-15","rpl-35","rpl-36","rpl-33","rpl-27")


dat <- read.csv("cpm.csv")

colnames(dat)

dat <- dat[,-1]

a1  <- NULL

aaa <- function(x){mean(x, na.rm = TRUE)}

for (i in x1) {
  a <- dat[which(dat[,"genename"]== i),-185]
  
  tmp1 <- apply(a, 2, aaa)
  
  a1 <-rbind(a1,tmp1)
}

rownames(a1) <- x1



a2 <- NULL

for (i in x2) {
  a<-dat[which(dat[,"genename"]== i),-185]
  
  tmp1 <- apply(a, 2, aaa)
  
  a2<-rbind(a2,tmp1)
}


rownames(a2) <- x2



#4.4相应提出来做成表格
a<-rbind(a1,a2)





cv <-function(x){
  cv=sd(x)/mean(x)
  return(cv)
}

cmp <- apply(a, 1, cv) 

sort(cmp)


############################Gini##########################################
#BiocManager::install("ineq")

library(ineq)

cmp1 <- apply(a, 1, ineq::Gini)

sort(cmp1)

uresult <- cbind( cmp1, cmp)

dat <- read.delim("clipboard", header = TRUE)

dat <- dat[1:22548, ]

rownames(dat) <- dat$y

length(unique(dat$y))

g <- as.character(unique(dat$y))


dat <- na.omit(dat)

gene.ex <- function(x) {
  
  tmp <- dat[as.character(dat$y) == x, ]
  
  tmp <- apply(tmp[, 2:13], 2, mean) 
  
  tmp2 <- paste("CV is", sd(tmp)/mean(tmp), sep = " ")
  
  tmp <- summary(tmp)
  
  return(tmp2)
  
}

gene.ex("rps-4")

x <- "col-119"


