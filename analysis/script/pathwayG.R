# author : lwm
# date : 2024 - 11 - 15
# #
library(argparse)
parser = argparse::ArgumentParser(description="Script for Time Series DEG Analyze ")
parser$add_argument('-I','--input',type = "character", help='input rds  object of seurat')
parser$add_argument('-N','--name',type = "character", help='input name of output file ')
parser$add_argument('-P','--padj',type = "integer", default = 0.05, help='pajust of DEG')
parser$add_argument('-L','--log2FC', type = "integer", default = 0.5, help='log2fc of DEG')
parser$add_argument('-G','--organism',type = "character",help = "index of Mouse or Human")
parser$add_argument('-M','--method',type = "character",help = "method of GO or KEGG")
parser$add_argument('-C','--condition',type = "character",help = "condition index to split")
parser$add_argument('-K','--ko',type = "character",help = "kegg organism , mmu or hsa ")
parser$add_argument('-S','--simplify',type = "logical",default = FALSE,help = "use a simplify function ")
parser$add_argument('-O','--out',type = "character", help='out pathway')
opts = parser$parse_args()
# parameter #
print("Parameter of opts")
print(opts)

## main ##
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(patchwork))
suppressMessages(library(gridExtra))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ggpubr))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(MAST))
suppressMessages(library(reshape2))
suppressMessages(library(tidyr))

# input #
markers <- read.csv(opts$input,sep = "\t")
print("input Done!") 

# function # 
GO <- function(markers,padj,log2FC,organism,simplify = FALSE){
  EnrichGeneGo <- function(gene_set){
    Go <- enrichGO(gene_set, OrgDb = organism, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
    if(simplify){
      Go <- simplify(Go, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL )
    }
    if(!is.null(Go)){
      Go  <- setReadable(Go,OrgDb = organism,keyType = "ENTREZID")  # head # 
      Go  <- as.data.frame(Go)
    } else {
      message("Go is NULL !")
      Go <- data.frame()
    }
    return(Go)
  }
  ids = bitr(markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism)
  markers = merge(markers, ids, by.x='gene', by.y='SYMBOL')
  markers <- markers[which((markers$p_val_adj < padj) & ((markers$avg_log2FC > log2FC ) | (markers$avg_log2FC < -log2FC ))),]
  markers$group <- factor(ifelse(markers$avg_log2FC < 0, -1, 1), levels = c(-1, 1))
  g <- c(1,-1)
  d <- c("Up","Down")
  dat <- data.frame()
  for(i in 1:length(g)){
    genes <- subset(markers, group == g[i])$ENTREZID
    if(length(genes) != 0){
      message(paste("Num of ",d[i]," genes :",length(genes),sep = ""))
      data <- EnrichGeneGo(gene_set = genes)
      if(nrow(data) != 0 ){
        data <- data %>%  dplyr::mutate(group=rep(g[i], n())) %>% dplyr::mutate(deg = rep(d[i],n()))
      } 
    } else {
      data <- data.frame()
    }
    dat <- rbind(dat,data)
  }
  return(dat)
}
CondAllEnrich <- function(markers,celltypes,padj,log2FC,orgdb,method,kegg_organism = NULL,
                          prefix = NULL,output = NULL,save = F,simplify = F,condition = NULL,
                          condlst){
  GetPathway <- function(data,padj,log2FC,organism,method,simplify,kegg_organism){
    if(method == "GO"){
      pathway <- GO(markers = data, padj = padj, log2FC = log2FC, organism = organism,
                    simplify = simplify)
      
    } else if(method == "KEGG"  && !is.null(kegg_organism)){
      pathway <- KEGG(markers = data, padj = padj, log2FC = log2FC, organism = organism,
                      kegg_organism = kegg_organism)
    } else {
      stop("Please check your parameter of method or kegg_organism!")
    }
    return(pathway)
  }
  Ndata <- data.frame()
  if(orgdb == "Mouse"){
    organism <- "org.Mm.eg.db"
  } else if(orgdb == "Human"){
    organism <- "org.Hs.eg.db"
  } else {
    stop("Please check your parameter of orgdb !")
  }
  for(i in 1:length(unique(celltypes))){
    print(i)
    cl <- celltypes[i]
    print(cl)
    try({
      dat <- markers %>% subset(celltypes == cl)
      # condition list # 
      pathlst <- list()
      if(!is.null(condition)){
        for(c in condlst){
          print(c)
          data <- dat[dat[,condition] == c,]
          pathlst[[c]] <- GetPathway(data = data,padj = padj , log2FC = log2FC,organism = organism,
                                     method = method,
                                     simplify = simplify,kegg_organism = kegg_organism )
          if(nrow(pathlst[[c]]) !=0){
            pathlst[[c]][,condition] = c
          }
        }
        pathway <- do.call(rbind,pathlst)
        if(nrow(pathway) !=0){ 
          pathway$celltypes <- celltypes[i]
        }
        if(save){
          if(!is.null(perfix) && !is.null(output)){
            write.table(pathway,file = paste(output,"/",prefix,"_",i,"_",method,".txt",sep = ""),sep = "\t",quote = F)
          }
        }
        Ndata <- rbind(Ndata,pathway)
      }
    })
  }
  return(Ndata)
}
# running # 
cl <- unique(markers$celltypes)
print("celltypes")
print(cl)
time <- unique(markers$time) 
print("time")
print(time)

print("GO analyze running")
# tissue # 
tissue <- unique(markers$tissue)
lst <- list()
for(i in tissue){
  message(paste("tissue : ", i , sep = ""))
  mk <- markers[markers[,'tissue'] == i,]
  lst[[i]]  <- CondAllEnrich(markers = mk, celltypes = cl,
                             padj = opts$padj, log2FC = opts$log2FC,
                             orgdb = opts$organism, method = opts$method,
                             kegg_organism = opts$ko,simplify = opts$simplify,
                             condition = opts$condition,condlst = time)
  if(nrow(lst[[i]]) !=0){
    lst[[i]][,'tissue'] = i
  }
}
pathway <- do.call(rbind,lst)

f <- paste(opts$out,"/",opts$name,"_",opts$method,"_ALL_pathway.xls",sep = "")
print(f)
write.table(pathway,file = f,quote = F,row.names = F,sep = "\t")
print(" Finished ! ")

a <- data.frame()
b <- as.data.frame(a)
