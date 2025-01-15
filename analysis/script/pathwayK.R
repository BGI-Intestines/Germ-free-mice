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
parser$add_argument('-K','--ko',type = "character",help = "kegg organism , mmu or hsa ")
parser$add_argument('-S','--simplify',type = "logical",default = FALSE,help = "use a simplify function ")
parser$add_argument('-C','--condition',type = "character",default = NULL,help = "which index of condition  to split") # bug ？
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
KEGG <- function(markers,padj,log2FC,organism,kegg_organism){
  EnrichGeneKegg <- function(gene_set){
    kegg <- enrichKEGG(gene_set, organism = kegg_organism, keyType = 'kegg', pvalueCutoff = 0.05,pAdjustMethod = 'BH',
                       minGSSize = 10,maxGSSize = 500, qvalueCutoff = 0.2,use_internal_data = FALSE)
    if(!is.null(kegg)){
      kegg <- setReadable(kegg,OrgDb = organism,keyType = "ENTREZID") 
      kegg <- as.data.frame(kegg)
    } else {
      message("kegg is NULL !")
      kegg <- data.frame()
    }
    return(kegg)
  }
  ids = bitr(markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=organism)
  markers = merge(markers, ids, by.x='gene', by.y='SYMBOL')
  markers <- markers[which((markers$p_val_adj < padj) & ((markers$avg_log2FC > log2FC ) | (markers$avg_log2FC < -log2FC ))),]
  markers$group <- factor(ifelse(markers$avg_log2FC < 0, -1, 1), levels = c(-1, 1))
  g <- c(1,-1)
  d <- c("Up","Down")
  dat <- data.frame()
  #return(markers)
  #stop('haha')
  for(i in 1:length(g)){
    genes <- subset(markers, group == g[i])$ENTREZID
    if(length(genes) != 0){
      message(paste("Num of ",d[i]," genes :",length(genes),sep = ""))
      data <- EnrichGeneKegg(gene_set = genes) 
      if(nrow(data) != 0 ){
        data <- data  %>% dplyr::mutate(group=rep(g[i], n())) %>% dplyr::mutate(deg = rep(d[i],n()))
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

print("KEGG analyze running")
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

# check KEGG # 
#data <- test %>% filter(time == "0h vs 24h") %>% filter(celltypes == "Villus SMC") %>% filter(tissue == "Duodenum")
#table(data$group)
#haha <- KEGG(markers = data,padj = 0.05,log2FC = 0.5,organism = "org.Mm.eg.db",kegg_organism = "mmu")
#haha
#head(haha)
#table(haha$group)
