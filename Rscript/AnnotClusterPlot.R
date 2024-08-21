# author : Lwm
# data : 2022-9-13
# read file #
# updata :  add a file of celltypes for DEG
parser = argparse::ArgumentParser(description="Script to celltypes annotation  and cluster plot")
parser$add_argument('-I','--input', help='input rds,No annotation or No Newcluster')
parser$add_argument('-C',"--CTlist",help = "file of celltypes to annotation")
parser$add_argument('--prefix',help = "tissue or integrate prefix")
parser$add_argument('-A','--annotation',help = "Do Annotation or not !")
parser$add_argument('-M',"--method",help = "NormScale  or SCT, a prefix  of  output name")
parser$add_argument("-T",'--tissueplot', help='tissue or tissue2')
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)


# Function #
Do_annotation <- function(integrate,filelist){
  for(i in 1:length(filelist$clusters)){
      print(i)
      print(filelist$clusters[i])
      print(filelist$annotation[i])
      cells.use <- WhichCells(integrate, idents = filelist$clusters[i])
      integrate <- SetIdent(integrate, cells = cells.use, value = filelist$annotation[i])
  }
  return(integrate)
}
Do_NewCluster <- function(integrate,filelist){
  for(i in 1:length(filelist$annotation)){
      print(i)
      print(filelist$annotation[i])
      print(filelist$newclusters[i])
      cells.use <- WhichCells(integrate, idents = filelist$annotation[i])
      integrate <- SetIdent(integrate, cells = cells.use, value = filelist$newcluster[i])
  }
  return(integrate)
}


# main #
# library #
library(Seurat)
library(ggplot2)
library(dplyr)
# Read RDS #
integrate  <- readRDS(opts$input)
print( " Read RDS  Finish !")

# Read file with celltypes ,newclusters #
filelist <- read.csv(opts$CTlist,header = T,sep = "\t")
print( " Read filelist Finish !")

# Do annotation #
if(opts$annotation == TRUE){
	integrate <- Do_annotation(integrate,filelist)
	print(unique(integrate@active.ident))
	integrate@meta.data$celltypes <- Idents(integrate)
	print("Do annotation Finish!")
} else {
	print("Don't do annotation !")
}

# 
file <- filelist %>% select(annotation,newclusters)
NewClusters <- file[!duplicated(file),]
NewClusters <- NewClusters %>% arrange(newclusters)
print("Newclsuters file ")
print(NewClusters)

# Do New Cluster #
integrate <- Do_NewCluster(integrate = integrate,filelist = NewClusters)
print(unique(integrate@active.ident))
integrate@meta.data$clusters <- Idents(integrate)
saveRDS(integrate,file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_NewCluster.rds",sep = ""))
print("Do New Cluster Finish !")
#

#testis #
#Idents(integrate) <- "clusters"
#intestines2 Epithelial
#Idents(integrate) <- "seurat_clusters"
# new cluster #
# plot #
pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_NewCluster.pdf",sep = ""),w = 10, h = 8)
p0 <- DimPlot(integrate,reduction = "umap",label = T, label.size = 6) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
print(p0)
dev.off()
#  mice Umap #
pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_mice.pdf",sep = ""),w = 10, h = 8)
p1 <- DimPlot(integrate,reduction = "umap",group.by = "mice") + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
print(p1)
dev.off()
pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_mice_split.pdf",sep = ""),w = 20, h = 8)
p2 <- DimPlot(integrate,reduction = "umap",split.by = "mice",label = T, label.size = 6) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
print(p2)
dev.off()

# platform #
len <- length(unique(integrate@meta.data$platform))
print(len)
if (len > 1) {
	pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_platform.pdf",sep = ""),w = 10, h = 8)
	p1 <- DimPlot(integrate,reduction = "umap",group.by = "platform") + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
	print(p1)
	dev.off()
	pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_platfrom_split.pdf",sep = ""),w = 20, h = 8)
	p2 <- DimPlot(integrate,reduction = "umap",split.by = "platform",label = T, label.size = 6) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
	print(p2)
	dev.off()
} else {
	print("Don't draw a plot of platform ! ")
}

# tissue #
if(opts$tissueplot == "tissue"){
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_group.pdf",sep = ""),w = 10, h = 8)
        p8 <- DimPlot(integrate,reduction = "umap",group.by = "tissue",label = F) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
        print(p8)
        dev.off()
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_split.pdf",sep = ""),w = 21, h = 7)
	p9 <- DimPlot(integrate,reduction = "umap",split.by = "tissue",label = T, label.size = 6) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
	print(p9)
	dev.off()
} else if (opts$tissueplot == "tissue2"){
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_group.pdf",sep = ""),w = 10, h = 8)
        p8 <- DimPlot(integrate,reduction = "umap",group.by = "tissue2",label = F) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
        print(p8)
        dev.off()
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_split.pdf",sep = ""),w = 21, h = 7)
        p9 <- DimPlot(integrate,reduction = "umap",split.by = "tissue2",label = T, label.size = 6 , ncol = 4) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20)) # ncol = 4
        print(p9)
        dev.off()
} else {
	print("Don't draw plot of tissue")
}

# old cluster #
pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_Seurat_clusters.pdf",sep = ""),w = 10, h = 8)
p0 <- DimPlot(integrate,reduction = "umap",group.by = "seurat_clusters",label = T, label.size = 6) + theme(axis.title = element_text(size = 20)) +theme(axis.text = element_text(size = 20))
print(p0)
dev.off()

# annotation #
pdf(file = paste(opts$out,"/",opts$prefix,"_",opts$method,"_Umap_celltypes.pdf",sep = ""),w = 10, h = 8)
p0 <- DimPlot(integrate,reduction = "umap",group.by = "celltypes",label = T, label.size = 6) + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 20))
print(p0)
dev.off()

#stat # 
stat <- table(integrate@meta.data$celltypes,integrate@meta.data$mice)
write.table(stat,file = paste0(opts$out,"/",opts$prefix,'.stat.txt'),sep = "\t",quote = F)

print("Finish!")
