#author : lwm
#Date: 2022-9-5
date()
# 
library(argparse, quietly = TRUE)
parser = argparse::ArgumentParser(description = 'Script for QC View  analysis of scRNA-seq ')
parser$add_argument('-i', '--input', dest = 'input', help = 'input a file of pathway of data')
parser$add_argument('-o', '--out', dest = 'out', help = 'a pathway of output')
parser$add_argument('-s', '--sample', dest = 'sample', help = 'name of tissue with GF or SPF')
parser$add_argument('-t', '--tissue', dest = 'tissue', help = 'name of tissue')
parser$add_argument('-m', '--mice', dest = "mice", help = 'GF or SPF index in meta.data')

opts = parser$parse_args()
print(opts)
# Function
# Read lsit of tissue pathway 
Read_List <- function(FileList) {
  object <- list() #  lwm 2023-10-12
  for (i in 1:length(FileList$tissue)){

    	object[[FileList$tissue[i]]] <- CreateSeuratObject(Read10X(FileList$path[i],gene.column = 1), # C4 gene.column = 1
           project = FileList$tissue[i],min.cells = 3,min.features = 200)
	object[[FileList$tissue[i]]] <- RenameCells(object[[FileList$tissue[i]]],add.cell.id = FileList$tissue[i])
    }
  return(object)
} 
# QC View
QC_View <- function(Object,species = "Human"){
    if (species == "Human"){
        ls <-c("^MT-","^RP[SL]","HB[^(P)]","PECAM1|PF4")
    }else if(species == "Mouse"){
        ls <- c("^mt-|^Mt-","^Rp[sl]","Hb[^(p)]","Pecam1|Pf4") 
    }else{
        stop("Check your species! Maybe add a pattern!")
    }
    Object[["percent.mt"]] <- PercentageFeatureSet(Object, pattern = ls[1])
    Object[["percent.ribo"]] <- PercentageFeatureSet(Object, pattern = ls[2])
    Object[["percent.hb"]] <- PercentageFeatureSet(Object, pattern = ls[3])
    Object[["percent.plat"]] <-  PercentageFeatureSet(Object,pattern = ls[4])
    return(Object)
}
# barplot #
hist_plot <- function(data) {
    pdf(file = paste(opts$out,"/",opts$sample,"_Histogram_Plot.pdf",sep = ""),w = 10,h = 8)
    p1 <- ggplot(data,aes(x = nCount_RNA)) + geom_histogram(binwidth = 200) + ggtitle("number of nCount_RNA ") + 
            theme(plot.title = element_text(hjust = 0.5,size = 25))
    p2 <- ggplot(data,aes(x = nFeature_RNA)) + geom_histogram(binwidth = 50) + ggtitle("number of nFeature_RNA") + 
            theme(plot.title = element_text(hjust = 0.5,size = 25))
    print(p1 + p2)
    dev.off()
}
# Anaylze #
# library #
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(DoubletFinder))
suppressMessages(library(patchwork))
# 
# read file of data pathway
filelist <- read.csv(opts$input,header=T,sep = "\t")
print(filelist)
print("Step1 : Read Finish!")

# Create some object of seurat
# object <- list() # 2023-10-11 # lwm 
obj.list <- Read_List(FileList = filelist)
print("Step2 : Craete object Finish!")

# Merge object of seurat 
print(length(obj.list))
Merge <- merge(obj.list[[1]],obj.list[2:length(obj.list)]) # C4

print(Merge)
print("Step3 : Merge object Finish!")

# QC View 
Merge <- QC_View(Merge,species = "Mouse")
pdf(file = paste(opts$out,"/",opts$sample,"_QC_ViolinPLot.pdf",sep = ''),w = 10,h = 8)
p3 <- VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.plat"), 
                  ncol = 3, pt.size = 0)
p3
dev.off()
pdf(file = paste(opts$out,"/",opts$sample,"_QC_Scatter.pdf",sep = ""),w = 10,h = 8)
p4 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
p5 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p4 + p5
dev.off()
print("Step4 : QC View Finish!")

# hist plot #
hist_plot(Merge@meta.data)
print("Step5 : Hist Plot Finish!")

# tissue, RNA, mice index in meda.data
Merge@meta.data$mice <- opts$mice
Merge@meta.data$tissue <- opts$tissue
print("Step6 : index Finish!")
# save 
saveRDS(Merge,paste(opts$out,"/",opts$sample,"_QC_View_seuratObject.rds",sep = ""))
print("Step7 : Save Finish!")
date()
