#author : lwm
#Date: 2022-9-5
date()
# 
library(argparse, quietly = TRUE)
parser = argparse::ArgumentParser(description = 'Script for QC View  analysis of scRNA-seq ')
parser$add_argument('-g', '--GF', dest = 'GF_input', help = 'input rds of GF')
parser$add_argument('-s', '--SPF', dest = 'SPF_input', help = 'input rds of SPF')
parser$add_argument('-o', '--out', dest = 'out', help = 'a pathway of output')
parser$add_argument('-f', '--prefix', dest = 'prefix', help = 'name of tissue')


opts = parser$parse_args()
print(opts)
# Function

# QC_View
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
    pdf(file = paste(opts$out,"/",opts$prefix,"_Histogram_Plot.pdf",sep = ""),w = 10,h = 8)
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
# read RDS
GF <-  readRDS(opts$GF_input)
SPF <- readRDS(opts$SPF_input)
print("Step1 : Read Finish!")

# Create some object of seurat


# Merge object of seurat 
Merge <- merge(GF,SPF)
print(Merge)
print("Step2 : Merge object Finish!")

# QC View 
Merge <- QC_View(Merge,species = "Mouse")
pdf(file = paste(opts$out,"/",opts$prefix,"_QC_ViolinPLot.pdf",sep = ''),w = 10,h = 8)
p3 <- VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.plat"), 
                  ncol = 3, pt.size = 0)
p3
dev.off()
pdf(file = paste(opts$out,"/",opts$prefix,"_QC_Scatter.pdf",sep = ""),w = 10,h = 8)
p4 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
p5 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p4 + p5
dev.off()


pdf(file = paste(opts$out,"/",opts$prefix,"_QC_Group_ViolinPLot.pdf",sep = ''),w = 10,h = 8)
p6 <- VlnPlot(Merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.plat"),
                  group.by = "mice",ncol = 3, pt.size = 0)
p6
dev.off()

pdf(file = paste(opts$out,"/",opts$prefix,"_QC_Group_Scatter.pdf",sep = ""),w = 10,h = 8)
p7 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "mice")
p8 <- FeatureScatter(Merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "mice")
p7 + p8
dev.off()
print("Step3 : QC View Finish!")

# hist plot #
hist_plot(Merge@meta.data)
print("Step4 : Hist Plot Finish!")

# save 
saveRDS(Merge,paste(opts$out,"/",opts$prefix,"_QC_View_seuratObject.rds",sep = ""))
print("Step5 : Save Finish!")
date()
