## author : lwm
# date : 2022-9-10
date()
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input filelist of rds and sample ID')
parser$add_argument('-T','--tissue',help = "tissue of mouse")
parser$add_argument('-W','--method',help = "NormScale or SCT")
parser$add_argument('--tissueplot',help = "tissue umap plot")
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)

# main #
library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(presto)

# Read RDS #
Integrate <- readRDS(opts$input)
print("Step 0 : Read RDS Finish!")
# PCA Plot#
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_PCA.pdf",sep = ""),w = 10,h = 8)
p1 <- DimPlot(Integrate,reduction = "pca",group.by = "orig.ident")
print(p1)
dev.off()
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_PCA_Group.pdf",sep = ""),w = 10, h = 8)
p2 <- DimPlot(Integrate,reduction = "pca",group.by = "mice")
print(p2)
dev.off()
print("Step1 : PCA Plot Finish!")

# Umap Plot #
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap.pdf",sep = ""),w = 10,h = 8)
p3 <- DimPlot(Integrate,reduction = "umap",group.by = "orig.ident",label = F)
print(p3)
dev.off()

pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_Group.pdf",sep = ""),w = 10, h = 8)
p4 <- DimPlot(Integrate,reduction = "umap",group.by = "mice")
print(p4)
dev.off()

pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_Cluster.pdf",sep = ""),w = 10, h = 8)
p5 <- DimPlot(Integrate,reduction = "umap",group.by = "seurat_clusters",label = T, label.size = 5)
print(p5)
dev.off()

pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_split.pdf",sep = ""),w = 18, h = 8)
p6 <- DimPlot(Integrate,reduction = "umap",split.by = "mice",label = T, label.size = 5)
print(p6)
dev.off()

pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_platform.pdf",sep = ""),w = 10, h = 8)
p7 <- DimPlot(Integrate,reduction = "umap",group.by = "platform",label = F)
print(p7)
dev.off()

if(opts$tissueplot == "tissue"){
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_group.pdf",sep = ""),w = 10, h = 8)
        p8 <- DimPlot(Integrate,reduction = "umap",group.by = "tissue",label = F)
        print(p8)
        dev.off()
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_split.pdf",sep = ""),w = 21, h = 7)
	p9 <- DimPlot(Integrate,reduction = "umap",split.by = "tissue",label = T, label.size = 5)
	print(p9)
	dev.off()
} else if (opts$tissueplot == "tissue2"){
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_group.pdf",sep = ""),w = 10, h = 8)
        p8 <- DimPlot(Integrate,reduction = "umap",group.by = "tissue2",label = F)
        print(p8)
        dev.off()
	pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_tissue_split.pdf",sep = ""),w = 21, h = 7)
        p9 <- DimPlot(Integrate,reduction = "umap",split.by = "tissue2",label = T, label.size = 5)
        print(p9)
        dev.off()
} else {
	print("Don't draw plot of tissue")
}
print("Step2 : Umap Plot Finish!")
date()
