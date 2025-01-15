# author : lwm
# date : 2022-9-7
date()
### Get the parameters
parser = argparse::ArgumentParser(description="Script to subcelltypes Clustering of cRNA data")
parser$add_argument('-I','--input', help='input a RDS')
parser$add_argument('-F','--subcelltypes', help = "input a file which sub celltypes do you want")
parser$add_argument('--prefix',help = "A name of output file")
parser$add_argument('-T','--tissue',help = "A name of tissue")
parser$add_argument('--celltypes',help = "A name of celltypes Idents")
parser$add_argument('-W','--method',help = "NormScale or SCT")
parser$add_argument('-D','--Dim',default = 30,type = 'integer',help = "A number of Dim to PCA and clsuter etc.")
parser$add_argument('-R','--resolution',default = 0.5,type = 'double',help = "resolution of cluster")
parser$add_argument('--tissueplot',default = NULL, help = "New name of tissue ! only a RDS including bm or blood !")
parser$add_argument('--topN', default = 1, type = 'integer', help = 'Top N marker genes in each cluster')
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)
set.seed(1234)

## Function ##
# Read filelist
Extract_SubCellType <- function(RDS,subcelltypes){
	Idents(RDS) <- opts$celltypes
	print(unique(Idents(RDS)))
	RDS <- subset(RDS,idents = subcelltypes)
	return(RDS)
}

# sub  to ScaleData,PCA,UNAP,FindNeighbors,FindCluster,TSNE
#Do_subClustering_Analyze <- function(RDS,method,dim,resolution){
#	if(method == "NormScale"){
#		RDS <- FindVariableFeatures(object = RDS,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
#  		RDS <- ScaleData(object = RDS, verbose = FALSE)
#	}
#  	else if(method == "SCT"){
#		RDS <- SCTransform(RDS, assay = "RNA", new.assay.name = "SCT", verbose = FALSE)	
#	}else {
#		stop("Please check your para!")
#	}
#	RDS <- RunPCA(object = RDS, npcs = dim, verbose = FALSE)
#	RDS <- RunUMAP(object = RDS, reduction = "pca",dims = 1:dim)
#	RDS <- FindNeighbors(object = RDS,k.param=20,dims = 1:dim)
#	RDS <- FindClusters(object = RDS,resolution = resolution)
#	RDS <- RunTSNE(object = RDS,dims = 1:dim,seed.use = 2)
#	DefaultAssay(RDS) <- "RNA"
#	return(RDS)
#}
# Integrate to ScaleData,PCA,UNAP,FindNeighbors,FindCluster,TSNE 
Do_Integrate_Analyze <- function(obj.list,dim,resolution){
  reference.list <- obj.list
  features <- SelectIntegrationFeatures(object.list = reference.list)
  anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = features)
  integrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(object = integrated) <- "integrated"
  integrated <- ScaleData(object = integrated, verbose = FALSE)
  integrated <- RunPCA(object = integrated, npcs = dim, verbose = FALSE)
  integrated <- RunUMAP(object = integrated, reduction = "pca",
                            dims = 1:dim)
  integrated = FindNeighbors(object = integrated,k.param=20,dims = 1:dim)
  integrated <- FindClusters(object = integrated,resolution = resolution)
  integrated <- RunTSNE(object = integrated,dims = 1:dim)
  DefaultAssay(integrated) <- "RNA"
  return(integrated)
}
# Marker Visualize #
Do_Marker_Visualize <- function(marker,RDS,topN){
	all.cluster.topn <- marker %>% group_by(cluster) %>% top_n(n = topN,wt = avg_log2FC)
	write.table(all.cluster.topn,file = paste(output,"/",opts$tissue,"_",opts$prefix,"_marker_Gene_top",opts$topN,".txt",sep  =""))
	target_gene <- unique(all.cluster.topn$gene)
	if(length(target_gene) > 40){
		target_gene <- target_gene[1:40]
		print("Warning! Only 40 target genes were selected!")
	}
	pdf(file = paste(output,'/',opts$tissue,"_", opts$prefix, '_cluster_markers_top',opts$topN,'_dotPlot.pdf', sep = ""))
	p <- DotPlot(RDS, features = target_gene) + RotatedAxis() +  theme(axis.text = element_text(size = 6))
	print (p)
	while (!is.null(dev.list()))  dev.off()
	pdf(file = paste(output,'/',opts$tissue,"_", opts$prefix, '_cluster_markers_top',opts$topN,'_VinPlot.pdf', sep = ""))
	p <- VlnPlot(RDS, features = target_gene, stack = TRUE) + theme(axis.text = element_text(size = 6))
	print (p)
	while (!is.null(dev.list()))  dev.off()
	pdf(file = paste(output,'/',opts$tissue,"_", opts$prefix, '_cluster_markers_top',opts$topN,'_FeaturePlot_UMAP.pdf', sep = ""),w = 16, h = 4*(length(target_gene)%/%3))
	p <- FeaturePlot(RDS, features = target_gene, reduction = "umap", pt.size = 0.5)
	print (p)
	while (!is.null(dev.list()))  dev.off()
	pdf(file = paste(output,'/',opts$tissue,"_", opts$prefix, '_cluster_markers_top',opts$topN,'_FeaturePlot_TSNE.pdf', sep = ""),w = 16, h = 4*(length(target_gene)%/%3))
	p <- FeaturePlot(RDS, features = target_gene, reduction = "tsne", pt.size = 0.5)
	print (p)
	while (!is.null(dev.list()))  dev.off()
}

Do_UmapPlot <- function(RDS){
	# PCA Plot#
	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_PCA.pdf",sep = ""),w = 10,h = 8)
	p1 <- DimPlot(RDS,reduction = "pca",group.by = "orig.ident")
	print(p1)
	while (!is.null(dev.list()))  dev.off()
	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_PCA_Group.pdf",sep = ""),w = 10, h = 8)
	p2 <- DimPlot(RDS,reduction = "pca",group.by = "mice")
	print(p2)
	while (!is.null(dev.list()))  dev.off()
	# Umap Plot #
	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap.pdf",sep = ""),w = 10,h = 8)
	p3 <- DimPlot(RDS,reduction = "umap",group.by = "orig.ident",label = F)
	print(p3)
	while (!is.null(dev.list()))  dev.off()

	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_Group.pdf",sep = ""),w = 10, h = 8)
	p4 <- DimPlot(RDS,reduction = "umap",group.by = "mice")
	print(p4)
	while (!is.null(dev.list()))  dev.off()

	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_Cluster.pdf",sep = ""),w = 10, h = 8)
	p5 <- DimPlot(RDS,reduction = "umap",group.by = "seurat_clusters",label = T, label.size = 5)
	print(p5)
	while (!is.null(dev.list()))  dev.off()

	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_split.pdf",sep = ""),w = 18, h = 8)
	p6 <- DimPlot(RDS,reduction = "umap",split.by = "mice",label = T, label.size = 5)
	print(p6)
	while (!is.null(dev.list()))  dev.off()

	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_platform.pdf",sep = ""),w = 10, h = 8)
	p7 <- DimPlot(RDS,reduction = "umap",group.by = "platform",label = F)
	print(p7)
	while (!is.null(dev.list()))  dev.off()
	
	#if(opts$tissueplot == "tissue"){
        #	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_tissue_group.pdf",sep = ""),w = 10, h = 8)
        #	p8 <- DimPlot(RDS,reduction = "umap",group.by = "tissue",label = F)
        #	print(p8)
       	#	while (!is.null(dev.list()))  dev.off()
        #	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_tissue_split.pdf",sep = ""),w = 21, h = 7)
        #	p9 <- DimPlot(RDS,reduction = "umap",split.by = "tissue",label = T, label.size = 5)
       	#	print(p9)
        #	while (!is.null(dev.list()))  dev.off()
	#} else if (opts$tissueplot == "tissue2"){
        #	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_tissue_group.pdf",sep = ""),w = 10, h = 8)
        #	p8 <- DimPlot(RDS,reduction = "umap",group.by = "tissue2",label = F)
       #	print(p8)
        #	while (!is.null(dev.list()))  dev.off()
        #	pdf(file = paste(output,"/",opts$tissue,"_",opts$prefix,"_",opts$method,"_Umap_tissue_split.pdf",sep = ""),w = 21, h = 7)
        #	p9 <- DimPlot(RDS,reduction = "umap",split.by = "tissue2",label = T, label.size = 5)
        #	print(p9)
        #	while (!is.null(dev.list()))  dev.off()
	#} else {
        #	print("Don't draw plot of tissue")
	#}
}

# main #
library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(presto)
library(ggplot2)
library(cowplot)

# Input RDS #
RDS <- readRDS(opts$input)
#DefaultAssay(RDS) <- "integrated" # 2022-10-20
print(DefaultAssay(RDS))
print("Step 1 : input  RDS Finish!")

# Input File which celltypes do you want 
#file <- read.csv(opts$subcelltypes) 
#subcelltype <- file$celltype
subcelltype <- c("1","10")
print(subcelltype)
print("Step 2 : input file Finish!")


# dir create 
output <- paste(opts$out,"/",opts$prefix,sep="")
dir.create(output)
print(output)
print("Step 3 : dir create Finish!")

# Extract subCelltype
subRDS <- Extract_SubCellType(RDS,subcelltypes = subcelltype)
print(unique(subRDS@meta.data$celltypes))

#saveRDS(subRDS,file=paste0(output,"/",opts$tissue,"_",opts$prefix,"_median.rds"))
print("Step 4 : Extract sub Celltypes")

# Extract plotform and mice
Extract_sub <- function(RDS,Ident,sub){
	Idents(RDS) <- Ident
	print(unique(Idents(RDS)))
	RDS <- subset(RDS,ident = sub)
}
#X10 <- Extract_sub(subRDS,Ident = "platform",sub = "10X")
C4 <- Extract_sub(subRDS,Ident = "platform",sub = "C4")
#GF_10x <- Extract_sub(X10,Ident = "mice" ,sub = "GF")
#SPF_10x <- Extract_sub(X10,Ident = "mice",sub = "SPF")
GF_C4 <- Extract_sub(C4,Ident = "mice" ,sub = "GF")
SPF_C4 <- Extract_sub(C4,Ident = "mice" ,sub = "SPF")
obj.list <- list(GF_C4 = GF_C4,SPF_C4 = SPF_C4)
saveRDS(obj.list,file=paste0(output,"/",opts$tissue,"_",opts$prefix,"_obj.list.rds"))

print("Extract")

# sub Clustering  #
#subRDS <- Do_subClustering_Analyze(subRDS,method = opts$method,dim = opts$Dim, resolution = opts$resolution)
#saveRDS(subRDS,file=paste0(output,"/",opts$tissue,"_",opts$prefix,"_Integrate.rds"))
#print(DefaultAssay(RDS))
print("Step 5 : sub Clustering Integrate Finish!")
# Integrate 
subRDS <- Do_Integrate_Analyze(obj.list = obj.list,dim = opts$Dim ,resolution = opts$resolution)
saveRDS(subRDS,file=paste0(output,"/",opts$tissue,"_",opts$prefix,"_Integrate.rds"))
print("Integrate")

# Do Umap
Do_UmapPlot(subRDS)
print("Step 6 : Do Umap Plot !")

# FindAllMarker
inputmarkers <- FindAllMarkers(subRDS)
write.table(inputmarkers,paste0(output,"/",opts$tissue,"_",opts$prefix,"_Marker.xls"),sep="\t",quote = FALSE)
print("Step 7 : Whole Find All maker Finish!")

# Marker visualize 
Do_Marker_Visualize(marker = inputmarkers,RDS = subRDS,topN = opts$topN)
print("Step 8 : Marker Visualize Done!")

# Count the number of cells in each cluster with GF and SPF
stat <- table(subRDS@meta.data$seurat_clusters,subRDS@meta.data$mice)
write.table(stat,file = paste0(output,"/",opts$tissue,"_",opts$prefix,"_stat.txt"),sep = "\t")
print("Step 9 : Stat Finish!")
#
date()
