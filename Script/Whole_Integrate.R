# author : lwm
# date : 2022-9-7
date()
### Get the parameters
parser = argparse::ArgumentParser(description="Script to a integarte of  Cluster scRNA data")
parser$add_argument('-I','--input', help='input filelist of rds and sample ID')
parser$add_argument('--prefix',help = "A name of integrate")
parser$add_argument('-W','--method',help = "NormScale or SCT")
parser$add_argument('--topN', default = 1, type = 'integer', help = 'Top N marker genes in each cluster')
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)

## Function ##
# Read filelist
ReadRDS_list <- function(filelist){
  for(i in 1:dim(filelist)[1]){ # 2022-9-19
      print(dim(filelist)[1]) 
      obj.list[[i]] <- readRDS(filelist$path[i]) # 2022-9-19
  }
  return(obj.list)
}

# Integrate to ScaleData,PCA,UNAP,FindNeighbors,FindCluster,TSNE
Do_Integrate_Analyze <- function(obj.list,method){
	reference.list <- obj.list
	if(method == "NormScale"){
  		features <- SelectIntegrationFeatures(object.list = reference.list) #default nfeatures = 2000
  		anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = features)
		integrated <- IntegrateData(anchorset = anchors)
		DefaultAssay(object = integrated) <- "integrated"
		integrated <- ScaleData(object = integrated, verbose = FALSE)
	}
  	else if(method == "SCT"){
		features <- SelectIntegrationFeatures(object.list = reference.list,nfeatures = 3000)
		reference.list <- PrepSCTIntegration(object.list = reference.list,anchor.features = features)
		anchors <- FindIntegrationAnchors(object.list = reference.list,normalization.method = "SCT",features = features)
		integrated <- IntegrateData(anchorset = anchors)
	}else {
		stop("Please check your para!")
	}
	integrated <- RunPCA(object = integrated, npcs = 30, verbose = FALSE)
	integrated <- RunUMAP(object = integrated, reduction = "pca",dims = 1:30)
	integrated = FindNeighbors(object = integrated,k.param=20,dims = 1:30)
	integrated <- FindClusters(object = integrated,resolution = 0.5)
	integrated <- RunTSNE(object = integrated,dims = 1:30,seed.use = 2)
	DefaultAssay(integrated) <- "RNA"
	return(integrated)
}
# Marker Visualize #
Do_Marker_Visualize <- function(Integrate,marker,topN){
  all.cluster.topn <- marker %>% group_by(cluster) %>% top_n(n = topN,wt = avg_log2FC)
  write.table(all.cluster.topn,file = paste(opts$out,"/",opts$prefix,"_marker_Gene_top",opts$topN,".txt",sep  =""),sep  ="\t")
  target_gene <- unique(all.cluster.topn$gene)
  if(length(target_gene) > 40){
    target_gene <- target_gene[1:40]
    print("Warning! Only 40 target genes were selected!")
  }
  pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_dotPlot.pdf', sep = ""))
  p <- DotPlot(Integrate, features = target_gene) + RotatedAxis()
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_VinPlot.pdf', sep = ""))
  p <- VlnPlot(Integrate, features = target_gene, stack = TRUE)
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_FeaturePlot_UMAP.pdf', sep = ""),w = 16, h = 4*(length(target_gene)%/%3))
  p <- FeaturePlot(Integrate, features = target_gene, reduction = "umap", pt.size = 0.5)
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_FeaturePlot_TSNE.pdf', sep = ""),w = 16, h = 4*(length(target_gene)%/%3))
  p <- FeaturePlot(Integrate, features = target_gene, reduction = "tsne", pt.size = 0.5)
  print (p)

  while (!is.null(dev.list()))  dev.off()
}
# main #
library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)
#library(presto)
library(ggplot2)
library(cowplot)

# Read filtlist of RDS #
filelist <- read.csv(opts$input,header = T, sep = "\t")
print("Step1 : Read filelist of RDS Finish!")

# Integrate #
obj.list <- list()
obj.list <- ReadRDS_list(filelist)
Integrate <- Do_Integrate_Analyze(obj.list,method = opts$method)
saveRDS(Integrate,file=paste0(opts$out,"/",opts$prefix,"_Integrate.rds"))
print("Step2 : Whole Integrate Finish!")

# FindAllMarker
inputmarkers <- FindAllMarkers(Integrate)
write.table(inputmarkers,paste0(opts$out,"/",opts$prefix,"_Marker.xls",sep = ""),sep="\t",quote = FALSE)
print("Step3 : Whole Find All maker Finish!")

# Marker visualize 
Do_Marker_Visualize(marker = inputmarkers,Integrate = Integrate,topN = opts$topN)
print("Step4 : Marker Visualize Done!")

# Count the number of cells in each cluster with GF and SPF
stat <- table(Integrate@meta.data$seurat_clusters,Integrate@meta.data$mice)
write.table(stat,file = paste0(opts$out,"/",opts$prefix,"stat.txt",sep = ""),sep = "\t")
print("Step5 : Stat Finish!")
#
date()
