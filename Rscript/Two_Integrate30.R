# author : lwm
# date : 2022-9-7
date()
### Get the parameters
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input filelist of rds and sample ID')
parser$add_argument('-C','--filterCounts', type = "double", help='counts number threhold')
parser$add_argument('-G','--filterGene', type = "double", help='gene number threhold')
parser$add_argument('-M','--filterMito', type = "double", help='mito ratio threhold')
parser$add_argument('-T','--tissue',help = "tissue of mouse")
parser$add_argument('-W','--method',help = "NormScale or SCT")
parser$add_argument('--topN', default = 1, type = 'integer', help = 'Top N marker genes in each cluster')
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)

## Function ##
# Read filelist
#ReadRDS_list <- function(filelist){
  #for(i in 1:length(filelist)){ # old 
#   for(i in 1:dim(filelist)[1])  #new 2022-9-18
#          obj.list[[filelist$mice[i]]] <- readRDS(filelist$path[i])
#  }
#  return(obj.list)
#}
ReadRDS_list <- function(filelist){ # nrow > 2
  print(dim(filelist)[1]) 
  for(i in 1:dim(filelist)[1]){  #new 2022-9-18
          obj.list[[i]] <- readRDS(filelist$path[i])
  }
  return(obj.list)
}
# Filter
QC_Filter <- function(obj.list,filterCounts,filterGene,filterMito){
    for(i in 1:length(obj.list)){
        obj.list[[i]] <- subset(obj.list[[i]],subset = nCount_RNA < filterCounts & nFeature_RNA > filterGene  & percent.mt < filterMito)
    }
  return(obj.list)
}
# Find Doublet
Find_doublet <- function(data,sct = FALSE){
  print(sct)
  sweep.res.list <- paramSweep_v3(data, PCs = 1:30, sct = sct) # use SCT : sct = TRUE
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  nExp_poi <- round(0.05*ncol(data))
  p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
  data <- doubletFinder_v3(data, PCs = 1:30, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data
}

# Norm,Scale,PCA,UMAP  to Find_Doublet
Do_Analyze <- function(obj.list,method){
  doublet=list()
  for (i in 1:length(x = obj.list)) {
    if(method == "NormScale"){
        obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
        obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                          selection.method = "vst", nfeatures = 2000, verbose = FALSE)
        obj.list[[i]] <- ScaleData(obj.list[[i]],features = rownames(obj.list[[i]]))
    } else if(method  == "SCT"){
        obj.list[[i]] <- SCTransform(obj.list[[i]],assay = "RNA",new.assay.name = "SCT",verbose = FALSE)
    } else {
        stop("Please check your parameter!")
    }
    obj.list[[i]] <- RunPCA(obj.list[[i]],features = VariableFeatures(object = obj.list[[i]]))
    obj.list[[i]] <- RunUMAP(obj.list[[i]], dims = 1:30)
    obj.list[[i]] <- Find_doublet(obj.list[[i]])
    doublet[[i]] <- obj.list[[i]]@meta.data[,-6]
    obj.list[[i]] <- subset(obj.list[[i]],subset=doublet_info=="Singlet")
  }
  
  #doublet_df=do.call(rbind,doublet)
  #write.table(doublet_df,file = paste0(opts$out,"/","Doublet_info_",opts$tissue,".txt"),sep="\t",quote=FALSE)
  saveRDS(doublet,file = paste0(opts$out,"/","Doublet_info_",opts$tissue,".rds"))
  return(obj.list)
}

# Integrate #
Do_Integrate_Analyze <- function(obj.list){
  reference.list <- obj.list
  features <- SelectIntegrationFeatures(object.list = reference.list)
  anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features = features)
  integrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(object = integrated) <- "integrated"
  integrated <- ScaleData(object = integrated, verbose = FALSE)
  integrated <- RunPCA(object = integrated, npcs = 30, verbose = FALSE)
  integrated <- RunUMAP(object = integrated, reduction = "pca",
                            dims = 1:30)
  integrated = FindNeighbors(object = integrated,k.param=20,dims = 1:30)
  integrated <- FindClusters(object = integrated,resolution = 0.5)
  integrated <- RunTSNE(object = integrated,dims = 1:30)
  DefaultAssay(integrated) <- "RNA"
  return(integrated)
}

# Find All marker #
Do_FindAllMarkers <- function(Integrate){
    Integrate.markers <- FindAllMarkers(Integrate, assay = "RNA", only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
    #Integrate.markers.auc <- wilcoxauc(Integrate.markers,seurat_assay = "RNA",'seurat_clusters')
    write.table(Integrate.markers,file = paste(opts$out,"/",opts$tissue,"_markerGene.txt",sep = ""),sep = "\t",quote = F,row.names = F)
    #write.table(Integrate.markers.auc,file = paste(opts$out,"/",opts$tissue,"_markerGene.txt",sep = ""),sep = "\t",quote = F,row.names = F)
    return(Integrate.markers)
}

# All marker visualize #
Do_Marker_Visualize <- function(Integrate,marker,topN){
  all.cluster.topn <- marker %>% group_by(cluster) %>% top_n(n = topN,wt = avg_log2FC)
  write.table(all.cluster.topn,file = paste(opts$out,"/",opts$tissue,"_marker_Gene_top",opts$topN,".txt",sep  =""),sep  ="\t")
  target_gene <- unique(all.cluster.topn$gene)
  if(length(target_gene) > 40){
    target_gene <- target_gene[1:40]
    print("Warning! Only 40 target genes were selected!")
  }
  #pdf(file = paste(opts$out,'/',opts$tissue,"_cluster_markers_top",opts$topN,"_DoHeatmap.pdf",sep = ""),w = 10,h = 4 * length(target_gene)/20)
  #p <- DoHeatmap(Integrate,features = target_gene,label = T) + NoLegend()
  #pritn(p)
  #while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$tissue, '_cluster_markers_top',opts$topN,'_dotPlot.pdf', sep = ""))
  p <- DotPlot(Integrate, features = target_gene) + RotatedAxis()
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$tissue, '_cluster_markers_top',opts$topN,'_VinPlot.pdf', sep = ""))
  p <- VlnPlot(Integrate, features = target_gene, stack = TRUE)
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$tissue, '_cluster_markers_top',opts$topN,'_FeaturePlot_UMAP.pdf', sep = ""),w = 16, h = 4*(length(target_gene)%/%3))
  p <- FeaturePlot(Integrate, features = target_gene, reduction = "umap", pt.size = 0.5)
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$tissue, '_cluster_markers_top',opts$topN,'_FeaturePlot_TSNE.pdf', sep = ""),w = 16, h = 4*(length(target_gene)%/%3))
  p <- FeaturePlot(Integrate, features = target_gene, reduction = "tsne", pt.size = 0.5)
  print (p)
  while (!is.null(dev.list()))  dev.off()
}
# main #
# library #
library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(presto)

# Read filtlist of RDS #
filelist <- read.csv(opts$input,header = T, sep = "\t")
print("Step1 : Read filelist of RDS Finish!")

# QC_Filter delete_doublet Intergrate #
obj.list <- list()
obj.list <- ReadRDS_list(filelist)
obj.list <- QC_Filter(obj.list,filterCounts = opts$filterCounts,filterGene = opts$filterGene , filterMito = opts$filterMito)
print(length(obj.list))
obj.list <- Do_Analyze(obj.list,method = opts$method)
# 
saveRDS(obj.list,file=paste0(opts$out,"/",opts$tissue,".Objectlist.rds"))
Integrate <- Do_Integrate_Analyze(obj.list)
print(Integrate)
print("Step2 : Filter data Finish!")

saveRDS(Integrate,file=paste0(opts$out,"/",opts$tissue,".Integrate.rds"))
print("Step6 : Save RDS Finish!")
# Find All Marker #
Integrate.markers <- Do_FindAllMarkers(Integrate)
print("Step3 : Find All Marker Finish!")

# "Marker Gene Visualize #
Do_Marker_Visualize(Integrate,marker = Integrate.markers,topN = opts$topN)
print("Step4 : Marker Visualize Finish!")

# Count the number of cells in each cluster with GF and SPF
stat <- table(Integrate@meta.data$seurat_clusters,Integrate@meta.data$mice)
write.table(stat,file = paste0(opts$out,"/",opts$tissue,'.stat.txt'),sep = "\t",quote = F)
print("Step5 : Stat Finish!")

# Save RDS #
#saveRDS(Integrate,file=paste0(opts$out,"/",opts$tissue,".Integrate.rds"))
#print("Step6 : Save RDS Finish!")

# PCA Plot#
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_PCA.pdf",sep = ""),w = 10,h = 8)
p1 <- DimPlot(Integrate,reduction = "pca",group.by = "orig.ident")
print(p1)
dev.off()
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_PCA_Group.pdf",sep = ""),w = 10, h = 8)
p2 <- DimPlot(Integrate,reduction = "pca",group.by = "mice")
print(p2)
dev.off()
print("Step7 : PCA Plot done!")

# Umap Plot #
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap.pdf",sep = ""),w = 10,h = 8)
p3 <- DimPlot(Integrate,reduction = "umap",group.by = "orig.ident")
print(p3)
dev.off()
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_Group.pdf",sep = ""),w = 10, h = 8)
p4 <- DimPlot(Integrate,reduction = "umap",group.by = "mice")
print(p4)
dev.off()
pdf(file = paste(opts$out,"/",opts$tissue,"_",opts$method,"_Umap_Cluster.pdf",sep = ""),w = 10, h = 8)
p5 <- DimPlot(Integrate,reduction = "umap",group.by = "seurat_clusters",label = T,label.size = 5)
print(p5)
dev.off()
print("Step8 : Umap Plot Finish!")
date()
