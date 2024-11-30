################################################################################
######## Import packages 
################################################################################
library(argparse)
library(Seurat)
library(ggdark)
library(ggsci)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
################################################################################
######## Args
################################################################################
parser = argparse::ArgumentParser(description = 'Script for a Umap plot of celltypes')
parser$add_argument('--input', required = T, help = 'a RDS after annotation ')
parser$add_argument('--output', required = T, help = 'output Path')
parser$add_argument('--annotation', default = NULL, help = 'cellcycle annotation file')
parser$add_argument('--saveannotation', default = FALSE, type = "logical", help = 'save annotation')
parser$add_argument('--umapcluster', choices = c("Cluster","Group","Tissue","Platform"), help = 'umapCluster')
parser$add_argument('--filtercluster', default = NULL, help = 'filtercluster')
parser$add_argument('--filtercluster_var', default = "Cluster", choices = c("Cluster","Group","Tissue","Platform"), help = 'filtercluster_var')
parser$add_argument('--clustercolour', default = NULL, help = 'ClusterColour')
parser$add_argument('--umapgroup', default = NULL, choices = c("Cluster","Group","Tissue","Platform"), help = 'umapGroup')
parser$add_argument('--specificgroup', default = NULL, help = 'specificgroup')
parser$add_argument('--labelfont', default = "sans",  choices = c("serif","sans"), help = 'labelfont')
parser$add_argument('--labelfontsize', default = 6,  type = "double", help = 'labelfontsize')
parser$add_argument('--axistitlefont', default = "sans",  choices = c("serif","sans"), help = 'axistitlefont')
parser$add_argument('--axistitlefontsize', default = 4,  type = "double", help = 'axistitlefontsize')
parser$add_argument('--axislinesize', default = 1,  type = "double", help = 'axislinesize')
parser$add_argument('--xaxislengthfold', default = 0.3,  type = "double", help = 'xaxislengthfold')
parser$add_argument('--yaxislengthfold', default = 0.3,  type = "double", help = 'yaxislengthfold')
parser$add_argument('--pointsize', default = 0.11, type = "double", help = 'pointsize')
parser$add_argument('--pointalpha', default = 1, type = "double", help = 'pointalpha')
parser$add_argument('--plottitlefont', default = "sans",  choices = c("serif","sans"), help = 'plottitlefont')
parser$add_argument('--plottitlefontsize', default = 30,  type = "double", help = 'plottitlefontsize')
parser$add_argument('--legend', default = FALSE, type = "logical", help = 'Legend')
parser$add_argument('--legendnumsize', default = 4, type = "double", help = 'Legend')
parser$add_argument('--legendnumfont', default = "sans",  choices = c("serif","sans"), help = 'labelfont')
parser$add_argument('--legendlabelfontsize', default = 5, type = "double", help = 'Legend')
parser$add_argument('--legendlabelfont', default = "sans",  choices = c("serif","sans"), help = 'labelfont')
parser$add_argument('--nrow_plot', default = NULL, type = 'integer', help = 'nrow_plot')
parser$add_argument('--ncol_plot', default = NULL, type = 'integer', help = 'ncol_plot')
parser$add_argument('--height', default = 10,  type = "double", help = 'height')
parser$add_argument('--width', default = 20,  type = "double", help = 'width')
parser$add_argument('--legendheight', default = 4,  type = "double", help = 'legendheight')
parser$add_argument('--legendwidth', default = 5,  type = "double", help = 'legendwidth')

opts = parser$parse_args()
# print(opts)

Seurat_rds_file <- opts$input
OutputPath <- opts$output
cellcycle_annotation <- opts$annotation
saveannotation <- opts$saveannotation
FilterCluster <- opts$filtercluster
if (!is.null(FilterCluster)){
  FilterCluster <- eval(str2lang(FilterCluster))
}
filtercluster_var <- opts$filtercluster_var
umapCluster <- opts$umapcluster
ClusterColour <- opts$clustercolour
if (!is.null(ClusterColour)){
  ClusterColour <- eval(str2lang(ClusterColour))
}
umapGroup <- opts$umapgroup
specificgroup <- opts$specificgroup
if (!is.null(specificgroup)){
  specificgroup <- eval(str2lang(specificgroup))
}
labelfont = opts$labelfont
labelfontsize = opts$labelfontsize
axistitlefont = opts$axistitlefont
axistitlefontsize = opts$axistitlefontsize
axislinesize = opts$axislinesize
xaxislengthfold = opts$xaxislengthfold
yaxislengthfold = opts$yaxislengthfold
pointsize = opts$pointsize
pointalpha = opts$pointalpha
plottitlefont = opts$plottitlefont
plottitlefontsize = opts$plottitlefontsize
Legend <- opts$legend
legendnumfont <- opts$legendnumfont
legendnumsize <- opts$legendnumsize
legendlabelfont <- opts$legendlabelfont
legendlabelfontsize <- opts$legendlabelfontsize
nrow_plot = opts$nrow_plot
ncol_plot = opts$ncol_plot
plot_height = opts$height
plot_width = opts$width
legendheight = opts$legendheight
legendwidth = opts$legendwidth

################################################################################
######## Function 
################################################################################
GetGroupCenter <- function(df_umap, group){
  GroupCenter <- df_umap %>% group_by(eval(str2lang(group))) %>% summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2)) %>% as.data.frame()
  colnames(GroupCenter)[1] <- "GroupCenter"
  return(GroupCenter)
}

UmapPlot <- function(df, cluster, df_center, labelfontsize, axistitlefontsize, pointsize, pointalpha){
  p <- ggplot() + 
    geom_point(df, mapping=aes(x = UMAP_1, y = UMAP_2, color = eval(str2lang(cluster))), size = pointsize, alpha = pointalpha) +
    dark_theme_void()+
    scale_x_continuous(limits = c(xmin-2, xmax))+
    scale_y_continuous(limits = c(ymin-2, ymax))+
    theme(axis.text = element_blank())+
    theme(axis.ticks = element_blank())+
    theme(panel.border = element_blank())+
    theme(axis.title = element_blank())+
    scale_color_manual(values=mycolor)+
    geom_segment(aes(x=xmin, y=ymin, xend=xmin+abs(xmin)*xaxislengthfold, yend=ymin), arrow = arrow(length=unit(0.2, "cm")), color="white", linewidth=axislinesize)+
    geom_segment(aes(x=xmin, y=ymin, xend=xmin, yend=ymin+abs(ymin)*yaxislengthfold), arrow = arrow(length=unit(0.2, "cm")), color="white", linewidth=axislinesize)+
    annotate("text", x = xmin + 1.5, y = ymin - 1, label = "UMAP_1",
             color="white",size = axistitlefontsize, fontface="bold", family = axistitlefont) +
    annotate("text", x = xmin - 1, y = ymin + 1.5, label = "UMAP_2",
             color="white",size = axistitlefontsize, fontface="bold" ,angle=90, family = axistitlefont) +
    # geom_point(mapping=aes(x=UMAP_1, y=UMAP_2), data = df_center, size = 1, alpha = 1) +
    # geom_text_repel(data = df_center[which(df_center$GroupCenter %in% unique(df[,cluster])),], mapping=aes(x = UMAP_1, y = UMAP_2, label=GroupCenter), fontface="bold",
    #                 point.padding=unit(2, "lines"),color = "white",size = labelfontsize)+
    geom_text(data = df_center[which(df_center$GroupCenter %in% unique(df[,cluster])),], mapping=aes(x = UMAP_1, y = UMAP_2, label=GroupCenter), 
              fontface="bold",color = "white",size = labelfontsize, family = labelfont)+
    theme(legend.position = "none")
  return(p)
}

MultiUmaplots <- function(df,group){
  umap_groupname_list <- NULL
  for (i in seq_len(nlevels(df[,group]))){
    groupnames <- levels(df[,group])
    df_tmp <- df[which(df[,group] == groupnames[i]),]
    umap_groupname_list[[i]] <-  UmapPlot(df_tmp, umapCluster, ct_center, labelfontsize, axistitlefontsize, pointsize, pointalpha)
  }
  return(umap_groupname_list)
}

################################################################################
######## Main 
################################################################################
dir.create(OutputPath)
setwd(OutputPath)

seuratobj <- readRDS(Seurat_rds_file)

df_umap_ct = seuratobj@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% 
  # cbind(Cluster = seuratobj@meta.data$seurat_clusters) %>%
  cbind(Cluster = seuratobj@meta.data$clusters) %>%
  cbind(Group = seuratobj@meta.data$mice) 
if (is.null(seuratobj@meta.data$tissue2)){
  df_umap_ct$Tissue <- seuratobj@meta.data$tissue
}else{
  df_umap_ct$Tissue <- seuratobj@meta.data$tissue2
}
if (!is.null(seuratobj@meta.data$platform)){
  df_umap_ct$Platform <- seuratobj@meta.data$platform
  df_umap_ct$Platform <- as.factor(df_umap_ct$Platform)
}

df_umap_ct$Cluster <- as.factor(df_umap_ct$Cluster)
df_umap_ct$Cluster <- factor(df_umap_ct$Cluster, levels=sort(as.numeric(levels(df_umap_ct$Cluster))))
df_umap_ct$Group <- as.factor(df_umap_ct$Group)
df_umap_ct$Tissue <- as.factor(df_umap_ct$Tissue)

if (!is.null(FilterCluster)){
  if (filtercluster_var == umapCluster){
    df_umap_ct <- df_umap_ct[which(!df_umap_ct[,umapCluster] %in% FilterCluster),]
  }else{
    df_umap_ct <- df_umap_ct[which(!df_umap_ct[,filtercluster_var] %in% FilterCluster),]
  }
}

xmin <- min(df_umap_ct$UMAP_1)-2
xmax <- max(df_umap_ct$UMAP_1)+2
ymin <- min(df_umap_ct$UMAP_2)-2
ymax <- max(df_umap_ct$UMAP_2)+2

ct_center <- GetGroupCenter(df_umap_ct, umapCluster)

nCluster = length(unique(df_umap_ct[,umapCluster]))

if (is.null(ClusterColour)){
  mycolor = colorRampPalette(brewer.pal(9, "Set1"))(nCluster)
  names(mycolor) <- levels(df_umap_ct[,umapCluster])[sort(unique(df_umap_ct[,umapCluster]))]
}else{
  if (is.null(names(ClusterColour))){
    mycolor <- ClusterColour
    names(mycolor) <- levels(df_umap_ct[,umapCluster])[sort(unique(df_umap_ct[,umapCluster]))]
  }else{
    mycolor <- ClusterColour
  }
}


if (is.null(specificgroup)){
  if (is.null(umapGroup)){
    umap_overall <- UmapPlot(df_umap_ct, umapCluster, ct_center, labelfontsize, axistitlefontsize, pointsize, pointalpha)
    ggsave(paste0("umap_overall_",umapCluster,".pdf"),umap_overall,device='pdf',width = plot_width, height = plot_height)
  }else{
    umap_groupname_list <- MultiUmaplots(df_umap_ct, umapGroup)
    umap_group <- ggarrange(plotlist=umap_groupname_list,
                            labels = levels(df_umap_ct[,umapGroup]),
                            font.label = list(size = plottitlefontsize, color = "white", face = "bold", family = plottitlefont),
                            ncol = ncol_plot, nrow = nrow_plot)
    ggsave(paste0("umap_", umapCluster, "_", umapGroup, ".pdf"), umap_group,device='pdf',width = plot_width, height = plot_height)
  }
}else{
  df_umap_ct <- df_umap_ct[which(df_umap_ct[,specificgroup[1]] == specificgroup[2]),]
  if (is.null(umapGroup)){
    umap_overall <- UmapPlot(df_umap_ct, umapCluster, ct_center, labelfontsize, axistitlefontsize, pointsize, pointalpha)
    ggsave(paste0(specificgroup[2],"_umap_overall_",umapCluster,".pdf"),umap_overall,device='pdf',width = plot_width, height = plot_height)
  }else{
    umap_groupname_list <- MultiUmaplots(df_umap_ct, umapGroup)
    umap_group <- ggarrange(plotlist=umap_groupname_list,
                            labels = levels(df_umap_ct[,umapGroup]),
                            font.label = list(size = plottitlefontsize, color = "white", face = "bold", family = plottitlefont),
                            ncol = ncol_plot, nrow = nrow_plot)
    ggsave(paste0(specificgroup[2],"_umap_", umapCluster, "_", umapGroup, ".pdf"),umap_group,device='pdf',width = plot_width, height = plot_height)
  }
}

if (is.null(cellcycle_annotation)){
  ct_anno <- unique(data.frame(Cluster=seuratobj@meta.data$cluster, ClusterName=seuratobj@meta.data$celltypes))
  ct_anno$Cluster <- factor(ct_anno$Cluster, levels=sort(as.numeric(levels(ct_anno$Cluster))))
  cellcycle_annotation <- ct_anno[order(ct_anno$Cluster, decreasing = F),]
  if (!is.null(FilterCluster)){
    cellcycle_annotation <- cellcycle_annotation[which(!cellcycle_annotation$Cluster %in% FilterCluster),]
  }
  if (saveannotation){
    write.table(cellcycle_annotation, "cellcycle_annotation.txt", sep = "\t", quote = F, row.names = F)
  }
}else{
  cellcycle_annotation <- read.table(cellcycle_annotation, sep="\t", header = T)
  cellcycle_annotation$Cluster <- as.factor(cellcycle_annotation$Cluster)
  cellcycle_annotation$Cluster <- factor(cellcycle_annotation$Cluster, levels=sort(as.numeric(levels(cellcycle_annotation$Cluster))))
  cellcycle_annotation <- cellcycle_annotation[order(cellcycle_annotation$Cluster),]
}

if (Legend){
  if (!is.null(FilterCluster)){
    cellcycle_annotation <- cellcycle_annotation[which(!cellcycle_annotation$Cluster %in% FilterCluster),]
  }
  
  xend <- ceiling(nCluster/15)
  xeach <- ceiling(nCluster/xend)
  xlist <- sort(rep(c(1:xend), each=xeach, length.out = nCluster))
  df_xlist <- table(xlist) %>% as.data.frame()
  ylist <- rep(c(max(df_xlist$Freq):1),nrow(df_xlist), length.out=nCluster)
  
  df_legend <- cellcycle_annotation
  df_legend$x <- xlist
  df_legend$y <- ylist
  
  p_lgd <- ggplot(df_legend, aes(x=x, y=y))+
    geom_point(alpha=1, size=6,aes(colour=as.factor(Cluster)))+
    labs(x='',y='',title="")+
    geom_text(aes(x=x, y=y,label=Cluster), fontface="bold",color = "white",size = legendnumsize, family = legendnumfont)+
    geom_text(aes(x=x+xend*0.03, y=y,label=ClusterName), fontface="bold",color = "white",size = legendlabelfontsize, family = legendlabelfont, hjust = 0)+
    dark_theme_void()+
    scale_x_continuous(limits = c(1, xend+1))+
    # scale_y_continuous(limits = c(ymin-2, ymax))+
    theme(axis.text = element_blank())+
    theme(axis.ticks = element_blank())+
    theme(panel.border = element_blank())+
    theme(axis.title = element_blank())+
    scale_color_manual(values=mycolor)+
    theme(legend.position = "none")
  ggsave("umap_legend.pdf",p_lgd,device='pdf',width = legendwidth, height = legendheight)
  
}

