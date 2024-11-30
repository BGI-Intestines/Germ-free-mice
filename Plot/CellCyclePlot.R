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
library(tibble)
################################################################################
######## Args
################################################################################
parser = argparse::ArgumentParser(description = 'Script for a Umap plot of celltypes')
parser$add_argument('--input_GF', required = T, help = '')
parser$add_argument('--input_SPF', required = T, help = '')
parser$add_argument('--output', required = T, help = 'output pathway')
parser$add_argument('--annotation', default = NULL, help = 'cellcycle annotation')
parser$add_argument('--filtercluster', default = NULL, help = 'filtercluster')
parser$add_argument('--clustercolour', default = NULL, help = 'ClusterColour')
parser$add_argument('--cellcyclecolour', default = NULL, help = 'cellcyclecolour')
parser$add_argument('--labelfont', default = "sans",  choices = c("serif","sans"), help = 'labelfont')
parser$add_argument('--labelfontsize', default = 6,  type = "double", help = 'labelfontsize')
parser$add_argument('--plottitlefont', default = "sans",  choices = c("serif","sans"), help = 'plottitlefont')
parser$add_argument('--plottitlefontsize', default = 20,  type = "double", help = 'plottitlefontsize')
parser$add_argument('--legend', default = TRUE, type = "logical", help = 'Legend')
parser$add_argument('--legendnumsize', default = 4, type = "double", help = 'Legend')
parser$add_argument('--legendnumfont', default = "sans",  choices = c("serif","sans"), help = 'labelfont')
parser$add_argument('--legendlabelfontsize', default = 5, type = "double", help = 'Legend')
parser$add_argument('--legendlabelfont', default = "sans",  choices = c("serif","sans"), help = 'labelfont')
parser$add_argument('--height', default = 10,  type = "double", help = 'height')
parser$add_argument('--width', default = 20,  type = "double", help = 'width')
parser$add_argument('--legendheight', default = 4,  type = "double", help = 'legendheight')
parser$add_argument('--legendwidth', default = 5,  type = "double", help = 'legendwidth')

opts = parser$parse_args()
# print(opts)

GF_cellcycle_file <- opts$input_GF
SPF_cellcycle_file <- opts$input_SPF
OutputPath <- opts$output
cellcycle_annotation <- opts$annotation
FilterCluster <- opts$filtercluster
if (!is.null(FilterCluster)){
  FilterCluster <- eval(str2lang(FilterCluster))
}
ClusterColour <- opts$clustercolour
if (!is.null(ClusterColour)){
  ClusterColour <- eval(str2lang(ClusterColour))
}
CellCycleColour <- opts$cellcyclecolour
if (!is.null(CellCycleColour)){
  CellCycleColour <- eval(str2lang(CellCycleColour))
}

labelfont = opts$labelfont
labelfontsize = opts$labelfontsize
plottitlefont = opts$plottitlefont
plottitlefontsize = opts$plottitlefontsize

Legend <- opts$legend
legendnumfont <- opts$legendnumfont
legendnumsize <- opts$legendnumsize
legendlabelfont <- opts$legendlabelfont
legendlabelfontsize <- opts$legendlabelfontsize
plot_height = opts$height
plot_width = opts$width
legendheight = opts$legendheight
legendwidth = opts$legendwidth

################################################################################
######## Function 
################################################################################
CellCyclePlot <- function(df, groupname){
  df$Overall <- apply(df,1,sum)
  SumAll <- sum(df$Overall)
  df$gap <- SumAll*0.01
  df$Overall <- NULL
  df <- rownames_to_column(df, var = "Cluster")
  df_long <- reshape2::melt(df, id.var="Cluster", variable.name="CellCycle", value.name = "Value")
  
  df_long$Cluster <- as.factor(df_long$Cluster)
  df_long$Cluster <- factor(df_long$Cluster, levels=sort(as.numeric(levels(df_long$Cluster))))
  df_long$CellCycle <- factor(df_long$CellCycle, levels=c("G1","G2M","S","gap"))
  
  df_long <- df_long[order(df_long[,"Cluster"],df_long[,"CellCycle"]),]
  
  df_long$fraction = df_long$Value / SumAll
  df_long$ymax = cumsum(df_long$fraction)
  df_long$ymin = c(0, head(df_long$ymax, n = -1))
  
  write.csv(df_long, paste0("df_long_",groupname,".csv"),row.names = F)
  
  df_long_center <- df_long[which(df_long$CellCycle!="gap"),c("Cluster","ymax","ymin")]
  df_center <- reshape2::melt(df_long_center, id.var="Cluster")
  df_center_y <- df_center %>% group_by(Cluster) %>% dplyr::summarise(ymin=min(value),ymax=max(value)) %>% as.data.frame()
  df_center_y$ymean <- (df_center_y$ymax + df_center_y$ymin)/2
  
  p <- ggplot() +
    geom_rect(data = df_long[which(df_long$CellCycle!="gap"),], aes(fill = Cluster, ymax = ymax, ymin = ymin, xmax = 7, xmin = 3),show.legend = F) +
    geom_rect(data = df_long[which(df_long$CellCycle=="gap"),], aes( ymax = ymax, ymin = ymin, xmax = 7, xmin = 3),fill = "black",show.legend = F) +
    geom_rect(data = df_long, aes(fill = CellCycle, ymax = ymax, ymin = ymin, xmax = 8, xmin = 7.1),show.legend = F) +
    # geom_point(data = df_center_y, aes(x=rep(5,nrow(df_center_y)),y=ymean), alpha=0, size=0.01)+
    geom_text(data = df_center_y, aes(x=rep(9,nrow(df_center_y)),y=ymean,label=Cluster), fontface="bold",
              color = "white",size = labelfontsize, family = labelfont)+
    annotate("text", x = 0.1, y = median(df_long$ymax), label = groupname,
             color="white",size = plottitlefontsize, family = plottitlefont, fontface="bold", angle=0)+
    scale_fill_manual(values=mycolor2)+
    coord_polar(theta = "y") +
    xlim(c(0, 10)) +
    # theme_light() +
    dark_theme_bw()+
    labs(x="",y="")+
    theme(panel.grid=element_blank()) + 
    theme(axis.text=element_blank()) + 
    theme(axis.ticks=element_blank()) + 
    theme(panel.border=element_blank()) 
  return(p)
}

df_rename <- function(df, annote_ct){
  df_merged <- merge(df, annote_ct, by.x="row.names", by.y="ClusterName")
  rownames(df_merged) <- df_merged$Cluster
  df_merged <- df_merged[order(df_merged$Cluster),]
  df_merged$Cluster <- NULL
  df_merged$Row.names <- NULL
  return(df_merged)
}

df_annot <- function(df){
  df_a_s <- NULL
  df <- rownames_to_column(df, var = "ClusterName")
  df <- rownames_to_column(df, var = "Cluster")
  ct_annotation <- df[,1:2]
  df_selected <- df %>% select(-c("ClusterName")) %>% column_to_rownames("Cluster")
  df_a_s[[1]] <- ct_annotation
  df_a_s[[2]] <- df_selected
  return(df_a_s)
}

################################################################################
######## Main 
################################################################################
dir.create(OutputPath)
setwd(OutputPath)

GF_Phase <- read.table(GF_cellcycle_file,sep="\t", header = T, row.names = 1)
SPF_Phase <- read.table(SPF_cellcycle_file,sep="\t", header = T, row.names = 1)

if (!is.null(FilterCluster)){
  GF_Phase <- GF_Phase[which(! rownames(GF_Phase) %in% FilterCluster),]
  SPF_Phase <- SPF_Phase[which(! rownames(SPF_Phase) %in% FilterCluster),]
}


if (!is.null(cellcycle_annotation)){
  cellcycle_annotation <- read.table(cellcycle_annotation, sep="\t", header = T)
  GF_Phase <- df_rename(GF_Phase, cellcycle_annotation)
  SPF_Phase <- df_rename(SPF_Phase, cellcycle_annotation)
}else{
  GF_annot_res <- df_annot(GF_Phase)
  GF_Phase <- GF_annot_res[[2]]
  SPF_annot_res <- df_annot(SPF_Phase)
  SPF_Phase <- SPF_annot_res[[2]]
  cellcycle_annotation <- SPF_annot_res[[1]]
}

nCluster = nrow(GF_Phase)

if (is.null(ClusterColour)){
  mycolor = colorRampPalette(brewer.pal(9, "Set1"))(nCluster)
  names(mycolor) <- rownames(GF_Phase)
}else{
  if (is.null(names(ClusterColour))){
    mycolor <- ClusterColour
    names(mycolor) <- rownames(GF_Phase)
  }else{
    mycolor <- ClusterColour
  }
}

if (is.null(CellCycleColour)){
  CellCycleColour <- c("#0073C2FF", "#EFC000FF", "#CD534CFF")
}
mycolor2 <- c(mycolor, "G1"=CellCycleColour[1], "G2M"=CellCycleColour[2], "S"=CellCycleColour[3],"gap"="black")

ccp_gf <- CellCyclePlot(GF_Phase, "GF")
ccp_spf <- CellCyclePlot(SPF_Phase, "SPF")

ccp_group <- ggarrange(ccp_gf, ccp_spf, 
                       ncol = 2, nrow = 1)
ggsave("CellCyclePlot_group.pdf",ccp_group,device='pdf',width = plot_width, height = plot_height)

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
    scale_color_manual(values=mycolor2)+
    theme(legend.position = "none")
  ggsave("CellCyclePlot_legend.pdf",p_lgd,device='pdf',width = legendwidth, height = legendheight)
  
  df_legend2 <- data.frame(Cluster=c("G1","G2M","S"))
  df_legend2$Cluster <- as.factor(df_legend2$Cluster)
  df_legend2$y <- c(rep(0,3))
  df_legend2$x <- c(seq(1,3))
  
  p_lgd2 <- ggplot(df_legend2, aes(x=x, y=y))+
    geom_rect(data = df_legend2, aes(fill = as.factor(Cluster), ymax = y+0.02, ymin = y-0.02, xmax = x+0.05, xmin = x-0.05),show.legend = F)+
    labs(x='',y='',title="")+
    geom_text(aes(x=x+0.1, y=y,label=Cluster), fontface="bold",color = "white",size = legendlabelfontsize, family = legendlabelfont, hjust = 0)+
    geom_text(aes(x=1, y=0.1, label="Cell Cycle"), fontface="bold",color = "white",size = 6, family = legendlabelfont, hjust = 0)+
    # theme_map()+
    dark_theme_void()+
    scale_x_continuous(limits = c(-0.5, 4))+
    scale_y_continuous(limits = c(-0.5, 0.5))+
    theme(axis.text = element_blank())+
    theme(axis.ticks = element_blank())+
    theme(panel.border = element_blank())+
    theme(axis.title = element_blank())+
    scale_fill_manual(values=CellCycleColour)+
    theme(legend.position = "none")
  ggsave("CellCyclePlot_legend2.pdf",p_lgd2,device='pdf',width = legendwidth, height = legendheight)
}



