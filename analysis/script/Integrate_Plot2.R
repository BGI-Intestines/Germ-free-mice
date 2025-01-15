## author : lwm
# date : 2022-9-10
date()
parser = argparse::ArgumentParser(description="Script to QC and Cluster scRNA data")
parser$add_argument('-I','--input', help='input filelist of rds and sample ID')
parser$add_argument('-T','--tissue',help = "tissue of mouse")
parser$add_argument('-M','--marker',help = "a file of ALL marker ")
parser$add_argument('--topN', default = 1, type = 'integer', help = 'Top N marker genes in each cluster')
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)

# main #
library(data.table)
library(dplyr)
library(Seurat)
library(DoubletFinder)
library(presto)
library(ggplot2)

# Read RDS #
Integrate <- readRDS(opts$input)
print("Step 0 : Read RDS Finish!")
marker <- read.csv(opts$marker,sep = "\t",header = T)
print(head(marker))
print("Stop 1 : Read marker file Finish!")

# marker visualize#
Do_Marker_Visualize <- function(Integrate,marker,topN){
  all.cluster.topn <- marker %>% group_by(cluster) %>% top_n(n = topN,wt = avg_log2FC)
  target_gene <- unique(all.cluster.topn$gene)
  if(length(target_gene) > 40){
    target_gene <- target_gene[1:40]
    print("Warning! Only 40 target genes were selected!")
  }
  pdf(file = paste(opts$out,'/', opts$tissue, '_cluster_markers_top',opts$topN,'_dotPlot.pdf', sep = ""))
  p <- DotPlot(Integrate, features = target_gene) + RotatedAxis() + theme(axis.text = element_text(size = 6))
  print (p)
  while (!is.null(dev.list()))  dev.off()
  pdf(file = paste(opts$out,'/', opts$tissue, '_cluster_markers_top',opts$topN,'_VinPlot.pdf', sep = ""))
  p <- VlnPlot(Integrate, features = target_gene, stack = TRUE) + theme(axis.title = element_text(size = 6))
  print (p)
  while (!is.null(dev.list()))  dev.off()
}
 
Do_Marker_Visualize(Integrate = Integrate , marker = marker , topN = opts$topN)
print("Step 2 : marker Visualize Finish!")
date()
