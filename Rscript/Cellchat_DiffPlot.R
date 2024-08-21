# author : lwm
# date : 2022-10-8
date()
### Get the parameter
parser = argparse::ArgumentParser(description="Script to do Diff CellChat Plot!")
parser$add_argument('--input1', help='input1 RDS after annotation of celltypes')
parser$add_argument('--input2', help='input2 RDS after annotation of celltypes')
parser$add_argument('--prefix',help = "GF_vs_SPF")
parser$add_argument('-T','--tissue',help = "tissue of mouse")
parser$add_argument('-O','--out', help='out directory')
opts = parser$parse_args()
print(opts)



# main #
#library#
library(CellChat)
library(ComplexHeatmap)
library(patchwork)
# Read RDS #
GF <- readRDS(opts$input1)
SPF <- readRDS(opts$input2)
# check celltypes overlap #
# merge #
cellchat.list <- list(GF = GF , SPF = SPF)
#cellchat <- mergeCellChat(cellchat.list ,add.names =  names(cellchat.list),cell.prefix = TRUE)
cellchat <- mergeCellChat(cellchat.list,add.names = names(cellchat.list))
outRDS <- paste(opts$out,"/",opts$prefix,"_",opts$tissue,"_cellchat.rds",sep = "")
saveRDS(cellchat,file = outRDS)
print("save RDS Done!")

# overlap celltypes for diff cellchat 
GF_celltypes <- as.vector(unique(GF@meta$celltypes))
print("GF")
print(GF_celltypes)
SPF_celltypes <- as.vector(unique(SPF@meta$celltypes))
print("SPF")
print(SPF_celltypes)
len <- length(SPF_celltypes) - length(GF_celltypes)
print(len)
if(len != 0){
	print(len)
	overlap_celltypes <- intersect(GF_celltypes,SPF_celltypes)
	print(overlap_celltypes)
	print(length(overlap_celltypes))
	subcellchat <- subsetCellChat(cellchat,idents = overlap_celltypes)
} else if(len == 0){
	print(len)
	subcellchat <- cellchat
} else{
	stop("Please check your parameter!")
}

out2 <- paste(opts$out,"/","DIFF",sep = "")
dir.create(out2)

pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_barplot.pdf', sep = ""),w = 7, h = 7)
gg1 <- compareInteractions(cellchat, show.legend = F, group  = c(1,2) , measure = "count",size.text = 20)
gg2 <- compareInteractions(cellchat, show.legend = F, group  = c(1,2) , measure = "weight",size.text =20)
p <- gg1 + gg2
print(p)
while (!is.null(dev.list()))  dev.off()

options(repr.plot.height = 7 ,repr.plot.width = 14)
pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_circle.pdf', sep = ""),w = 14, h = 7)
par(mfrow = c(1,2),xpd = TRUE)
p1 <- netVisual_diffInteraction(subcellchat,weight.scale = T)
p2 <- netVisual_diffInteraction(subcellchat,weight.scale = T, measure = "weight")
print(p1)
print(p2)
while (!is.null(dev.list()))  dev.off()

pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_heatmap.pdf', sep = ""),w = 14, h = 7)
par(mfrow = c(1,2))
h1 <- netVisual_heatmap(subcellchat)
h2 <- netVisual_heatmap(subcellchat,measure = "weight")
p <- h1 + h2
print(p)
while (!is.null(dev.list()))  dev.off()

pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_number_interaction_circle.pdf', sep = ""),w = 14, h = 7)
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cellchat.list,attribute = c("ident","count"))
for (i in 1:length(cellchat.list)){
     netVisual_circle(cellchat.list[[i]]@net$count,weight.scale  = T , label.edge = F , edge.weight.max = weight.max[2] , edge.width.max = 12,
                    title.name = paste0("Number of interactions - " , names(cellchat.list)[i]))
}
while (!is.null(dev.list()))  dev.off()

# sub set #
#pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_subcelltypes_number_interaction_circle.pdf', sep = ""),w = 14, h = 7)
#par(mfrow  = c(1,2),xpd = TRUE)
#s.cell <- c("Neutrophil(Chil3 high)","Neutrophil(Mpo high)","Neutrophil(Mmp8 high)","Neutrophil(Ltf high)","Neutrophil(Elane high)")
#counts1 <- cellchat.list[[1]]@net$count[s.cell,s.cell]
#counts2 <- cellchat.list[[2]]@net$count[s.cell,s.cell]
#p1 <- weight.max <- max(max(counts1),max(counts2))
#p2 <- netVisual_circle(counts1,weight.scale  = T, label.edge = T, edge.weight.max = weight.max , edge.width.max = 12 , 
#                title.name = paste0("Number of interactions- ", names(cellchat)[[1]]))
#netVisual_circle(counts2,weight.scale  = T, label.edge = T, edge.weight.max = weight.max , edge.width.max = 12 , 
#                title.name = paste0("Number of interactions- ", names(cellchat)[[2]]))
#print(p1)
#print(p2)
#while (!is.null(dev.list()))  dev.off()

#pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_rankNet.pdf', sep = ""),w = 14, h = 7)
#gg1 <- rankNet(cellchat, mode = "comparison",stacked = T,  do.stat = TRUE)
#gg2 <- rankNet(cellchat, mode = "comparison",stacked = F, do.stat  = TRUE)
#p <- gg1 + gg2
#print(p)
#while (!is.null(dev.list()))  dev.off()

# structual #
options(repr.plot.height = 7 ,repr.plot.width = 9)
pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_structural_clustering.pdf', sep = ""),w = 7, h = 7)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat,type = "structural",umap.method = "uwot")
cellchat <- netClustering(cellchat,type = "structural")
netVisual_embeddingPairwise(cellchat,type = "structural", label.size = 5)
while (!is.null(dev.list()))  dev.off()

pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_structural_pathways.pdf', sep = ""),w = 7, h = 7)
p <- rankSimilarity(cellchat,type = "structural",font.size = 15) + ggtitle("Structural similarity of pathway")
print(p)
dev.off()

options(repr.plot.height = 7 ,repr.plot.width = 14)
pathway.union <- union(cellchat.list[[1]]@netP$pathways, cellchat.list[[2]]@netP$pathways)
pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_all_signalingRole.pdf', sep = ""),w = 14, h = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.list[[1]], pattern = "all" , signaling  = pathway.union,title = names(cellchat.list)[1],width = 8,height = 14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.list[[2]], pattern = "all" , signaling  = pathway.union,title = names(cellchat.list)[2],width = 8,height = 14)
draw(ht1 + ht2, ht_gap = unit(0.5,"cm"))
while (!is.null(dev.list()))  dev.off()

pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_incoming_signalingRole.pdf', sep = ""),w = 14, h = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.list[[1]], pattern = "incoming" , signaling  = pathway.union,title = names(cellchat.list)[1],width = 8,height = 14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.list[[2]], pattern = "incoming" , signaling  = pathway.union,title = names(cellchat.list)[2],width = 8,height = 14)
draw(ht1 + ht2, ht_gap = unit(0.5,"cm"))
while (!is.null(dev.list()))  dev.off()

pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_Diff_outgoing_signalingRole.pdf', sep = ""),w = 14 , h = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.list[[1]], pattern = "outgoing" , signaling  = pathway.union,title = names(cellchat.list)[1],width = 8,height = 14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.list[[2]], pattern = "outgoing" , signaling  = pathway.union,title = names(cellchat.list)[2],width = 8,height = 14)
draw(ht1 + ht2, ht_gap = unit(0.5,"cm"))
while (!is.null(dev.list()))  dev.off()

#for(i in 1:length(pathway.union)){
#	pathways.show <- pathway.union[1]
#	out3 <- paste(out2,'/',pathways.show,sep = "")
#	dir.create(out3)
#	weight.max <- getMaxWeight(cellchat.list, slot.name  = c("netP"), attribute = pathways.show)
#	pdf(file = paste(out3,'/', opts$prefix,'_',opts$tissue, '_' ,pathways.show,'_diff_pathway_circle.pdf', sep = ""),w = 14 , h = 7)
#	par(mfrow = c(1,2) , xpd = TRUE)
#	for( i in 1:length(cellchat.list)){
#        	netVisual_aggregate(cellchat.list[[i]],signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1],edge.width.max = 10,signaling.name = paste(pathways.show, names(cellchat.list)[i]))	
#	}
#	while (!is.null(dev.list()))  dev.off()
#	pdf(file = paste(out3,'/', opts$prefix,'_',opts$tissue, '_' ,pathways.show,'_diff_pathway_Reds.pdf', sep = ""),w = 7, h = 7)
#	par(mfrow = c(1,2) , xpd = TRUE)
#	ht <- list()
#	for( i in 1:length(cellchat.list)){
#       		ht[i] <- netVisual_heatmap(cellchat.list[[i]],signaling = pathways.show, color.heatmap = "Reds",
#                        title.name = paste(pathways.show, "signaling",names(cellchat.list)[i]))
#	}
#	ComplexHeatmap::draw(ht[[1]] + ht[[2]] , ht_gap = unit(0.5,"cm"))
#	while (!is.null(dev.list()))  dev.off()
#	pdf(file = paste(out3,'/', opts$prefix,'_',opts$tissue, '_' ,pathways.show,'_diff_pathway_chord.pdf', sep = ""))
#	par(mfrow = c(1,2) , xpd = TRUE)
#	for( i in 1:length(cellchat.list)){
#        	p <- netVisual_aggregate(cellchat.list[[i]],signaling = pathways.show, layout = "chord", pt.title = 3, title.space = 0.05, vertex.label.cex = 0.6,signaling.name = paste(pathways.show, names(cellchat.list)[i]))
#		print(p)
#	}
#	while (!is.null(dev.list()))  dev.off()	
#	stop("test")
}

#options(repr.plot.height = 7 ,repr.plot.width = 7)
#pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_' ,pathways.show,'_diff_pathway_bubble.pdf', sep = ""),w = 10, h = 8)
#print(levels(cellchat@idents$joint))
#p <- netVisual_bubble(cellchat, sources.use = c(4,5,7),targets.use = c(1,2,3,6), comparison = c(1,2),angle.x = 45)
#print(p)
#while (!is.null(dev.list()))  dev.off()

#options(repr.plot.height = 7 ,repr.plot.width = 14)
#pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_' ,pathways.show,'_diff_pathway_signaling_bubble.pdf', sep = ""),w = 10, h = 8)
#p1 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1,2), angle.x = 45,max.dataset = 1,
#                       title.name = paste("Increased signaling in ",pathways.show,sep = ""), remove.isolate = T)
#p2 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6),comparison = c(1,2), angle.x = 45,max.dataset = 2,
#                       title.name = paste("Decreased signaling in ",pathways.show,sep = ""), remove.isolate = T)
#pc <- p1 + p2
#print(pc)
#while (!is.null(dev.list()))  dev.off()

#pdf(file = paste(out2,'/', opts$prefix,'_',opts$tissue, '_' ,pathways.show,'_diff_chord_gene.pdf', sep = ""),w = 10, h = 8)
#par(mfrow = c(1,2),xpd = TRUE)
#for( i in 1:length(cellchat.list)){
#        netVisual_chord_gene(cellchat.list[[i]],sources.use = c(4,5), targets.use = c(1,2,3,6), lab.cex = 0.6 , legend.pos.x = 10, legend.pos.y = 20,title.name = paste0("Signaling from Treg - " , names(cellchat.list)[i]))
#}
#while (!is.null(dev.list()))  dev.off()
# save RDS #
saveRDS(cellchat,file = paste0(out2,"/",opts$tissue,"_diffcellchat.rds"))
print("Step : save RDS finish!")
print("haha")

