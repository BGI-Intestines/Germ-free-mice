library(ggsci)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(reshape2)
library(pheatmap)

#dir.create(OutputPath)
#setwd(OutputPath)

wd <- '/Users/dengysh/Work/Projects/SpatialTranscriptomics/GermFreeMiceProject/GetFigures/GOPlots/tissue_GO/'
# list.files(wd)

getdf <- function(df){
  df <- read.delim(df,sep="\t", header = T)
  generatio <- sapply(df$GeneRatio,function(x){eval(str2lang(x))})
  df$generatio <- round(generatio, 3)
  return(df)
}

all_file <- '/Users/dengysh/Work/Projects/SpatialTranscriptomics/GermFreeMiceProject/GetFigures/GOPlots/tissue_GO/Epithelial_DEG_ALL_GO.txt'
df_all <- getdf(all_file)

celist <- c("EC (Car1 high)","EC (Hmgcs2 high)","Goblet1")
colist <- c("EC (Saa1 high)","EC (Cmss1 high)","Goblet2")
illist <- c("EC (Apoa4 high)","EC (Reg3g high)","EC (Olfm4 high)")

co_aph <- paste0(wd,"cecum_EP_EC (Cmss1 high)_GO.txt")
co_csh <- paste0(wd,"cecum_EP_EC (Saa1 high)_GO.txt")
co_gl <- paste0(wd,"cecum_EP_Goblet2_GO.txt")

ce_car1h <- paste0(wd,"colon_EP_EC (Car1 high)_GO.txt")
ce_hch <- paste0(wd,"colon_EP_EC (Hmgcs2 high)_GO.txt")
ce_gl <- paste0(wd,"colon_EP_Goblet1_GO.txt")

il_aph <- paste0(wd,"ileum_EP_EC (Apoa4 high)_GO.txt")
il_omh <- paste0(wd,"ileum_EP_EC (Olfm4 high)_GO.txt")
il_rgh <- paste0(wd,"ileum_EP_EC (Reg3g high)_GO.txt")

ce_car1h <- getdf(ce_car1h)
ce_hch <- getdf(ce_hch)
ce_gl <- getdf(ce_gl)

co_aph <- getdf(co_aph)
co_csh <- getdf(co_csh)
co_gl <- getdf(co_gl)

il_aph <- getdf(il_aph)
il_omh <- getdf(il_omh)
il_rgh <- getdf(il_rgh)

ce_car1h_diff <- ce_car1h %>% arrange(p.adjust) %>% do(head(., n = 10))
ce_hch_diff <- ce_hch %>% arrange(p.adjust) %>% do(head(., n = 10))
ce_gl_diff <- ce_gl  %>% arrange(p.adjust) %>% do(head(., n = 10))

co_aph_diff <- co_aph %>% arrange(p.adjust) %>% do(head(., n = 10))
co_csh_diff <- co_csh %>% arrange(p.adjust) %>% do(head(., n = 10))
co_gl_diff <- co_gl %>% arrange(p.adjust) %>% do(head(., n = 10))

il_aph_diff <- il_aph %>% arrange(p.adjust) %>% do(head(., n = 10))
il_omh_diff <- il_omh %>% arrange(p.adjust) %>% do(head(., n = 10))
il_rgh_diff <- il_rgh %>% arrange(p.adjust) %>% do(head(., n = 10))

df_diff_top10 <- bind_rows(ce_car1h_diff, ce_hch_diff, ce_gl_diff, co_aph_diff, co_csh_diff, co_gl_diff, il_aph_diff, il_omh_diff, il_rgh_diff)

diff_go <- unique(df_diff_top10$ID)
targetct <- c(celist, colist, illist)
df_all_diffgo <- df_all[which(df_all$ID %in% c(diff_go)),] %>% filter(celltypes %in% targetct)

df_diff <- df_all_diffgo %>% mutate(CelltypeGroup = recode(celltypes,
                                            "EC (Cmss1 high)" = "colon_EP_EC (Cmss1 high)",
                                            "EC (Saa1 high)" = "colon_EP_EC (Saa1 high)",
                                            "Goblet2" = "colon_EP_Goblet2",
                                            "EC (Car1 high)" = "cecum_EP_EC (Car1 high)",
                                            "EC (Hmgcs2 high)" = "cecum_EP_EC (Hmgcs2 high)",
                                            "Goblet1" = "cecum_EP_Goblet1",
                                            "EC (Apoa4 high)" = "ileum_EP_EC (Apoa4 high)",
                                            "EC (Olfm4 high)" = "ileum_EP_EC (Olfm4 high)",
                                            "EC (Reg3g high)" = "ileum_EP_EC (Reg3g high)"
                                            ))




df_diff$p.adjust <- round(df_diff$p.adjust,3)
# df_diff_p <- reshape2::dcast(df_diff[,c("ID", "CelltypeGroup", "GeneRatio")], ID~CelltypeGroup, value.var = "GeneRatio")
write.csv(df_diff,"/Users/dengysh/Work/Projects/SpatialTranscriptomics/GermFreeMiceProject/GetFigures/GOPlots/Plots_20230830/df_diff.csv")

df_diff <- df_diff %>% mutate(idct = paste0(ID, CelltypeGroup))

df_rep_index <- df_diff %>% filter(duplicated(idct)) %>% pull(idct)
df_rep <- df_diff %>% filter(idct %in% df_rep_index)
df_uniq <- df_diff %>% filter(!idct %in% df_rep_index)

df_rep_corrt <- NULL
for (i in df_rep_index){
  df_tmp <- df_rep[df_rep$idct==i,] %>% filter(qvalue == max(qvalue))
  df_rep_corrt <- rbind(df_rep_corrt, df_tmp)
}
df_rep_corrt$group <- 2

df_diff_plot_corrected <- rbind(df_rep_corrt, df_uniq)

# df_diff_plot_corrected <- df_diff_plot_corrected %>% mutate(padj_sign = case_when((group == 1 & p.adjust < 0.001) ~ 1,
#                                                                                   (group == 1 & p.adjust < 0.01) ~ 2,
#                                                                                   (group == 1 & p.adjust < 0.1) ~ 3,
#                                                                                   (group == -1 & p.adjust < 0.001) ~ 4,
#                                                                                   (group == -1 & p.adjust < 0.01) ~ 5,
#                                                                                   (group == -1 & p.adjust < 0.1) ~ 6,
#                                                                                   (group == 2 & p.adjust < 0.001) ~ 7,
#                                                                                   (group == 2 & p.adjust < 0.01) ~ 8,
#                                                                                   (group == 2 & p.adjust < 0.1) ~ 9))

df_diff_plot_corrected <- df_diff_plot_corrected %>% mutate(padj_sign = case_when( p.adjust < 0.001 ~ 1,
                                                                                   p.adjust < 0.01 ~ 2,
                                                                                   p.adjust < 0.1 ~ 3))


df_diff_heatmap <- reshape2::dcast(df_diff_plot_corrected[,c("Description", "CelltypeGroup", "padj_sign")], Description~CelltypeGroup, value.var = "padj_sign")
df_diff_heatmap <- df_diff_heatmap %>% tibble::column_to_rownames(var="Description")
df_diff_heatmap[is.na(df_diff_heatmap)] <- 0

P <- reshape2::dcast(df_diff_plot_corrected[,c("Description", "CelltypeGroup", "group")], Description~CelltypeGroup, value.var = "group")
P <- P %>% tibble::column_to_rownames(var="Description")
P[is.na(P)] <- 0

d <- dist(P, method = "euclidean")
fit2 <- hclust(d, method="ward.D")
PB <- P[fit2$order,]

PB <- matrix(ifelse(PB == -1, "↓", ifelse(PB == 1, "↑",ifelse(PB == 2,"→"," "))), nrow(PB))
df_diff_heatmap <- df_diff_heatmap[fit2$order,]

rownames(PB) <- rownames(df_diff_heatmap)
colnames(PB) <- colnames(df_diff_heatmap)
write.csv(PB,"/Users/dengysh/Work/Projects/SpatialTranscriptomics/GermFreeMiceProject/GetFigures/GOPlots/Plots_20230830/df_diff_heatmap.csv")

library(ComplexHeatmap)
library(circlize)
col_fun1 = colorRamp2(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), c("lightgray", "#D62728FF","#FF7F0EFF", "#FFBB78FF", "blue", "lightblue", "#c9c1ff", "green", "#09d15a", "#86dda4"))
# cell_fun <- function(j, i, x, y, width, height, fill) {
#   grid.text(PB[i, j], x, y,gp = gpar(fontsize = 15))
# }
arrowcolor <- c("black", "white", "steelblue")
# arrowcolor <- c("red", "blue", "green")
cell_fun <- function(j, i, x, y, width, height, fill) {
  if (PB[i, j] == '↑'){
    grid.text(PB[i, j], x, y,gp = gpar(col = arrowcolor[1], fontsize = 15))
  }else if(PB[i, j] == '↓'){
    grid.text(PB[i, j], x, y,gp = gpar(col = arrowcolor[2], fontsize = 15))
  }else if(PB[i, j] == '→'){
    grid.text(PB[i, j], x, y,gp = gpar(col = arrowcolor[3], fontsize = 15))
  }else{
    grid.text(PB[i, j], x, y,gp = gpar(fontsize = 15))
  }
  
}

library(showtext)
pdf("/Users/dengysh/Work/Projects/SpatialTranscriptomics/GermFreeMiceProject/GetFigures/GOPlots/Plots_20230830/pheatmap.pdf",width = 40,height = 25)
showtext_begin() 
Heatmap(df_diff_heatmap, 
        col = col_fun1,
        cell_fun = cell_fun,
        rect_gp = gpar(col = "gray", lwd = 0.5),
        column_labels = c("EC (Car1 high)","EC (Hmgcs2 high)","Goblet1","EC (Cmss1 high)","EC (Saa1 high)","Goblet2","EC (Apoa4 high)","EC (Olfm4 high)","EC (Reg3g high)"),
        show_row_names = T,
        row_names_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize = 20),
        row_title = NULL,
        column_title_gp = gpar(fontsize = 20),
        cluster_columns = F,
        cluster_rows = F,
        column_split=rep(c("Cecum","Colon","Ileum"), each = 3),
        border="gray",
        show_row_dend=F,
        show_heatmap_legend = F,
        width = ncol(df_diff_heatmap)*unit(8, "mm"),
        height = nrow(df_diff_heatmap)*unit(8, "mm"))

at1 <- c(0, 1, 2, 3)
lgd1 = Legend(at = at1, labels =c("None","p<0.001","p<0.01","p<0.1"), border="gray",labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 25, fontface = "bold"),
              title = "", legend_gp = gpar(fill = col_fun1(at1)), grid_height = unit(10, "mm"), grid_width = unit(10, "mm"))

# at2 <- c(0, 4, 5, 6)
# lgd2 = Legend(at = at2, labels =c("None","p<0.001","p<0.01","p<0.1"), border="gray",labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 25, fontface = "bold"),
#               title = "Down", legend_gp = gpar(fill = col_fun1(at2)), grid_height = unit(10, "mm"), grid_width = unit(10, "mm"))
# 
# at3 <- c(0, 7, 8, 9)
# lgd3 = Legend(at = at3, labels =c("None","p<0.001","p<0.01","p<0.1"), border="gray",labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 25, fontface = "bold"),
#               title = "Both", legend_gp = gpar(fill = col_fun1(at3)), grid_height = unit(10, "mm"), grid_width = unit(10, "mm"))
# 
at4 <- c("↑","↓","→")
lgd4 = Legend(at = at4, labels = c("All upregulated","All downregulated","50% up, 50% down"), labels_gp = gpar(fontsize = 20),
              type = "points", legend_gp = gpar(col = arrowcolor), title = "", pch = at4, size = 20, background = "lightgray", grid_height = unit(10, "mm"), grid_width = unit(10, "mm"))

# lgd_list_vertical <- packLegend(lgd1, lgd2, lgd3, lgd4, direction = "horizontal")
lgd_list_vertical <- packLegend(lgd1, lgd4, direction = "horizontal")
draw(
  lgd_list_vertical,
  x = unit(0.60, "npc"), 
  y = unit(0.1, "npc"),
)

showtext_end()
dev.off() 


# library(showtext)
# library(ComplexHeatmap)
# library(circlize)
# pdf("/Users/dengysh/Desktop/pheatmap.pdf",width = 4,height = 15)
# showtext_begin() 
# pheatmap(df_diff_heatmap, display_numbers = matrix(ifelse(P == -1, "↓", ifelse(P == 1, "↑",ifelse(P == 2,"-"," "))), nrow(P)),
#          cluster_cols = F, cluster_rows = F, gaps_col = c(3, 6), color = colorRampPalette(colors = c("white","green","yellow","red"))(4),
#          fontsize_row = 5,fontsize_col=5,cellheight = 7,cellwidth = 7)
# 
# # col_fun1 = colorRamp2(c(-1, 0, 1, 2), c("white", "#E64B35FF","#4DBBD5FF", "#FFCD00FF"))
# at2 <- rev(c(-1, 0, 1, 2))
# lgd1 = Legend(at = at2, labels =c("Both","GF","SPF","None"), border="gray",
#               title = "EnrichedGroup", legend_gp = gpar(fill = col_fun1(at2)))
# showtext_end()
# dev.off()




