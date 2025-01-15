###############################################################
### functions for Layer score and plot #######################
### Layers 
layerScore <- function(name,object,n_bin=15){
  # gene list to calculate CV/PV, CYP2E1 为CV特异，CDH1（E-CAD)为PV特异
  CV=c("Glul","Cyp2e1","Cyp1a2")
  PV=c("Alb","Ass1","Asl","Cyp2f2") #"ALB",
  #CV=c("GLUL","CYP2E1")
  #PV=c("CDH1")
  # AddModulescore, CV-PV 
  object <-AddModuleScore(object,features=CV,name="CV",nbin=n_bin)  # default 24. came out errors then set to 15
  object <-AddModuleScore(object,features=PV,name="PV",nbin=n_bin)
  
  object$CVmean <- rowMeans(object@meta.data[,c("CV1","CV2","CV3")]) # # numbers is CV gene types
  object$PVmean <- rowMeans(object@meta.data[,c("PV1","PV2","PV3","PV4")])
  object$CV_PV <- object$CVmean-object$PVmean
  
  # seperate into 9 layers
  #Calculate breaks to divide into 9 layers
  breaks <- quantile(object$CV_PV, probs = seq(0, 1, by = 1/9), na.rm = TRUE)
  # Create Group column based on quantiles
  object$Layer <- cut(object$CV_PV, breaks = breaks, labels = FALSE, include.lowest = TRUE) # cv和pv顺序相反
  object$Layer <- abs(10-object$Layer)
  object$Layer <- paste0("L", object$Layer)
  layer_order <- c("L1","L2","L3","L4","L5","L6","L7","L8","L9")
  object@meta.data$Layer <- factor(object@meta.data$Layer, levels=layer_order)
  # Save the result
  nameSave <- paste0("meta.data.",name,".layers.tsv")
  write.table(object@meta.data, file=nameSave, col.names=NA, sep="\t", quote=F)
  
  # return object
  return(object)
}

layerPlot <- function(name,object,pt.size=NA){
  # Scatter plot 
  pLa <- ggplot(object@meta.data, aes(x = x, y = y, color = Layer)) +
    geom_point_rast(size = pt.size,shape = 16) +  # Customize point size if needed, default 0.1
    scale_color_manual(values = c("gray10","gray20","gray30","gray40", "gray50", "gray60", "gray70", 
                                  "gray80", "gray90")) +
    
    coord_equal() + theme_void() + 
    geom_raster(alpha=0.1) + # 未起作用
    guides(color = guide_legend(override.aes = list(shape = 15, size = 8))) 
  
  ggsave(paste0("Layers.9.", name, ".pdf"), pLa, width=6, height=6)
  ggsave(paste0("Layers.9.", name, ".png"), pLa, width=6, height=6, units="in", dpi=300)
  #theme_minimal()  
  return(pLa)
}
# # 
gf_liv<-readRDS(file="GF_liver_bin20.rds")
spf_liv<-readRDS(file="SPF_liver_bin20.rds")

GF <- layerScore(name="GF",gf_liv,n_bin=5)
layerPlot(name="GF",GF,pt.size=0.01)