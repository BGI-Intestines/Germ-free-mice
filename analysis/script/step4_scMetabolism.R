.libPaths(c("/software/anaconda3/envs/R411/lib/R/library", .libPaths()))
# author : lwm 
# data : 2023-2-17

### Get the parameters
parser = argparse::ArgumentParser(description="Script to subcelltypes Clustering of scRNA data")
parser$add_argument('-I','--input', help='input a RDS')
parser$add_argument('-M','--method',help="which method do you choose?,VISION,AUCell,ssGSEA,GSVA")
parser$add_argument('-O','--output', help='out directory')
opts = parser$parse_args()
print(opts)

library(homologene) # 人鼠同源转换 mouse2human(genelist2)
library(scMetabolism)
library(ggplot2)
library(rsvd)


#' scMetabolism
#' scMetabolism
#' @param obj
#' @keywords scMetabolism
#' @examples
#' sc.metabolism.Seurat()
#' @export sc.metabolism.Seurat
sc.metabolism.Seurat <- function(obj, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG",threshold=0.01,projection_genes="threshold") {
  countexp<-obj@assays$RNA@counts
  countexp<-data.frame(as.matrix(countexp))
  signatures_KEGG_metab <- system.file("data", "mouse_KEGG_metabolism_nc.gmt", package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "mouse_REACTOME_metabolism.gmt", package = "scMetabolism")
  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
  #imputation
  if (imputation == F) {
    countexp2<-countexp
  }
  if (imputation == T) {
    cat("Start imputation...\n")
    #Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588
    #Github: https://github.com/KlugerLab/ALRA
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(countexp))
    countexp2 <- result.completed[[3]]; row.names(countexp2) <- row.names(countexp)
  }
  #signature method
  cat("Start quantify the metabolism activity...\n")
  #VISION
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2) / n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile,threshold = threshold,projection_genes=c(projection_genes))
    options(mc.cores = 10)
    vis <- analyze(vis)
    signature_exp<-data.frame(t(vis@SigScores))
  }
  #AUCell
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F) #rank
    geneSets <- getGmt(gmtFile) #signature read
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) #calc
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  #ssGSEA
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  #GSVA
  if (method == "GSVA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }
  #obj@assays$METABOLISM$score<-signature_exp
  signature_exp
}

####################################################################   run scMetabolism  #######################################################
obj=readRDS(opts$input)

kegg=sc.metabolism.Seurat(obj, method = opts$method, imputation = F, ncores = 2, metabolism.type = "KEGG",threshold=0.01,projection_genes="threshold")
reactome=sc.metabolism.Seurat(obj, method = opts$method, imputation = F, ncores = 2, metabolism.type = "REACTOME",threshold=0.01,projection_genes="threshold")


outdir=opts$output
write.table(kegg,file=paste(outdir,"/kegg.xls",sep=""),row.names=rownames(kegg),col.names=colnames(kegg),sep="\t",quote=F)
write.table(reactome,file=paste(outdir,"/reactome.xls",sep=""),row.names=rownames(reactome),col.names=colnames(reactome),sep="\t",quote=F)

