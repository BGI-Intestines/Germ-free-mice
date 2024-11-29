library(argparse, quietly = TRUE)
parser = argparse::ArgumentParser(description = 'Script for make pathway list with tissue')
parser$add_argument('-i', '--input', dest = 'input', help = 'input file of list of pathway')
parser$add_argument('-o', '--out', dest = 'out', help = 'output file of list of pathway with a list of library names')
parser$add_argument('-s', '--sample', dest = 'sample', help = 'Input smapeloutput file of list of pathway')
parser$add_argument('-a', '--assay', dest = 'assay' , help = "Which assay do you choose!") # 2023 - 1 - 12
opts = parser$parse_args()
print(opts)

library(SeuratDisk)
library(Seurat)

# Read file #
RDS <- readRDS(opts$input)
file <- paste(opts$out,"/",opts$sample,".h5Seurat",sep = "")
print(file)
SaveH5Seurat(RDS, filename = file)
Convert(file, dest = 'h5ad', assay = opts$assay, overwrite = T) # 2023  - 1 -12
print("to H5ad Done!")

