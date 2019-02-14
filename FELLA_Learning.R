#Title: learning how to use fella packages
#Author: xiangruikong
#E-mail:xiangruikong@bioxxx.cn
#Time: 2019-02-13
##########################################
#1.Building the database
setwd("D:\\learngit\\bioinformatics-R")
library(FELLA)
set.seed(1)
graph <- buildGraphFromKEGGREST(organism = "ath")
rm(list=ls())
memory.limit(size = NA)
memory.size()
memory.size(max=T)
gc()
memory.limit(15000)
memory.size(max=F)
tmpdir <- paste0(tempdir(), "/my_database")
unlink(tmpdir, recursive = TRUE)
buildDataFromGraph(keggdata.graph = graph,
  databaseDir = tmpdir,
  internalDir = FALSE,
  matrices = "diffusion",
  normality = "diffusion",
  niter = 50)
fella.data <- loadKEGGdata(
  databaseDir = tmpdir,
  internalDir = FALSE,
  loadMatrix = "diffusion"
  )
fella.data
cat(getInfo(fella.data))
#2. Mapping the input metabolites
compounds.epithelial <- c(
  "C02862", "C00487", "C00025", "C00064",
  "C00670", "C00073", "C00588", "C00082", "C00043")
analysis.epithelial <- defineCompounds(
  compounds = compounds.epithelial,
  data = fella.data)
getInput(analysis.epithelial)
getExcluded(analysis.epithelial)
analysis.epithelial
#3.Enriching using diffusion
analysis.epithelial <- runDiffusion(
  object = analysis.epithelial,
  data = fella.data,
  approx = "normality")
analysis.epithelial
#4.plot
pdf("mygraph.pdf")
nlimit <- 220
vertex.label.cex <- .5
plot(
  analysis.epithelial,
  method = "diffusion",
  data = fella.data,
  nlimit = nlimit,
  vertex.label.cex = vertex.label.cex) 
  dev.off()
  getwd()
#5.Exporting the results


