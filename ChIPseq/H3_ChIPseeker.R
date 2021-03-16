library(ChIPseeker)
library(ggplot2)
library(biomaRt)
library("GenomicFeatures")
library(clusterProfiler)
files <- getSampleFiles()
print(files)
list.files(path = "/Library/Frameworks/R.framework/Versions/3.5/Resources/library/ChIPseeker/extdata/GEO_sample_data")

peak=GenomicRanges::GenomicRangesList(FB=readPeakFile(files[[2]]),
                                      FA=readPeakFile(files[[1]]),
                                      MB=readPeakFile(files[[9]]),
                                      MA=readPeakFile(files[[8]]),
                                      SB=readPeakFile(files[[11]]),
                                      SA=readPeakFile(files[[10]]),
                                      VB=readPeakFile(files[[13]]),
                                      VA=readPeakFile(files[[12]]))

peak_FA <- readPeakFile(files[[1]])
peak_FA
p <- covplot(peak, weightCol=5, title = "H3 Peaks over Chromosomes", chrs=c("DDB0232428", "DDB0232429", "DDB0232430", "DDB0232431", "DDB0232432", "DDB0232433" ))


print(p)
col <- c(VA='dodgerblue1', MA='deepskyblue1', SA='skyblue1', FA='turquoise1',
         VB='dodgerblue3', MB='deepskyblue3', SB='skyblue3', FB='turquoise3')
p + facet_grid(chr ~ .id) + scale_color_manual(values=col) + scale_fill_manual(values=col)
dicty_txdb_ensemble <-  makeTxDbFromBiomart(biomart="protists_mart",
                                            dataset="ddiscoideum_eg_gene",
                                            
                                            host="protists.ensembl.org")
dicty_txdb_ensemble1 <-  makeTxDbFromBiomart(biomart="protists_mart",
                                            dataset="ddiscoideum_eg_gene",
                                            
                                            host="protists.ensembl.org")

dropSeqlevels(dicty_txdb_ensemble, "CH709182", pruning.mode="coarse")
seqlevels(dicty_txdb_ensemble)
print(files)
promoter1 <- getPromoters(TxDb=dicty_txdb_ensemble1, upstream=1500, downstream=1500)
promoter <- renameSeqlevels(promoter, c("6" = "DDB0232433"))
promoter <- renameSeqlevels(promoter, c("5" = "DDB0232432"))
promoter <- dropSeqlevels(promoter, "DDB0237465", pruning.mode="coarse")
promoter

tagMatrixList_new <- lapply(peak, getTagMatrix, windows=promoter)

tagMatrixList1
getTagMatrix
plotAvgProf(tagMatrixList1, xlim=c(-1500, 1500))
tagHeatmap(tagMatrixList1, xlim=c(-1500, 1500), color = "blue")

peakAnnoList <- lapply(files[[1]], annotatePeak, TxDb=dicty_txdb_ensemble,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
renameSeqlevels(dicty_txdb_ensemble, c("6"="DDB0232433"))

seqlevels(dicty_txdb_ensemble1)
promoter



files[1]
#there is a bug here, have to do the following two step inorder to annotatePeak
peakAnno <- annotatePeak(files[[1]], tssRegion=c(-1500, 1500),
                         TxDb=dicty_txdb_ensemble)



peakAnnoList <- lapply(peak, annotatePeak, TxDb=dicty_txdb_ensemble,
                       tssRegion=c(-500, 500), verbose=FALSE)


plotAnnoBar(peakAnnoList)
