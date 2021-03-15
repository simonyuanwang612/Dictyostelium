library(ChIPseeker)
library(clusterProfiler)
library("ChIPpeakAnno")
library("GenomicFeatures")

dicty_txdb_ensemble <-  makeTxDbFromBiomart(biomart="protists_mart",
                                            dataset="ddiscoideum_eg_gene",
                                             host="protists.ensembl.org")

dicty_txdb_ensemble

# It's better to rename and use the dropseqlevels here for the consistency and conveniency later
renameSeqlevels(dicty_txdb_ensemble, c("6"="DDB0232433"))
seqlevels(dicty_txdb_ensemble)
dropSeqlevels(dicty_txdb_ensemble, "CH709182", pruning.mode="coarse")

files <- getSampleFiles()
files
write.table(files, "file_number9.txt", sep="\t")

# Read the merged peak files
peak_merge_compare=GenomicRanges::GenomicRangesList(V=readPeakFile(files[[100]]),
                                                S=readPeakFile(files[[99]]),
                                                M=readPeakFile(files[[98]]),
                                                F=readPeakFile(files[[97]]))

# Change to readable chromosome ID, as dicty has an average of 1kb plus promoter, we define the upstream and downstream as 1500 bp
promoter <- getPromoters(TxDb=dicty_txdb_ensemble, upstream=1500, downstream=1500)
promoter <- dropSeqlevels(promoter, "DDB0169550", pruning.mode="coarse")
promoter <- dropSeqlevels(promoter, "DDB0215018", pruning.mode="coarse")
promoter <- dropSeqlevels(promoter, "DDB0215151", pruning.mode="coarse")
promoter <- dropSeqlevels(promoter, "DDB0220052", pruning.mode="coarse")
promoter <- dropSeqlevels(promoter, "DDB0237465", pruning.mode="coarse")

# Plot the coverage plot of peak regions over chromosome and generate a figure to visualize
covplot(peak_merge_compare, weightCol=5, chrs=c("DDB0232428", "DDB0232429", "DDB0232430", "DDB0232431", "DDB0232432", "DDB0232433" ))


# calculate the tag matrix 
tagMatrixList_merge_compare <- lapply(peak_merge_compare, getTagMatrix, windows=promoter)

# plot the profile of peaks
plotAvgProf(tagMatrixList_merge_compare, xlim=c(-1500, 1500))

# plot the heatmap of tagmatrix
tagHeatmap(tagMatrixList_merge_compare, xlim=c(-1500, 1500), color="orange")

# Annotate peaks
peakAnnoList <- lapply(peak_merge_compare, annotatePeak, TxDb=dicty_txdb_ensemble,
                       tssRegion=c(-500, 500), verbose=FALSE)

# Visualize genomic annotation in a bar plot
plotAnnoBar(peakAnnoList)

V_annot <- as.data.frame(peakAnnoList[[1]]@anno)
S_annot <- as.data.frame(peakAnnoList[[2]]@anno)
M_annot <- as.data.frame(peakAnnoList[[3]]@anno)
F_annot <- as.data.frame(peakAnnoList[[4]]@anno)


table(peakAnnoList[[4]]@detailGenomicAnnotation$genic)

# plotDistToTss to calculate the percentage of open chromatin sites upstream and downstream from the TSS of the nearest genes
plotDistToTSS(peakAnnoList, ylab = "Open sites (%) (5'->3')",
              title = "Distribution of open chromatin regions relative to TSS")


