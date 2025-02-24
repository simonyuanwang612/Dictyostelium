files <- getSampleFiles()
list.files(path = "/Library/Frameworks/R.framework/Versions/3.5/Resources/library/ChIPseeker/extdata/GEO_sample_data")
peak=GenomicRanges::GenomicRangesList(FB=readPeakFile(files[[10]]),
                                      FA=readPeakFile(files[[3]]),
                                      MB=readPeakFile(files[[25]]),
                                      MA=readPeakFile(files[[18]]),
                                      SB=readPeakFile(files[[40]]),
                                      SA=readPeakFile(files[[33]]),
                                      VB=readPeakFile(files[[56]]),
                                      VA=readPeakFile(files[[48]]))
p <- covplot(peak, weightCol=5, title = "H3K27me3 Peaks over Chromosomes", chrs=c("DDB0232428", "DDB0232429", "DDB0232430", "DDB0232431", "DDB0232432", "DDB0232433" ))


print(p)
col <- c(VA='dodgerblue1', MA='deepskyblue1', SA='skyblue1', FA='turquoise1',
         VB='dodgerblue3', MB='deepskyblue3', SB='skyblue3', FB='turquoise3')
p + facet_grid(chr ~ .id) + scale_color_manual(values=col) + scale_fill_manual(values=col)
promoter <- getPromoters(TxDb=dicty_txdb_ensemble, upstream=1500, downstream=1500)
tagMatrixList <- lapply(peak, getTagMatrix, windows=promoter)
plotAvgProf(tagMatrixList, xlim=c(-1500, 1500))
peakAnnoList <- lapply(peak, annotatePeak, TxDb=dicty_txdb_ensemble,
                       tssRegion=c(-500, 500), verbose=FALSE)
plotAnnoBar(peakAnnoList)
