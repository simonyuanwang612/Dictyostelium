graphics.off()
library("DiffBind")
#cd /Library/Frameworks/R.framework/Versions/3.5/Resources/library/DiffBind
samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "H3K4me3_0619new.csv"))
names(samples)
samples
H3K4me3 <- dba(sampleSheet = "H3K4me3_0619new.csv", dir=system.file("extra", package="DiffBind"))

plot(H3K4me3)
H3K4me3 <- dba.count(H3K4me3, summits = 250)
H3K4me3
plot(H3K4me3)

H3K4me3

H3K4me3$contrasts=NULL
H3K4me3 <- dba.contrast(H3K4me3, H3K4me3$masks$Unicellular, H3K4me3$masks$Multicellular, "Uni", "Multi")
H3K4me3
H3K4me3 <- dba.analyze(H3K4me3)
plot(H3K4me3, contrast=1)
H3K4me3.DB <- dba.report(H3K4me3, DataType=DBA_DATA_FRAME)
H3K4me3.DB
write.table(H3K4me3.DB,"H3K4me3_0619.txt",sep="\t",row.names=FALSE)

dba.plotPCA(H3K4me3, contrast=1, attributes=DBA_FACTOR, vColors = c("red", "skyblue1", "lightskyblue3", "turquoise3"), components=1:2)
                                                                    

dba.plotMA(H3K4me3)
dba.plotVolcano(H3K4me3)
sum(H3K4me3.DB$Fold<0)
sum(H3K4me3.DB$Fold>0)
pvals <- dba.plotBox(H3K4me3)
pvals
# correlation heatmaps
corvals <- dba.plotHeatmap(H3K4me3, contrast=1)
# Binding affinity heatmaps
corvals <- dba.plotHeatmap(H3K4me3, contrast=1, correlation = FALSE, scale="row")
# Type the following two commands together, accounting for batch effects(replicates)
H3K4me3$contrasts=NULL
H3K4me3 <- dba.contrast(H3K4me3, H3K4me3$masks$Unicellular, H3K4me3$masks$Multicellular, "Uni", "Multi", block = DBA_REPLICATE)

H3K4me3 <- dba.analyze(H3K4me3, method=DBA_ALL_METHODS)
H3K4me3
dba.plotMA(H3K4me3, method = DBA_DESEQ2_BLOCK)
H3K4me3.DB_block <- dba.report(H3K4me3, method = DBA_DESEQ2_BLOCK)
H3K4me3.DB_block
corvals_block <- dba.plotHeatmap(H3K4me3, contrast=1, correlation = FALSE, scale="row", method = DBA_DESEQ2_BLOCK)
dba.plotPCA(H3K4me3, contrast = 1, method = DBA_DESEQ2_BLOCK)
#compare four methods
H3K4me3 <- dba.analyze(H3K4me3, method=DBA_ALL_METHODS)
dba.show(H3K4me3,bContrasts=T)[9:12]
tam.block<- dba.report(H3K4me3,  method = DBA_ALL_METHODS_BLOCK, bDB = TRUE, bAll = TRUE)
tam.block                     
dba.plotVenn(tam.block, 1:4, label1="edgeR",label2="DESeq2", label3="edgeR Blocked", label4="DESeq2 Blocked")            
#dba.overlap(H3K4me3, H3K4me3$masks$Unicellular, mode = DBA_OLAP_RATE, contrast = 1)

#olap.rate <- dba.overlap(H3K4me3, mode = DBA_OLAP_RATE)
#H3K4me3 <- dba.peakset(H3K4me3, consensus=DBA_CONDITION, minOverlap = 0.33)

# This step is to create a rda file 
H3K4me3 <- dba(sampleSheet = "H3K4me3_0601.csv", dir=system.file("extra", package="DiffBind"))
save(H3K4me3, file = "H3K4me3.rda")
data(H3K4me3)
olap.rate <- dba.overlap(H3K4me3, mode = DBA_OLAP_RATE)
olap.rate
plot(olap.rate,type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')
names(H3K4me3$masks)
dba.overlap(H3K4me3, H3K4me3$masks$Unicellular, mode=DBA_OLAP_RATE)
dba.plotVenn(H3K4me3, H3K4me3$masks$Unicellular)
dba.overlap(H3K4me3, H3K4me3$masks$Multicellular, mode=DBA_OLAP_RATE)
dba.plotVenn(H3K4me3, H3K4me3$masks$Multicellular)
H3K4me3_consensus <- dba.peakset(H3K4me3, consensus = c(DBA_CONDITION), minOverlap = 0.66)
H3K4me3_consensus
H3K4me3_consensus <- dba(H3K4me3_consensus, mask= H3K4me3_consensus$masks$Consensus, minOverlap = 1)
H3K4me3_consensus
consensus_peaks <- dba.peakset(H3K4me3_consensus, bRetrieve=TRUE)
H3K4me3 <- dba.count(H3K4me3, peaks=consensus_peaks)
H3K4me3
data(H3K4me3)
H3K4me3 <- dba.peakset(H3K4me3, consensus = DBA_CONDITION, minOverlap = 0.66)
dba.plotVenn(H3K4me3, H3K4me3$masks$Consensus)
#6.3
data(H3K4me3)
dba.overlap(H3K4me3, H3K4me3$masks$Unicellular, mode=DBA_OLAP_RATE)
dba.overlap(H3K4me3, H3K4me3$masks$Multicellular, mode=DBA_OLAP_RATE)
H3K4me3 <- dba.peakset(H3K4me3, consensus = DBA_CONDITION, minOverlap = 0.33)
dba.plotVenn(H3K4me3, H3K4me3$masks$Consensus)
H3K4me3.OL <- dba.overlap(H3K4me3, H3K4me3$masks$Consensus)
H3K4me3.OL$onlyA
H3K4me3.OL$onlyB
# Comparison of occupancy and affinity based analysis
H3K4me3 <- dba.peakset(H3K4me3, H3K4me3$masks$Consensus, minOverlap = 1, sampID = "OL Consensus")
H3K4me3 <- dba.peakset(H3K4me3, !H3K4me3$masks$Consensus, minOverlap = 3, sampID = "Consensus_3")
H3K4me3
dba.plotVenn(H3K4me3, 9:10)
