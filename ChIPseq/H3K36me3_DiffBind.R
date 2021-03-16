library("DiffBind")
#cd /Library/Frameworks/R.framework/Versions/3.5/Resources/library/DiffBind
samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "H3K36me3_0619.csv"))
names(samples)
samples
H3K36me3 <- dba(sampleSheet = "H3K36me3_0619.csv", dir=system.file("extra", package="DiffBind"))

plot(H3K36me3)
H3K36me3 <- dba.count(H3K36me3, summits = 250)
H3K36me3
plot(H3K36me3)

H3K36me3$contrasts=NULL
H3K36me3 <- dba.contrast(H3K36me3, H3K36me3$masks$Unicellular, H3K9me3$masks$Multicellular, "Uni", "Multi")
H3K36me3
H3K36me3 <- dba.analyze(H3K36me3)
plot(H3K36me3, contrast=1)
H3K36me3.DB <- dba.report(H3K36me3, DataType=DBA_DATA_FRAME)
H3K36me3.DB
write.table(H3K36me3.DB,"H3K36me3_0619.txt",sep="\t",row.names=FALSE)

dba.plotPCA(H3K36me3, contrast=1, attributes=DBA_FACTOR, vColors = c("red", "skyblue1", "lightskyblue3", "turquoise3"), components=1:2)

dba.plotMA(H3K36me3)
dba.plotVolcano(H3K36me3)
sum(H3K36me3.DB$Fold<0)
sum(H3K36me3.DB$Fold>0)
pvals <- dba.plotBox(H3K36me3)
# correlation heatmaps
corvals <- dba.plotHeatmap(H3K36me3, contrast=1)
# Binding affinity heatmaps
corvals <- dba.plotHeatmap(H3K36me3, contrast=1, correlation = FALSE, scale="row")
