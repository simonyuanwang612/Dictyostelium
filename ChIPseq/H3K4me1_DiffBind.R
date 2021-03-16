library("DiffBind")
#cd /Library/Frameworks/R.framework/Versions/3.5/Resources/library/DiffBind
samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "H3K4me1_0619.csv"))
names(samples)
samples
H3K4me1 <- dba(sampleSheet = "H3K4me1_0619.csv", dir=system.file("extra", package="DiffBind"))

plot(H3K4me1)
H3K4me1 <- dba.count(H3K4me1, summits = 250)
H3K4me1
plot(H3K4me3)


H3K4me1$contrasts=NULL
H3K4me1 <- dba.contrast(H3K4me1, H3K4me1$masks$Unicellular, H3K4me1$masks$Multicellular, "Uni", "Multi")
H3K4me1
H3K4me1 <- dba.analyze(H3K4me1)
plot(H3K4me1, contrast=1)
H3K4me1.DB <- dba.report(H3K4me1, DataType=DBA_DATA_FRAME)
H3K4me1.DB
write.table(H3K4me1.DB,"H3K4me1_0619.txt",sep="\t",row.names=FALSE)

dba.plotPCA(H3K4me1, contrast=1, attributes=DBA_FACTOR, vColors = c("red", "skyblue1", "lightskyblue3", "turquoise3"), components=1:2)

dba.plotMA(H3K4me1)
dba.plotVolcano(H3K4me1)
sum(H3K4me1.DB$Fold<0)
sum(H3K4me1.DB$Fold>0)
pvals <- dba.plotBox(H3K4me1)
pvals
# correlation heatmaps
corvals <- dba.plotHeatmap(H3K4me1, contrast=1)
# Binding affinity heatmaps
corvals <- dba.plotHeatmap(H3K4me1, contrast=1, correlation = FALSE, scale="row")
