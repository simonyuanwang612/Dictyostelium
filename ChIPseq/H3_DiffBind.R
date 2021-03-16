library("DiffBind")
#cd /Library/Frameworks/R.framework/Versions/3.5/Resources/library/DiffBind
samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "H3_062418.csv"))
names(samples)
samples
H3 <- dba(sampleSheet = "H3_062418.csv", dir=system.file("extra", package="DiffBind"))

plot(H3)
H3 <- dba.count(H3, summits = 250)
H3
plot(H3)


H3$contrasts=NULL
H3 <- dba.contrast(H3, H3$masks$Unicellular, H3$masks$Multicellular, "Uni", "Multi")
H3
H3 <- dba.analyze(H3)
plot(H3)
H3.DB <- dba.report(H3, DataType=DBA_DATA_FRAME)
H3.DB
write.table(H3.DB,"H3_0624.txt",sep="\t",row.names=FALSE)

dba.plotPCA(H3, attributes=DBA_FACTOR, vColors = c("red", "skyblue1", "lightskyblue3", "turquoise3"), components=1:2)

dba.plotMA(H3)
dba.plotVolcano(H3)
sum(H3.DB$Fold<0)
sum(H3.DB$Fold>0)
pvals <- dba.plotBox(H3)
pvals
# correlation heatmaps
corvals <- dba.plotHeatmap(H3, contrast=1)
# Binding affinity heatmaps
corvals <- dba.plotHeatmap(H3, contrast=1, correlation = FALSE, scale="row")
