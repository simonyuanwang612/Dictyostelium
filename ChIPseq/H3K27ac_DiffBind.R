library("DiffBind")
#cd /Library/Frameworks/R.framework/Versions/3.5/Resources/library/DiffBind
samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "H3K27ac_0619_new.csv"))
names(samples)
samples
H3K27ac <- dba(sampleSheet = "H3K27ac_0619_new.csv", dir=system.file("extra", package="DiffBind"))

plot(H3K27ac)
H3K27ac <- dba.count(H3K27ac, summits = 250)
H3K27ac
plot(H3K27ac)

H3K27ac$contrasts=NULL
H3K27ac <- dba.contrast(H3K27ac, H3K27ac$masks$Unicellular, H3K27ac$masks$Multicellular, "Uni", "Multi")
H3K27ac
H3K27ac <- dba.analyze(H3K27ac)
plot(H3K27ac, contrast=1)
H3K27ac.DB <- dba.report(H3K27ac, DataType=DBA_DATA_FRAME)
H3K27ac.DB
write.table(H3K27ac.DB,"H3K27ac_0619.txt",sep="\t",row.names=FALSE)

dba.plotPCA(H3K27ac, contrast=1, attributes=DBA_FACTOR, vColors = c("red", "skyblue1", "lightskyblue3", "turquoise3"), components=1:2)

dba.plotMA(H3K27ac)
dba.plotVolcano(H3K27ac)
sum(H3K27ac.DB$Fold<0)
sum(H3K27ac.DB$Fold>0)
pvals <- dba.plotBox(H3K27ac)
pvals
# correlation heatmaps
corvals <- dba.plotHeatmap(H3K27ac, contrast=1)
# Binding affinity heatmaps
corvals <- dba.plotHeatmap(H3K27ac, contrast=1, correlation = FALSE, scale="row")
