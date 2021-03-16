library("DiffBind")
#cd /Library/Frameworks/R.framework/Versions/3.5/Resources/library/DiffBind
samples <- read.csv(file.path(system.file("extra", package="DiffBind"), "H3K27me3_0619.csv"))
names(samples)
samples
H3K27me3 <- dba(sampleSheet = "H3K27me3_0619.csv", dir=system.file("extra", package="DiffBind"))

plot(H3K27me3)
H3K27me3 <- dba.count(H3K27me3, summits = 250)
H3K27me3
plot(H3K27me3)

H3K27me3$contrasts=NULL
H3K27me3 <- dba.contrast(H3K27me3, H3K27me3$masks$Unicellular, H3K27me3$masks$Multicellular, "Uni", "Multi")
H3K27me3
H3K27me3 <- dba.analyze(H3K27me3)
plot(H3K27me3, contrast=1)
H3K27me3.DB <- dba.report(H3K27me3, DataType=DBA_DATA_FRAME)
H3K27me3.DB

dba.plotPCA(H3K27me3, attributes=DBA_FACTOR, vColors = c("red", "skyblue1", "lightskyblue3", "turquoise3"), components=2:3)

dba.plotMA(H3K27me3)
