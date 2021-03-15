library("DiffBind")
library("ChIPpeakAnno")

# For the pooled samples
library("DiffBind")

ATAC
ATAC <- dba(sampleSheet = "ATAC_0606.csv", dir=system.file("extra", package="DiffBind"))
plot(ATAC) 
ATAC1 <- dba.count(ATAC, summits = 200, minOverlap=0.1)

plot(ATAC1)
ATAC1
counts <- dba.peakset(ATAC1, bRetrieve=TRUE, writeFile="readscores_pool_0624.txt")
counts

ATAC2 <- dba.contrast(ATAC1, categories = DBA_CONDITION, minMembers = 3)
ATAC2
ATAC3 <- dba.analyze(ATAC2)


plot(ATAC3)  

dba.plotPCA(ATAC3,DBA_FACTOR,label=DBA_CONTROL)

dba.show(ATAC3, bContrast=T)
ATAC3
plot(ATAC3, contrast=1)

ATAC4 <- dba.report(ATAC3, th=0.05, bUsePval=FALSE,  bNormalized = TRUE, bCounts = TRUE)
ATAC4
write.table(ATAC4,"ATAC_DB_Pool_0624.txt", sep="\t",row.names=FALSE)
dba.plotPCA(ATAC3, contrast = 1,attributes=DBA_FACTOR, vColors = c("blue", "orange", "red", "grey"), label=DBA_CONTROL)

dba.plotMA(ATAC3)

crukBlue.save <- getFromNamespace("crukBlue", "DiffBind")
crukMagenta.save <- getFromNamespace("crukMagenta", "DiffBind")
assignInNamespace("crukBlue","black","DiffBind")
assignInNamespace("crukMagenta","orange","DiffBind")
dba.plotMA(ATAC3)

dba.plotVolcano(ATAC3)

sum(ATAC4$Fold<0)
sum(ATAC4$Fold>0)

pvals <- dba.plotBox(ATAC3)

# correlation heatmaps of pooled samples
corvals <- dba.plotHeatmap(ATAC3, contrast=1, ColAttributes = DBA_FACTOR, colSideCols = c("mediumpurple1", "palegreen2", "lightskyblue", "orange", "orchid1"), colScheme="Blues")

# Open chromatin regions heatmaps of pooled samples 
corvals <- dba.plotHeatmap(ATAC3, contrast=1, correlation = FALSE, scale="row", ColAttributes = DBA_FACTOR, colSideCols = c("mediumpurple1", "palegreen2", "lightskyblue", "orange", "orchid1"), colScheme="Blues")

#Start with Granges for ATAC3
counts <- dba.peakset(ATAC3, bRetrieve=TRUE)
counts
counts_gene <- annotatePeakInBatch(counts, AnnotationData=annoData)
counts_gene

#Start with Granges for ATAC.DB

# Give ranges numeric names in order
names(counts_gene) <- c(1:length(counts_gene))
counts_gene
write.table(counts_gene, "counts_gene_pool_06242020.txt", sep = "\t", col.names = T, row.names = F, quote = F)


merge_file <- merge(ATAC_DB_Pool_0624, counts_gene_pool_06242020, by =c("seqnames", "start", "end"))

write.table(merge_file, "merge_table_06242020_pool.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# For the fresh samples
library("DiffBind")

ATAC
ATAC <- dba(sampleSheet = "ATAC_0622.csv", dir=system.file("extra", package="DiffBind"))
plot(ATAC) 
ATAC1 <- dba.count(ATAC, summits = 200, minOverlap=0.1)

plot(ATAC1)
ATAC1
counts <- dba.peakset(ATAC1, bRetrieve=TRUE, writeFile="readscores_fresh_0624.txt")
counts

ATAC2 <- dba.contrast(ATAC1, categories = DBA_CONDITION, minMembers = 3)
ATAC2
ATAC3 <- dba.analyze(ATAC2)


plot(ATAC3)  

dba.plotPCA(ATAC3,DBA_FACTOR,label=DBA_CONTROL)

dba.show(ATAC3, bContrast=T)
ATAC3
plot(ATAC3, contrast=1)

ATAC4 <- dba.report(ATAC3, th=0.05, bUsePval=FALSE,  bNormalized = TRUE, bCounts = TRUE)
ATAC4
write.table(ATAC4,"ATAC_DB_Fresh_0624.txt", sep="\t",row.names=FALSE)
dba.plotPCA(ATAC3, contrast = 1,attributes=DBA_FACTOR, vColors = c("blue", "orange", "red", "grey"), label=DBA_CONTROL)

dba.plotMA(ATAC3)

crukBlue.save <- getFromNamespace("crukBlue", "DiffBind")
crukMagenta.save <- getFromNamespace("crukMagenta", "DiffBind")
assignInNamespace("crukBlue","black","DiffBind")
assignInNamespace("crukMagenta","orange","DiffBind")
dba.plotMA(ATAC3)

dba.plotVolcano(ATAC3)

sum(ATAC4$Fold<0)
sum(ATAC4$Fold>0)

pvals <- dba.plotBox(ATAC3)

# correlation heatmaps of fresh samples
corvals <- dba.plotHeatmap(ATAC3, contrast=1, ColAttributes = DBA_FACTOR, colSideCols = c("mediumpurple1", "palegreen2", "lightskyblue", "orange", "orchid1"), colScheme="Blues")

# Open chromatin regions heatmaps of fresh samples
corvals <- dba.plotHeatmap(ATAC3, contrast=1, correlation = FALSE, scale="row", ColAttributes = DBA_FACTOR, colSideCols = c("mediumpurple1", "palegreen2", "lightskyblue", "orange", "orchid1"), colScheme="Blues")

#Start with Granges for ATAC3 of fresh samples
counts <- dba.peakset(ATAC3, bRetrieve=TRUE)

counts
counts_gene <- annotatePeakInBatch(counts, AnnotationData=annoData)
counts_gene

#Start with Granges for ATAC.DB of fresh samples

# Give ranges numeric names in order
names(counts_gene) <- c(1:length(counts_gene))
counts_gene
write.table(counts_gene, "counts_gene_fresh_06242020.txt", sep = "\t", col.names = T, row.names = F, quote = F)


merge_file <- merge(ATAC_DB_Fresh_0624, counts_gene_fresh_06242020, by =c("seqnames", "start", "end"))

write.table(merge_file, "merge_table_06242020_fresh.txt", sep = "\t", col.names = T, row.names = F, quote = F)

# For the frozen samples

library("DiffBind")

ATAC
ATAC <- dba(sampleSheet = "ATAC_0624_frozen.csv", dir=system.file("extra", package="DiffBind"))
plot(ATAC) 
ATAC1 <- dba.count(ATAC, summits = 200, minOverlap=0.1)

plot(ATAC1)
ATAC1
counts <- dba.peakset(ATAC1, bRetrieve=TRUE, writeFile="readscores_frozen_0624.txt")
counts

ATAC2 <- dba.contrast(ATAC1, categories = DBA_CONDITION, minMembers = 2)
ATAC2
ATAC3 <- dba.analyze(ATAC2)


plot(ATAC3)  

dba.plotPCA(ATAC3,DBA_FACTOR,label=DBA_CONTROL)

dba.show(ATAC3, bContrast=T)
ATAC3
plot(ATAC3, contrast=1)

ATAC4 <- dba.report(ATAC3, th=0.05, bUsePval=FALSE,  bNormalized = TRUE, bCounts = TRUE)
ATAC4
write.table(ATAC4,"ATAC_DB_Frozen_0624.txt", sep="\t",row.names=FALSE)
dba.plotPCA(ATAC3, contrast = 1,attributes=DBA_FACTOR, vColors = c("blue", "orange", "red", "grey"), label=DBA_CONTROL)

dba.plotMA(ATAC3)

crukBlue.save <- getFromNamespace("crukBlue", "DiffBind")
crukMagenta.save <- getFromNamespace("crukMagenta", "DiffBind")
assignInNamespace("crukBlue","black","DiffBind")
assignInNamespace("crukMagenta","orange","DiffBind")
dba.plotMA(ATAC3)

dba.plotVolcano(ATAC3)

sum(ATAC.DB$Fold<0)
sum(ATAC.DB$Fold>0)

pvals <- dba.plotBox(ATAC3)

# correlation heatmaps of frozen samples
corvals <- dba.plotHeatmap(ATAC3, contrast=1, ColAttributes = DBA_FACTOR, colSideCols = c("mediumpurple1", "palegreen2", "lightskyblue", "orange", "orchid1"), colScheme="Blues")

# Open chromatin regions heatmaps of frozen samples
corvals <- dba.plotHeatmap(ATAC3, contrast=1, correlation = FALSE, scale="row", ColAttributes = DBA_FACTOR, colSideCols = c("mediumpurple1", "palegreen2", "lightskyblue", "orange", "orchid1"), colScheme="Blues")

#Start with Granges for ATAC3
counts <- dba.peakset(ATAC3, bRetrieve=TRUE)
?dba.peakset
counts
counts_gene <- annotatePeakInBatch(counts, AnnotationData=annoData)
counts_gene

#Start with Granges for ATAC.DB of frozen samples

# Give ranges numeric names in order
names(counts_gene) <- c(1:length(counts_gene))
counts_gene
write.table(counts_gene, "counts_gene_frozen_06242020.txt", sep = "\t", col.names = T, row.names = F, quote = F)


merge_file <- merge(ATAC_DB_Frozen_0624, counts_gene_frozen_06242020, by =c("seqnames", "start", "end"))

write.table(merge_file, "merge_table_06242020_frozen.txt", sep = "\t", col.names = T, row.names = F, quote = F)







