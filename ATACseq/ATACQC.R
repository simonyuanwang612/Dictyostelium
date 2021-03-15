library(BiocManager)
library(ATACseqQC)
library(ChIPpeakAnno)
library(MotifDb)
library(GenomicAlignments)
library(Rsamtools)
library(BSgenome.Ddiscoideum.ensembl.27)

seqlev <- "DDB0232428"
Ddiscoideum
which <- as(seqinfo(Ddiscoideum)[seqlev], "GRanges")
which
gal <- readBamFile(bamFile = "ATAC_bam_files_2nd/I4.trim.sort.bam", tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
objs <- shiftGAlignmentsList(gal)
shiftedBamfile <- file.path(outPath, "shifted.bam")
shiftedBamfile
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)

#Produce fragment length and read density files
V1_second <- c("ATAC_bam_files_2nd/G1.bam")
bamfile.labels <- gsub(".bam", "", basename(V1_second))
fragSize <- fragSizeDist(V1_second, bamfile.labels)

V2_second <- c("ATAC_bam_files_2nd/G2.bam")
bamfile.labels <- gsub(".bam", "", basename(V2_second))
fragSize <- fragSizeDist(V2_second, bamfile.labels)

V3_second <- c("ATAC_bam_files_2nd/G3.bam")
bamfile.labels <- gsub(".bam", "", basename(V3_second))
fragSize <- fragSizeDist(V3_second, bamfile.labels)

St1_second <- c("ATAC_bam_files_2nd/H1.bam")
bamfile.labels <- gsub(".bam", "", basename(St1_second))
fragSize <- fragSizeDist(St1_second, bamfile.labels)

St2_second <- c("ATAC_bam_files_2nd/H2.bam")
bamfile.labels <- gsub(".bam", "", basename(St2_second))
fragSize <- fragSizeDist(St2_second, bamfile.labels)

St3_second <- c("ATAC_bam_files_2nd/H3.bam")
bamfile.labels <- gsub(".bam", "", basename(St3_second))
fragSize <- fragSizeDist(St3_second, bamfile.labels)

M1_second <- c("ATAC_bam_files_2nd/I1.bam")
bamfile.labels <- gsub(".bam", "", basename(M1_second))
fragSize <- fragSizeDist(M1_second, bamfile.labels)

M2_second <- c("ATAC_bam_files_2nd/I2.bam")
bamfile.labels <- gsub(".bam", "", basename(M2_second))
fragSize <- fragSizeDist(M2_second, bamfile.labels)

M3_second <- c("ATAC_bam_files_2nd/I3.bam")
bamfile.labels <- gsub(".bam", "", basename(M3_second))
fragSize <- fragSizeDist(M3_second, bamfile.labels)

M4_second <- c("ATAC_bam_files_2nd/I4.bam")
bamfile.labels <- gsub(".bam", "", basename(M4_second))
fragSize <- fragSizeDist(M4_second, bamfile.labels)

F1_second <- c("ATAC_bam_files_2nd/J1.bam")
bamfile.labels <- gsub(".bam", "", basename(F1_second))
fragSize <- fragSizeDist(F1_second, bamfile.labels)

F2_second <- c("ATAC_bam_files_2nd/J2.bam")
bamfile.labels <- gsub(".bam", "", basename(F2_second))
fragSize <- fragSizeDist(F2_second, bamfile.labels)


F3_second <- c("ATAC_bam_files_2nd/J3.bam")
bamfile.labels <- gsub(".bam", "", basename(F3_second))
fragSize <- fragSizeDist(F3_second, bamfile.labels)


F4_second <- c("ATAC_bam_files_2nd/J4.bam")
bamfile.labels <- gsub(".bam", "", basename(F4_second))
fragSize <- fragSizeDist(F4_second, bamfile.labels)

