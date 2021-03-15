new <- merge(mart_export, kallisto_dicty_matrix, by = "target_id")
head(new)
dim(new)
write.table(new, "kallisto_dicty_matrix_geneID.txt", sep="\t")
head(new)
library(edgeR)
#read in entire counts data .txt file
setwd("/Users/xm/Desktop/Dicty_RNA_seq_kallisto/")
ls()
dir()

#counts_strand <- read.delim("~/Desktop/edgeR/counts_strand.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
counts_strand_SV <- read.delim("~/Desktop/Dicty_RNA_seq_kallisto/kallisto_dicty_matrix.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
x <-  counts_strand_SV[c( "lib13", "lib14", "lib15", "lib16")]
head(x)
tail(x)
dim(x)

#create groups
group <- c(rep("V", 2), rep("S", 2))
group

#Building the edgeR Object
#create DGEList object to work on, where are the columns are the counts
#Put the counts and other information into a DGEList object
y <- DGEList(counts=x[,1:4], group=group)
y
dim(y)

# Human introduced error
#Filter out low count reads since it would be impossible to detect differential expression. The method used in the
#edgeR vignette is to keep only those genes that have at least 1 read per million in at least 2 samples. Therefore, Filter out very lowly expressed tags, 
#keeping genes that are expressed at a reasonable level in at least one condition. Since the group sizes are 2, keep genes that achieve at least 
#one count per million (cpm) in 2 samples
#Once this filtering is done, calculate the normalization factors which correct for the different compositions of the samples. The effective
#library sizes are then the product of the actual library sizes and these factors.

keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep,]
dim(y)
#Re-compute the library sizes
y$samples$lib.size <- colSums(y$counts)
y$samples

##TMM Normalization
#Compute effective library sizes using Trimmed Mean of M Values (TMM) normalization
#This accounts for differences in transcriptome sizes
y <- calcNormFactors(y)
y$samples

#effective library sizes
# This is relative normalization
y$samples$lib.size*y$samples$norm.factors

##Data exploration
#Generate a Multi-dimensional Scaling (MDS) plot
#An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions.
#An MDS plot show distances, in terms of biological coefficient of variation (BCV), between samples
#Dimension 1 separates the Normal from TNBC samples. The MDS plot shows the replicates are consistent. 
#This is a good indication edgeR analysis will produce adequate DE genes.
#To view MDS plot
plotMDS(y, main ="MDS Plot for TMM Normalized Count Data V_S")

#Write plot as a pdf
pdf("MDS_Plot_V_S.pdf" , width = 11 , height = 7) # in inches
plotMDS(y, main ="MDS Plot for TMM Normalized Count Data")
dev.off() # this tells [R] to close and stop writing to the pdf

#Generate the heatmap
install.packages("mixOmics")
library(mixOmics)
sampleDists <- as.matrix(dist(t(y$pseudo.counts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)

##Estimating common and tagwise dispersion
#The common dispersion estimates the overall BCV of the dataset, averaged over all genes with the qCML method
y <- estimateCommonDisp(y, verbose=TRUE)
# Disp = 0.01977, BCV = 0.1406

#Now estimate gene-specific (tagwise) dispersions with the qCML method
y <- estimateTagwiseDisp(y)
summary(y$tagwise.dispersion)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0003089 0.0003089 0.0003089 0.0064673 0.0054748 0.3500564

#Plot the estimated dispersions 
#In general, BCV is higher with lower average log CPM
#To view estimated dispersions plot
plotBCV(y, main = "Plot of Estimated Dispersions_S_V")

#Plot the estimated dispersions 
#In general, BCV is higher with lower average log CPM
#To view estimated dispersions plot
plotBCV(y, main = "Plot of Estimated Dispersions_S_V")

#Write plot as a pdf
pdf("Est_Dispersions_Plot_S_V.pdf" , width = 11 , height = 7) # in inches
plotBCV(y, main = "Plot of Estimated Dispersions_S_V")
dev.off() # this tells [R] to close and stop writing to the pdf

##Differential expression
#FDR is less significant than the original P value, correcting for the P value. Set false discovery treshold.
#Compute NB exact genewise tests for differential expression between samples
et <- exactTest(y, pair=2:1)
options(digits = 3) #print only 3 digits
top <- topTags(et)
top

#Check the individual cpm values for the top genes
#Check up/down expression to  *VERY IMPORTANT*
cpm(y)[rownames(top), ]

#print the total number of DE genes at 5% FDR
summary(de <- decideTestsDGE(et))

#Plot the log-fold-changes with a smear plot, highlighting the DE genes
#Blue lines indicate 2-fold changes
#to view Smear Plot
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags, main = "Smear Plot S-V")

abline(h=c(-1, 1), col="blue")

#Write plot as a pdf
pdf("Smear_Plot_S_V.pdf" , width = 7 , height = 7) # in inches
plotSmear(et, de.tags=detags, main = "Smear Plot S-V")
abline(h=c(-1, 1), col="blue")
dev.off() # this tells [R] to close and stop writing to the pdf

?plotSmear
#Create an object called "resultsTbl.tagwise" containing the four output parameters: logFC, logCPM, PValue,FDR
resultsTbl.tagwise <- topTags(et, n = nrow(et$table))$table


#Check resultsTbl.tagwise
head(resultsTbl.tagwise)
dim(resultsTbl.tagwise)

#Generate objects according to three significance levels for tagwise dispersion estimates: 0.05, 0.001, 1e-06
#These will be used for GOseq analyses in classes 5 and 6.
#Names/IDs of DE genes
de.genes.tagwise_0.05 <- rownames(resultsTbl.tagwise)[resultsTbl.tagwise$FDR <= 0.05]
de.genes.tagwise_0.001 <- rownames(resultsTbl.tagwise)[resultsTbl.tagwise$FDR <= 0.001]
de.genes.tagwise_1e_06 <- rownames(resultsTbl.tagwise)[resultsTbl.tagwise$FDR <= 1e-06]

#Check number of genes at each significance level
length(de.genes.tagwise_0.05) #3021 DE genes
length(de.genes.tagwise_0.001) #1936 DE genes
length(de.genes.tagwise_1e_06) #1174 DE genes

#Check percentage of total genes
length(de.genes.tagwise_0.05)/nrow(resultsTbl.tagwise)*100 #17.8
length(de.genes.tagwise_0.001)/nrow(resultsTbl.tagwise)*100 #8.14
length(de.genes.tagwise_1e_06)/nrow(resultsTbl.tagwise)*100 #2.98

#Check Up/Down regulated summaries for tagwise results
summary(de <- decideTestsDGE(et, p.value=0.05))
#Down    1529
#NotSig 6872
#Up      1492
summary(de <- decideTestsDGE(et, p.value=0.001))
#Down    1028
#NotSig 7957
#Up      908
summary(de <- decideTestsDGE(et, p.value=1e-06))
#Down    674
#NotSig 8719
#Up      107
head(de.genes.tagwise_1e_06)

#Double check
is.vector(de.genes.tagwise_1e_06)

##Format Results and count data form DGEList object(y)
colnames(resultsTbl.tagwise) <- c("logFC", "logCPM","PValue", "FDR")
wh.rows.tagwise <- match(rownames(resultsTbl.tagwise), rownames(y$counts))


#Check for all 9893 genes
length(wh.rows.tagwise)

#Combine exact test results, count data from DGEList object(y) and tagwise dispersion
combResults.tagwise <- cbind(resultsTbl.tagwise,
                             "Tgw.Disp" = y$tagwise.dispersion[wh.rows.tagwise],
                             "UpDown.Tgw" = decideTestsDGE(et, p.value = 0.05)[wh.rows.tagwise],
                             y$counts[wh.rows.tagwise,])

#Check combResults.tagwise 
head(combResults.tagwise)
dim(combResults.tagwise)

#Create significance level objects from original counts data (x) and FDR vectors
de.tag_05<-x[which(row.names(x)%in%(de.genes.tagwise_0.05)),]
de.tag_001<-x[which(row.names(x)%in%(de.genes.tagwise_0.001)),]
de.tag_1e_06<-x[which(row.names(x)%in%(de.genes.tagwise_1e_06)),]

#Check
head(de.tag_05)
dim(de.tag_05)
head(de.tag_001)
dim(de.tag_001)
head(de.tag_1e_06)
dim(de.tag_1e_06)

#Check again
head(combResults.tagwise)
head(resultsTbl.tagwise)

dim(combResults.tagwise)

##Merge all output into FDR specific objects
combResults.tagwise_05<-merge(combResults.tagwise [,1:6], de.tag_05, by="row.names")
head(combResults.tagwise_05)

#Rename column 1 as "Gene"
colnames(combResults.tagwise_05)[1] ="target_id"
head(combResults.tagwise_05)
dim(combResults.tagwise_05)

#Create edgeR_05 output object and reorder columns 
edgeR_05=(combResults.tagwise_05)
head(edgeR_05)
dim(edgeR_05)

#Repeat for FDR of 0.001
combResults.tagwise_001<-merge(combResults.tagwise [,1:6], de.tag_001, by="row.names")
head(combResults.tagwise_001)

#Rename column 1 as "Gene"
colnames(combResults.tagwise_001)[1] ="target_id"
head(combResults.tagwise_001)
dim(combResults.tagwise_001)

#Create edgeR_001 output object and reorder columns
edgeR_001=combResults.tagwise_001
head(edgeR_001, 20)
dim(edgeR_001)

#Repeat for FDR of 1e06
combResults.tagwise_1e_06<-merge(combResults.tagwise [,1:6], de.tag_1e_06, by="row.names")
head(combResults.tagwise_1e_06)

#Rename column 1 as "Gene"
colnames(combResults.tagwise_1e_06)[1] ="target_id"
head(combResults.tagwise_1e_06)
dim(combResults.tagwise_1e_06)


#Create edgeR_1e_06 output object and reorder columns
edgeR_1e_06 = combResults.tagwise_1e_06
head(edgeR_1e_06)
dim(edgeR_1e_06)

#Sort results by FDR as a check for row names.
edgeR_05_out<-edgeR_05[order(edgeR_05[,7],decreasing=FALSE),]
head(edgeR_05_out)

edgeR_001_out<-edgeR_001[order(edgeR_001[,7],decreasing=FALSE),]
head(edgeR_001_out)

edgeR_1e_06_out<-edgeR_1e_06[order(edgeR_1e_06[,7],decreasing=FALSE),]
head(edgeR_1e_06_out)

#Write all three edgeR FDR specific objects to .txt files - row.names=FALSE to maintain appropriate column order 
#write.table(edgeR_05_out, file="edgeR_tagwise_05_output_S_V.txt", row.names=FALSE, sep="\t")
#write.table(edgeR_001_out, file="edgeR_tagwise_001_output_S_V.txt", row.names=FALSE, sep="\t")
#write.table(edgeR_1e_06_out, file="edgeR_tagwise_1e_06_output_S_V.txt", row.names=FALSE, sep="\t")

edgeR_05_out_geneid <- merge(edgeR_05_out, mart_export, by = "target_id")
edgeR_001_out_geneid <- merge(edgeR_001_out, mart_export, by = "target_id")
edgeR_1e_06_out_geneid <- merge(edgeR_1e_06_out, mart_export, by = "target_id")

write.table(edgeR_05_out_geneid, file="edgeR_tagwise_05_output_S_V_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_001_out_geneid, file="edgeR_tagwise_001_output_S_V_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_1e_06_out_geneid, file="edgeR_tagwise_1e_06_output_S_V_kali.txt", row.names=FALSE, sep="\t")


dim(edgeR_05_out_geneid)
dim(edgeR_001_out_geneid)
dim(edgeR_1e_06_out_geneid)
