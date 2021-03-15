#counts_strand <- read.delim("~/Desktop/edgeR/counts_strand.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
counts_strand_MS <- read.delim("~/Desktop/Dicty_RNA_seq_kallisto/kallisto_dicty_matrix.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
x_MS <-  counts_strand_MS[c( "lib15", "lib16", "lib17", "lib18")]
head(x_MS)
tail(x_MS)
dim(x_MS)

#create groups
group_MS <- c(rep("S", 2), rep("M", 2))
group_MS

#Building the edgeR Object
#create DGEList object to work on, where are the columns are the counts
#Put the counts and other information into a DGEList object
y_MS <- DGEList(counts=x_MS[,1:4], group=group_MS)
y_MS
dim(y_MS)

# Human introduced error
#Filter out low count reads since it would be impossible to detect differential expression. The method used in the
#edgeR vignette is to keep only those genes that have at least 1 read per million in at least 2 samples. Therefore, Filter out very lowly expressed tags, 
#keeping genes that are expressed at a reasonable level in at least one condition. Since the group sizes are 2, keep genes that achieve at least 
#one count per million (cpm) in 2 samples
#Once this filtering is done, calculate the normalization factors which correct for the different compositions of the samples. The effective
#library sizes are then the product of the actual library sizes and these factors.

keep_MS <- rowSums(cpm(y_MS)>1) >= 2
y_MS <- y_MS[keep_MS,]
dim(y_MS)
#Re-compute the library sizes
y_MS$samples$lib.size <- colSums(y_MS$counts)
y_MS$samples

##TMM Normalization
#Compute effective library sizes using Trimmed Mean of M Values (TMM) normalization
#This accounts for differences in transcriptome sizes
y_MS <- calcNormFactors(y_MS)
y_MS$samples

#effective library sizes
# This is relative normalization
y_MS$samples$lib.size*y_MS$samples$norm.factors

##Data exploration
#Generate a Multi-dimensional Scaling (MDS) plot
#An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions.
#An MDS plot show distances, in terms of biological coefficient of variation (BCV), between samples
#Dimension 1 separates the Normal from TNBC samples. The MDS plot shows the replicates are consistent. 
#This is a good indication edgeR analysis will produce adequate DE genes.
#To view MDS plot
plotMDS(y_MS, main ="MDS Plot for TMM Normalized Count Data M_S")

#Write plot as a pdf
pdf("MDS_Plot_M_S.pdf" , width = 11 , height = 7) # in inches
plotMDS(y_MS, main ="MDS Plot for TMM Normalized Count Data_MS")
dev.off() # this tells [R] to close and stop writing to the pdf

#Generate the heatmap
install.packages("mixOmics")
library(mixOmics)
sampleDists <- as.matrix(dist(t(y$pseudo.counts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)

##Estimating common and tagwise dispersion
#The common dispersion estimates the overall BCV of the dataset, averaged over all genes with the qCML method
y_MS <- estimateCommonDisp(y_MS, verbose=TRUE)
# Disp = 0.0138 , BCV = 0.117 

#Now estimate gene-specific (tagwise) dispersions with the qCML method
y_MS <- estimateTagwiseDisp(y_MS)
summary(y_MS$tagwise.dispersion)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0002  0.0002  0.0002  0.0049  0.0071  0.2938

#Plot the estimated dispersions 
#In general, BCV is higher with lower average log CPM
#To view estimated dispersions plot
plotBCV(y_MS, main = "Plot of Estimated Dispersions_M_S")


#Write plot as a pdf
pdf("Est_Dispersions_Plot_M_S.pdf" , width = 11 , height = 7) # in inches
plotBCV(y_MS, main = "Plot of Estimated Dispersions_M_S")
dev.off() # this tells [R] to close and stop writing to the pdf

##Differential expression
#FDR is less significant than the original P value, correcting for the P value. Set false discovery treshold.
#Compute NB exact genewise tests for differential expression between samples
et_MS <- exactTest(y_MS, pair=2:1)
options(digits = 3) #print only 3 digits
top_MS <- topTags(et_MS)
top_MS

#Check the individual cpm values for the top genes
#Check up/down expression to  *VERY IMPORTANT*
cpm(y_MS)[rownames(top_MS), ]

#print the total number of DE genes at 5% FDR
summary(de_MS <- decideTestsDGE(et_MS))

#Plot the log-fold-changes with a smear plot, highlighting the DE genes
#Blue lines indicate 2-fold changes
#to view Smear Plot
detags_MS <- rownames(y_MS)[as.logical(de_MS)]
plotSmear(et_MS, de.tags=detags_MS, main = "Smear Plot M-S")

abline(h=c(-1, 1), col="blue")

#Write plot as a pdf
pdf("Smear_Plot_M_S.pdf" , width = 7 , height = 7) # in inches
plotSmear(et_MS, de.tags=detags_MS, main = "Smear Plot M-S")
abline(h=c(-1, 1), col="blue")
dev.off() # this tells [R] to close and stop writing to the pdf

?plotSmear
#Create an object called "resultsTbl.tagwise" containing the four output parameters: logFC, logCPM, PValue,FDR
resultsTbl.tagwise_MS <- topTags(et_MS, n = nrow(et_MS$table))$table


#Check resultsTbl.tagwise
head(resultsTbl.tagwise_MS)
dim(resultsTbl.tagwise_MS)

#Generate objects according to three significance levels for tagwise dispersion estimates: 0.05, 0.001, 1e-06
#These will be used for GOseq analyses in classes 5 and 6.
#Names/IDs of DE genes
de.genes.tagwise_0.05_MS <- rownames(resultsTbl.tagwise_MS)[resultsTbl.tagwise_MS$FDR <= 0.05]
de.genes.tagwise_0.001_MS <- rownames(resultsTbl.tagwise_MS)[resultsTbl.tagwise_MS$FDR <= 0.001]
de.genes.tagwise_1e_06_MS <- rownames(resultsTbl.tagwise_MS)[resultsTbl.tagwise_MS$FDR <= 1e-06]

#Check number of genes at each significance level
length(de.genes.tagwise_0.05_MS) #2174 DE genes
length(de.genes.tagwise_0.001_MS) #1245 DE genes
length(de.genes.tagwise_1e_06_MS) #736 DE genes

#Check percentage of total genes
length(de.genes.tagwise_0.05_MS)/nrow(resultsTbl.tagwise_MS)*100 #21.5
length(de.genes.tagwise_0.001_MS)/nrow(resultsTbl.tagwise_MS)*100 #12.3
length(de.genes.tagwise_1e_06_MS)/nrow(resultsTbl.tagwise_MS)*100 #7.29

#Check Up/Down regulated summaries for tagwise results
summary(de_MS <- decideTestsDGE(et_MS, p.value=0.05))
#Down    854
#NotSig 7927
#Up      1320
summary(de_MS <- decideTestsDGE(et_MS, p.value=0.001))
#Down    440
#NotSig 8856
#Up      805
summary(de_MS <- decideTestsDGE(et_MS, p.value=1e-06))
#Down    201
#NotSig 9365
#Up      535
head(de.genes.tagwise_1e_06_MS)

#Double check
is.vector(de.genes.tagwise_1e_06_MS)

##Format Results and count data form DGEList object(y)
colnames(resultsTbl.tagwise_MS) <- c("logFC", "logCPM","PValue", "FDR")
wh.rows.tagwise_MS <- match(rownames(resultsTbl.tagwise_MS), rownames(y_MS$counts))


#Check for all 10101 genes
length(wh.rows.tagwise_MS)

#Combine exact test results, count data from DGEList object(y) and tagwise dispersion
combResults.tagwise_MS <- cbind(resultsTbl.tagwise_MS,
                             "Tgw.Disp" = y_MS$tagwise.dispersion[wh.rows.tagwise_MS],
                             "UpDown.Tgw" = decideTestsDGE(et_MS, p.value = 0.05)[wh.rows.tagwise_MS],
                             y_MS$counts[wh.rows.tagwise_MS,])

#Check combResults.tagwise 
head(combResults.tagwise_MS)
dim(combResults.tagwise_MS)

#Create significance level objects from original counts data (x) and FDR vectors
de.tag_05_MS<-x_MS[which(row.names(x_MS)%in%(de.genes.tagwise_0.05_MS)),]
de.tag_001_MS<-x_MS[which(row.names(x_MS)%in%(de.genes.tagwise_0.001_MS)),]
de.tag_1e_06_MS<-x_MS[which(row.names(x_MS)%in%(de.genes.tagwise_1e_06_MS)),]

#Check
head(de.tag_05_MS)
dim(de.tag_05_MS)
head(de.tag_001_MS)
dim(de.tag_001_MS)
head(de.tag_1e_06_MS)
dim(de.tag_1e_06_MS)

#Check again
head(combResults.tagwise_MS)
head(resultsTbl.tagwise_MS)

dim(combResults.tagwise_MS)

##Merge all output into FDR specific objects
combResults.tagwise_05_MS<-merge(combResults.tagwise_MS [,1:6], de.tag_05_MS, by="row.names")
head(combResults.tagwise_05_MS)

#Rename column 1 as "target_id"
colnames(combResults.tagwise_05_MS)[1] ="target_id"
head(combResults.tagwise_05_MS)
dim(combResults.tagwise_05_MS)

#Create edgeR_05 output object and reorder columns 
edgeR_05_MS=(combResults.tagwise_05_MS)
head(edgeR_05_MS)
dim(edgeR_05_MS)

#Repeat for FDR of 0.001
combResults.tagwise_001_MS<-merge(combResults.tagwise_MS [,1:6], de.tag_001_MS, by="row.names")
head(combResults.tagwise_001_MS)

#Rename column 1 as "Gene"
colnames(combResults.tagwise_001_MS)[1] ="target_id"
head(combResults.tagwise_001_MS)
dim(combResults.tagwise_001_MS)

#Create edgeR_001 output object and reorder columns
edgeR_001_MS=combResults.tagwise_001_MS
head(edgeR_001_MS, 20)
dim(edgeR_001_MS)

#Repeat for FDR of 1e06
combResults.tagwise_1e_06_MS<-merge(combResults.tagwise_MS [,1:6], de.tag_1e_06_MS, by="row.names")
head(combResults.tagwise_1e_06_MS)

#Rename column 1 as "target_id"
colnames(combResults.tagwise_1e_06_MS)[1] ="target_id"
head(combResults.tagwise_1e_06_MS)
dim(combResults.tagwise_1e_06_MS)


#Create edgeR_1e_06 output object and reorder columns
edgeR_1e_06_MS = combResults.tagwise_1e_06_MS
head(edgeR_1e_06_MS)
dim(edgeR_1e_06_MS)

#Sort results by FDR as a check for row names.
edgeR_05_out_MS<-edgeR_05_MS[order(edgeR_05_MS[,7],decreasing=FALSE),]
head(edgeR_05_out_MS)

edgeR_001_out_MS<-edgeR_001_MS[order(edgeR_001_MS[,7],decreasing=FALSE),]
head(edgeR_001_out_MS)

edgeR_1e_06_out_MS<-edgeR_1e_06_MS[order(edgeR_1e_06_MS[,7],decreasing=FALSE),]
head(edgeR_1e_06_out_MS)

#Write all three edgeR FDR specific objects to .txt files - row.names=FALSE to maintain appropriate column order 
#write.table(edgeR_05_out, file="edgeR_tagwise_05_output_S_V.txt", row.names=FALSE, sep="\t")
#write.table(edgeR_001_out, file="edgeR_tagwise_001_output_S_V.txt", row.names=FALSE, sep="\t")
#write.table(edgeR_1e_06_out, file="edgeR_tagwise_1e_06_output_S_V.txt", row.names=FALSE, sep="\t")

edgeR_05_out_geneid_MS <- merge(edgeR_05_out_MS, mart_export, by = "target_id")
edgeR_001_out_geneid_MS <- merge(edgeR_001_out_MS, mart_export, by = "target_id")
edgeR_1e_06_out_geneid_MS <- merge(edgeR_1e_06_out_MS, mart_export, by = "target_id")

write.table(edgeR_05_out_geneid_MS, file="edgeR_tagwise_05_output_M_S_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_001_out_geneid_MS, file="edgeR_tagwise_001_output_M_S_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_1e_06_out_geneid_MS, file="edgeR_tagwise_1e_06_output_M_S_kali.txt", row.names=FALSE, sep="\t")

dim(edgeR_05_out_geneid_MS)
dim(edgeR_001_out_geneid_MS)
dim(edgeR_1e_06_out_geneid_MS)
