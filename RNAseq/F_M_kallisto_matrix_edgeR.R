counts_strand_FM <- read.delim("~/Desktop/Dicty_RNA_seq_kallisto/kallisto_dicty_matrix.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
x_FM <-  counts_strand_FM[c( "lib17", "lib18", "lib19", "lib20")]
head(x_FM)
tail(x_FM)
dim(x_FM)

#create groups
group_FM <- c(rep("M", 2), rep("F", 2))
group_FM

#Building the edgeR Object
#create DGEList object to work on, where are the columns are the counts
#Put the counts and other information into a DGEList object
y_FM <- DGEList(counts=x_FM[,1:4], group=group_FM)
y_FM
dim(y_FM)

# Human introduced error
#Filter out low count reads since it would be impossible to detect differential expression. The method used in the
#edgeR vignette is to keep only those genes that have at least 1 read per million in at least 2 samples. Therefore, Filter out very lowly expressed tags, 
#keeping genes that are expressed at a reasonable level in at least one condition. Since the group sizes are 2, keep genes that achieve at least 
#one count per million (cpm) in 2 samples
#Once this filtering is done, calculate the normalization factors which correct for the different compositions of the samples. The effective
#library sizes are then the product of the actual library sizes and these factors.

keep_FM <- rowSums(cpm(y_FM)>1) >= 2
y_FM <- y_FM[keep_FM,]
dim(y_FM)
#Re-compute the library sizes
y_FM$samples$lib.size <- colSums(y_FM$counts)
y_FM$samples

##TMM Normalization
#Compute effective library sizes using Trimmed Mean of M Values (TMM) normalization
#This accounts for differences in transcriptome sizes
y_FM <- calcNormFactors(y_FM)
y_FM$samples

#effective library sizes
# This is relative normalization
y_FM$samples$lib.size*y_FM$samples$norm.factors

##Data exploration
#Generate a Multi-dimensional Scaling (MDS) plot
#An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions.
#An MDS plot show distances, in terms of biological coefficient of variation (BCV), between samples
#Dimension 1 separates the Normal from TNBC samples. The MDS plot shows the replicates are consistent. 
#This is a good indication edgeR analysis will produce adequate DE genes.
#To view MDS plot
plotMDS(y_FM, main ="MDS Plot for TMM Normalized Count Data F_M")

#Write plot as a pdf
pdf("MDS_Plot_F_M.pdf" , width = 11 , height = 7) # in inches
plotMDS(y_FM, main ="MDS Plot for TMM Normalized Count Data_FM")
dev.off() # this tells [R] to close and stop writing to the pdf

#Generate the heatmap
install.packages("mixOmics")
library(mixOmics)
sampleDists <- as.matrix(dist(t(y$pseudo.counts)))
sampleDists
install.packages("RColorBrewer")
library("RColorBrewer") 
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cim(sampleDists, col=cimColor, symkey=FALSE)

##Estimating common and tagwise dispersion
#The common dispersion estimates the overall BCV of the dataset, averaged over all genes with the qCML method
y_FM <- estimateCommonDisp(y_FM, verbose=TRUE)
# Disp = 0.0138 , BCV = 0.117 

#Now estimate gene-specific (tagwise) dispersions with the qCML method
y_FM <- estimateTagwiseDisp(y_FM)
summary(y_FM$tagwise.dispersion)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0002  0.0002  0.0002  0.0049  0.0071  0.2938

#Plot the estimated dispersions 
#In general, BCV is higher with lower average log CPM
#To view estimated dispersions plot
plotBCV(y_FM, main = "Plot of Estimated Dispersions_F_M")


#Write plot as a pdf
pdf("Est_Dispersions_Plot_F_M.pdf" , width = 11 , height = 7) # in inches
plotBCV(y_FM, main = "Plot of Estimated Dispersions_F_M")
dev.off() # this tells [R] to close and stop writing to the pdf

##Differential expression
#FDR is less significant than the original P value, correcting for the P value. Set false discovery treshold.
#Compute NB exact genewise tests for differential expression between samples
et_FM <- exactTest(y_FM, pair=2:1)
options(digits = 3) #print only 3 digits
top_FM <- topTags(et_FM)
top_FM

#Check the individual cpm values for the top genes
#Check up/down expression to  *VERY IMPORTANT*
cpm(y_FM)[rownames(top_FM), ]

#print the total number of DE genes at 5% FDR
summary(de_FM <- decideTestsDGE(et_FM))

#Plot the log-fold-changes with a smear plot, highlighting the DE genes
#Blue lines indicate 2-fold changes
#to view Smear Plot
detags_FM <- rownames(y_FM)[as.logical(de_FM)]
plotSmear(et_FM, de.tags=detags_FM, main = "Smear Plot F-M")

abline(h=c(-1, 1), col="blue")

#Write plot as a pdf
pdf("Smear_Plot_F_M.pdf" , width = 7 , height = 7) # in inches
plotSmear(et_FM, de.tags=detags_FM, main = "Smear Plot F-M")
abline(h=c(-1, 1), col="blue")
dev.off() # this tells [R] to close and stop writing to the pdf

?plotSmear
#Create an object called "resultsTbl.tagwise" containing the four output parameters: logFC, logCPM, PValue,FDR
resultsTbl.tagwise_FM <- topTags(et_FM, n = nrow(et_FM$table))$table


#Check resultsTbl.tagwise
head(resultsTbl.tagwise_FM)
dim(resultsTbl.tagwise_FM)

#Generate objects according to three significance levels for tagwise dispersion estimates: 0.05, 0.001, 1e-06
#These will be used for GOseq analyses in classes 5 and 6.
#Names/IDs of DE genes
de.genes.tagwise_0.05_FM <- rownames(resultsTbl.tagwise_FM)[resultsTbl.tagwise_FM$FDR <= 0.05]
de.genes.tagwise_0.001_FM <- rownames(resultsTbl.tagwise_FM)[resultsTbl.tagwise_FM$FDR <= 0.001]
de.genes.tagwise_1e_06_FM <- rownames(resultsTbl.tagwise_FM)[resultsTbl.tagwise_FM$FDR <= 1e-06]

#Check number of genes at each significance level
length(de.genes.tagwise_0.05_FM) #2174 DE genes
length(de.genes.tagwise_0.001_FM) #1245 DE genes
length(de.genes.tagwise_1e_06_FM) #736 DE genes

#Check percentage of total genes
length(de.genes.tagwise_0.05_FM)/nrow(resultsTbl.tagwise_FM)*100 #22
length(de.genes.tagwise_0.001_FM)/nrow(resultsTbl.tagwise_FM)*100 #12.6
length(de.genes.tagwise_1e_06_FM)/nrow(resultsTbl.tagwise_FM)*100 #7.44

#Check Up/Down regulated summaries for tagwise results
summary(de_FM <- decideTestsDGE(et_FM, p.value=0.05))
#Down    854
#NotSig 7927
#Up      1320
summary(de_FM <- decideTestsDGE(et_FM, p.value=0.001))
#Down    440
#NotSig 8856
#Up      805
summary(de_FM <- decideTestsDGE(et_FM, p.value=1e-06))
#Down    201
#NotSig 9365
#Up      535
head(de.genes.tagwise_1e_06_FM)

#Double check
is.vector(de.genes.tagwise_1e_06_FM)

##Format Results and count data form DGEList object(y)
colnames(resultsTbl.tagwise_FM) <- c("logFC", "logCPM","PValue", "FDR")
wh.rows.tagwise_FM <- match(rownames(resultsTbl.tagwise_FM), rownames(y_FM$counts))


#Check for all 10101 genes
length(wh.rows.tagwise_FM)

#Combine exact test results, count data from DGEList object(y) and tagwise dispersion
combResults.tagwise_FM <- cbind(resultsTbl.tagwise_FM,
                                "Tgw.Disp" = y_FM$tagwise.dispersion[wh.rows.tagwise_FM],
                                "UpDown.Tgw" = decideTestsDGE(et_FM, p.value = 0.05)[wh.rows.tagwise_FM],
                                y_FM$counts[wh.rows.tagwise_FM,])

#Check combResults.tagwise 
head(combResults.tagwise_FM)
dim(combResults.tagwise_FM)

#Create significance level objects from original counts data (x) and FDR vectors
de.tag_05_FM<-x_FM[which(row.names(x_FM)%in%(de.genes.tagwise_0.05_FM)),]
de.tag_001_FM<-x_FM[which(row.names(x_FM)%in%(de.genes.tagwise_0.001_FM)),]
de.tag_1e_06_FM<-x_FM[which(row.names(x_FM)%in%(de.genes.tagwise_1e_06_FM)),]

#Check
head(de.tag_05_FM)
dim(de.tag_05_FM)
head(de.tag_001_FM)
dim(de.tag_001_FM)
head(de.tag_1e_06_FM)
dim(de.tag_1e_06_FM)

#Check again
head(combResults.tagwise_FM)
head(resultsTbl.tagwise_FM)

dim(combResults.tagwise_FM)

##Merge all output into FDR specific objects
combResults.tagwise_05_FM<-merge(combResults.tagwise_FM [,1:6], de.tag_05_FM, by="row.names")
head(combResults.tagwise_05_FM)

#Rename column 1 as "target_id"
colnames(combResults.tagwise_05_FM)[1] ="target_id"
head(combResults.tagwise_05_FM)
dim(combResults.tagwise_05_FM)

#Create edgeR_05 output object and reorder columns 
edgeR_05_FM=(combResults.tagwise_05_FM)
head(edgeR_05_FM)
dim(edgeR_05_FM)

#Repeat for FDR of 0.001
combResults.tagwise_001_FM<-merge(combResults.tagwise_FM [,1:6], de.tag_001_FM, by="row.names")
head(combResults.tagwise_001_FM)

#Rename column 1 as "Gene"
colnames(combResults.tagwise_001_FM)[1] ="target_id"
head(combResults.tagwise_001_FM)
dim(combResults.tagwise_001_FM)

#Create edgeR_001 output object and reorder columns
edgeR_001_FM=combResults.tagwise_001_FM
head(edgeR_001_FM, 20)
dim(edgeR_001_FM)

#Repeat for FDR of 1e06
combResults.tagwise_1e_06_FM<-merge(combResults.tagwise_FM [,1:6], de.tag_1e_06_FM, by="row.names")
head(combResults.tagwise_1e_06_FM)

#Rename column 1 as "target_id"
colnames(combResults.tagwise_1e_06_FM)[1] ="target_id"
head(combResults.tagwise_1e_06_FM)
dim(combResults.tagwise_1e_06_FM)


#Create edgeR_1e_06 output object and reorder columns
edgeR_1e_06_FM = combResults.tagwise_1e_06_FM
head(edgeR_1e_06_FM)
dim(edgeR_1e_06_FM)

#Sort results by FDR as a check for row names.
edgeR_05_out_FM<-edgeR_05_FM[order(edgeR_05_FM[,7],decreasing=FALSE),]
head(edgeR_05_out_FM)

edgeR_001_out_FM<-edgeR_001_FM[order(edgeR_001_FM[,7],decreasing=FALSE),]
head(edgeR_001_out_FM)

edgeR_1e_06_out_FM<-edgeR_1e_06_FM[order(edgeR_1e_06_FM[,7],decreasing=FALSE),]
head(edgeR_1e_06_out_FM)



edgeR_05_out_geneid_FM <- merge(edgeR_05_out_FM, mart_export, by = "target_id")
edgeR_001_out_geneid_FM <- merge(edgeR_001_out_FM, mart_export, by = "target_id")
edgeR_1e_06_out_geneid_FM <- merge(edgeR_1e_06_out_FM, mart_export, by = "target_id")

write.table(edgeR_05_out_geneid_FM, file="edgeR_tagwise_05_output_F_M_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_001_out_geneid_FM, file="edgeR_tagwise_001_output_F_M_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_1e_06_out_geneid_FM, file="edgeR_tagwise_1e_06_output_F_M_kali.txt", row.names=FALSE, sep="\t")

dim(edgeR_05_out_geneid_FM)
dim(edgeR_001_out_geneid_FM)
dim(edgeR_1e_06_out_geneid_FM)
